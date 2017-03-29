/*
Nexus

Copyright(C) 2012 - Federico Ponchio
ISTI - Italian National Research Council - Visual Computing Lab

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
for more details.
*/

BitStream = function(array) {
	var t = this;
	t.a = array;
	t.current = array[0];
	t.position = 0; //position in the buffer
	t.pending = 32;  //bits still to read
	t.mask = new Uint32Array([0x00, 0x01, 0x03, 0x07, 0x0f, 0x01f, 0x03f, 0x07f, 0xff, 
		0x01ff, 0x03ff, 0x07ff, 0x0fff, 0x01fff, 0x03fff, 0x07fff, 0xffff, 0x01ffff,
		0x03ffff, 0x07ffff, 0x0fffff, 0x01fffff, 0x03fffff, 0x07fffff, 0xffffff, 0x01ffffff,
		0x03ffffff,  0x07ffffff,  0x0fffffff, 0x01fffffff, 0x03fffffff, 0x7fffffff]);
};
 
BitStream.prototype = { 
	read: function(bits) {
		var t = this;
		var result = 0;
		if(bits > t.pending) {
			result |= (t.current << (bits - t.pending))>>>0;
			bits -= t.pending;
			t.current = t.a[++t.position];
			t.pending = 32;
		}
		t.pending -= bits;
		result |= (t.current >>> t.pending);
		t.current = (t.current & t.mask[t.pending])>>>0;
		return result;
	}
};

Stream = function(buffer) {
	var t = this;
	t.data = buffer;
	t.buffer = new Uint8Array(buffer);
	t.pos = 0;
	t.view = new DataView(buffer);
}

Stream.prototype = {
	readChar: function() {
		var c = this.buffer[this.pos++];
		if(c > 127) c -= 256;
		return c;
	},
	readUChar: function() {
		return this.buffer[this.pos++];
	},	
	readShort: function() {
		this.pos += 2;
		return this.view.getInt16(this.pos-2, true);
	},
	readFloat: function() {
		this.pos += 4;
		return this.view.getFloat32(this.pos-4, true);
	},
	readInt: function() {
		this.pos += 4;
		return this.view.getInt32(this.pos-4, true);
	},
	readArray: function(n) {
		var a = this.buffer.subarray(this.pos, this.pos+n);
		this.pos += n;
		return a;
	},
	readString: function() {
		var n = this.readShort();
		var s = String.fromCharCode.apply(null, this.readArray(n-1));
		this.pos++; //null terminator of string.
		return s;
	},
	readBitStream:function() {
		var n = this.readInt();
		var pad = this.pos & 0x3;
		if(pad != 0)
			this.pos += 4 - pad;
		var b = new BitStream(new Uint32Array(this.data, this.pos, n));
		this.pos += n*4;
		return b;
	},
	//make decodearray2,3 later //TODO faster to create values here or passing them?
	decodeArray: function(N, values) {
		var t = this;
		var bitstream = t.readBitStream();

		var tunstall = new Tunstall;
		var logs = tunstall.decompress(this);

		for(var i =0; i < logs.length; i++) {
			var diff = logs[i];
			if(diff == 0) {
				for(var c = 0; c < N; c++)
					values[i*N + c] = 0;
				continue;
			}
			//making a single read is 2/3 faster
			//uint64_t &max = bmax[diff];
			var max = ((1<<diff)>>1)>>>0;
			for(var c = 0; c < N; c++)
				values[i*N + c] = bitstream.read(diff) - max;
		}
		return logs.length;
	},

	decodeValues: function(N, values) {
		var t = this;
		var bitstream = t.readBitStream();
		var tunstall = new Tunstall;

		for(var c = 0; c < N; c++) {
			var logs = tunstall.decompress(this);

			for(var i = 0; i < logs.length; i++) {
				var diff = logs[i];
				if(diff == 0) {
					values[i*N + c] = 0;
					continue;
				}

				var val = bitstream.read(diff);
				var middle = 1<<(diff-1);
				if(val < middle)
					val = -val -middle;
				values[i*N + c] = val;
			}
		}
		return logs.lenght;
	}
};



function Tunstall(wordsize, lookup_size) {
	this.wordsize = wordsize? wordsize : 8;
	this.lookup_size = lookup_size? lookup_size : 8;
}

Tunstall.prototype = {
	decompress: function(stream) {
		var nsymbols = stream.readUChar();
		this.probabilities = stream.readArray(nsymbols*2);
		this.createDecodingTables();
		var size = stream.readInt();
		if(size > 100000000) throw("TOO LARGE!");
		var data = new Uint8Array(size);
		var compressed_size = stream.readInt();
		if(size > 100000000) throw("TOO LARGE!");
		var compressed_data = stream.readArray(compressed_size);
		if(size)
			this._decompress(compressed_data, compressed_size, data, size);
		return data;
	}, 

	createDecodingTables: function() {
		//read symbol,prob,symbol,prob as uchar.
		//Here probabilities will range from 0 to 0xffff for better precision

		var n_symbols = this.probabilities.length/2;
		if(n_symbols <= 1) return;

		var queues = []; //array of arrays
		var buffer = []; 

		//initialize adding all symbols to queues
		for(var i = 0; i < n_symbols; i++) {
			var symbol = this.probabilities[i*2];
			var s = [(this.probabilities[i*2+1])<<8, buffer.length, 1]; //probability, position in the buffer, length
			queues[i] = [s];
			buffer.push(this.probabilities[i*2]); //symbol
		}
		var dictionary_size = 1<<this.wordsize;
		var n_words = n_symbols;
		var table_length = n_symbols;

		//at each step we grow all queues using the most probable sequence
		while(n_words < dictionary_size - n_symbols +1) {
			//Should use a stack or something to be faster, but we have few symbols
			//find highest probability word
			var best = 0;
			var max_prob = 0;
			for(var i = 0; i < n_symbols; i++) {
				var p = queues[i][0][0]; //front of queue probability.
				if(p > max_prob) {
					best = i;
					max_prob = p;
				}
			}
			var symbol = queues[best][0];
			var pos = buffer.length;
			
			for(var i = 0; i < n_symbols; i++) {
				var sym = this.probabilities[i*2];
				var prob = this.probabilities[i*2+1]<<8;
				var s = [((prob*symbol[0])>>>16), pos, symbol[2]+1]; //combine probabilities, keep track of buffer, keep length of queue

				for(var k  = 0; k < symbol[2]; k++)
					buffer[pos+k] = buffer[symbol[1] + k]; //copy sequence of symbols

				pos += symbol[2];
				buffer[pos++] = sym; //append symbol
				queues[i].push(s);
			}
			table_length += (n_symbols-1)*(symbol[2] + 1) +1; 
			n_words += n_symbols -1;
			queues[best].shift(); //remove first thing
		}

		this.index = new Uint32Array(n_words);
		this.lengths = new Uint32Array(n_words);
		this.table = new Uint8Array(table_length);
		var word = 0;
		var pos = 0;
		for(i = 0; i < queues.length; i++) {
			var queue = queues[i];
			for(var k = 0; k < queue.length; k++) {
				var s = queue[k];
				this.index[word] = pos;
				this.lengths[word] = s[2]; //length
				word++;

				for(var j = 0; j < s[2]; j++)
					this.table[pos + j] = buffer[s[1] + j]; //buffer of offset
				pos += s[2]; //length
			}
		}
	},
	_decompress: function(input, input_size, output, output_size) {
		//TODO optimize using buffer arrays
		var input_pos = 0;
		var output_pos = 0;
		if(this.probabilities.length == 2) {
			var symbol = this.probabilities[0];
			for(var i = 0; i < output_size; i++)
				output[i] = symbol;
			return;
		}

		while(input_pos < input_size-1) {
			var symbol = input[input_pos++];
			var start = this.index[symbol];
			var end = start + this.lengths[symbol];
			for(var i = start; i < end; i++) 
				output[output_pos++] = this.table[i];
		}

		//last symbol might override so we check.
		var symbol = input[input_pos];
		var start = this.index[symbol];
		var end = start + output_size - output_pos;
		var length = output_size - output_pos;
		for(var i = start; i < end; i++)
			output[output_pos++] = this.table[i];

		return output;
	}
}
