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
	t.max = new Uint32Array(32);
	for(var i = 0; i < 32; i++)
		t.max[i] = (1<<i)>>>1;
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
//			var max = t.max[diff];
			var max = (1<<diff)>>>1;
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
				var middle = (1<<(diff-1))>>>0;
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

/*function Tunstall() {
}

Tunstall.prototype.wordsize = 8;
Tunstall.prototype.lookup_size = 8; */



Tunstall.prototype = {
	decompress: function(stream) {
		var nsymbols = stream.readUChar();
		this.probs = stream.readArray(nsymbols*2);
		this.createDecodingTables2();
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
		//Here probs will range from 0 to 0xffff for better precision

		var n_symbols = this.probs.length/2;
		if(n_symbols <= 1) return;

		var queues = []; //array of arrays
		var buffer = []; 


		//initialize adding all symbols to queues
		for(var i = 0; i < n_symbols; i++) {
			var symbol = this.probs[i*2];
			var s = [(this.probs[i*2+1])<<8, buffer.length, 1]; //probability, position in the buffer, length
			queues[i] = [s];
			buffer.push(this.probs[i*2]); //symbol
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
				var sym = this.probs[i*2];
				var prob = this.probs[i*2+1]<<8;
				var s = [((prob*symbol[0])>>>16), pos, symbol[2]+1]; //combine probs, keep track of buffer, keep length of queue

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
	createDecodingTables2: function() {
		//Here probs will range from 0 to 0xffff for better precision

		var t = this;
		var n_symbols = t.probs.length/2;
		if(n_symbols <= 1) return;

/*      QUEUE as column 

		A  x  AA
		C  C  CA
		G  G  GA
		T  T  TA
*/

		var dictionary_size = 1<<this.wordsize;

		//TODO: fill index and length: saves *3, and memory (compact in place, later)

		//(max length is word dictionary? nope the removed ones dont count, 1 per step so add max 255 potentially deleted)
		var queues = new Uint32Array(2*dictionary_size*3); //array of n symbols array of prob, position in buffer, length //limit to 8k*3
		var starts = new Uint32Array(n_symbols); //starts of each queue
		var end = 0; //keep track of queue end
		var buffer = t.table = new Uint8Array(dictionary_size*dictionary_size); //worst case for 2 symbols. 1 with 254 prob and other qwith 1 prob
		var pos = 0; //keep track of buffer first free space
		var table_length = 0;
		var n_words = 0;

		var count = 2;
		var p0 = (t.probs[1] << 8)>>>0; //first symbol
		var p1 = (t.probs[3]<<8)>>>0;    //second symbol.
		var prob = (p0*p0)>>>16;
		var max_count = (dictionary_size - 1)/(n_symbols - 1);
		while(prob > p1 && count < max_count) {
			prob = (prob*p0) >>> 16;
			count++;
		}

		if(count >= 16) { //Very low entropy results in large tables > 8K.
			buffer[pos++] = t.probs[0];
			for(var k = 1; k < n_symbols; k++) {
				for(var i = 0; i < count-1; i++)
					buffer[pos++] = t.probs[0];
				buffer[pos++] = t.probs[2*k];
			}
			starts[0] = (count-1)*n_symbols*3;
			for(var k = 1; k < n_symbols; k++)
				starts[k] = k*3;

			prob = 0xffff;
			for(var col = 0; col < count; col++) {
				for(var row = 1; row < n_symbols; row++) {
					var off = (row + col*n_symbols)*3;
					queues[off] = (prob * (t.probs[row*2+1]<<8)) >> 16;
					queues[off+1] = row*count - col;
					queues[off+1] = col+1;
				}
				prob = (prob*p0) >>> 16;
			}
			var first = ((count-1)*n_symbols)*3;
			queues[first] = prob;
			queues[first+1] = 0;
			queues[first+1] = count;
			n_words = 1 + count*(n_symbols - 1);
			table_length = n_words;
			end = count*n_symbols;

		} else {
			//initialize adding all symbols to queues
			for(var i = 0; i < n_symbols; i++) {
				var symbol = t.probs[i*2];
				var prob = t.probs[i*2+1];
				queues[i*3] = prob<<8;
				queues[i*3+1] = i;
				queues[i*3+2] = 1;
	
				starts[i] = i*3;
				buffer[i] = symbol;
			}
			pos = n_symbols;
			end = n_symbols*3;

			n_words = n_symbols;
			table_length = n_symbols;
		}

		//at each step we grow all queues using the most probable sequence
		while(n_words < dictionary_size - n_symbols +1) {
			//find highest probability word
			var best = 0;
			var max_prob = 0;
			for(var i = 0; i < n_symbols; i++) {
				var p = queues[starts[i]]; //front of queue probability.
				if(p > max_prob) {
					best = i;
					max_prob = p;
				}
			}
			var start = starts[best];
			var offset = queues[start+1];
			var len = queues[start+2];
			for(var i = 0; i < n_symbols; i++) {
				var sym = t.probs[i*2];
				var prob = t.probs[i*2+1]<<8;
				queues[end] = ((prob*queues[start])>>>16);
				queues[end+1] = pos;
				queues[end+2] = len + 1;
				end += 3;

				for(var k  = 0; k < len; k++)
					buffer[pos + k] = buffer[offset + k]; //copy sequence of symbols

				pos += len;
				buffer[pos++] = sym; //append symbol
			}
			starts[best] += n_symbols*3; //move one column

			table_length += (n_symbols-1)*(len + 1) +1;
			n_words += n_symbols -1;
		}

		t.index = new Uint32Array(n_words);
		t.lengths = new Uint32Array(n_words);
		var word = 0;
		var count = 0;
		for(i = 0, row = 0; i < end; i += 3, row++) {
			if(row >= n_symbols)
				row  = 0;
			if(starts[row] > i) continue; //skip deleted words

			t.index[word] = queues[i+1];
			t.lengths[word] = queues[i+2];
			word++;
		}
	},
	_decompress: function(input, input_size, output, output_size) {
		//TODO optimize using buffer arrays
		var input_pos = 0;
		var output_pos = 0;
		if(this.probs.length == 2) {
			var symbol = this.probs[0];
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
