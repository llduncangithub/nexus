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

//TODO rename format in tyep

var Format

function Attribute(name, q, components, format, strategy) {
	var t = this;
	t.name = name;
	t.q = q; //float
	t.components = components; //integer
	t.format = format;
	t.strategy = strategy;
}

Attribute.prototype = {
init: function(nvert, nface) {
	var t = this;
	var n = nvert*t.components;
	t.values = new Int32Array(n);  //local workspace 

	//init output buffers
	switch(t.format) {
	case 0:
	case 1: t.values = t.buffer = new Int32Array(n); break; //no point replicating.
	case 2:
	case 3: t.buffer = new Int16Array(n); break;
	case 4: t.buffer = new Int8Array(n); break;
	case 5: t.buffer = new Uint8Array(n); break;
	case 6:
	case 7: t.buffer = new Float32Array(n); break;
	default: throw "Error if reading";
	}
},

decode: function(nvert, stream) {
	var t = this;
	if(t.strategy & 0x2) //correlated
		stream.encodeArray(t.values, t.components);
	else
		stream.decodeValues(t.values, t.components);
},

deltaDecode: function(nvert, context) {
	var values = t.values;
	var N = t.N;
	var n = context.length/3;
	if(strategy & 0x1) { //parallel
		for(var i = 1; i < n; i++)
			for(var c = 0; c < N; c++)
				values[i*N + c] = values[context[i*3]*N + c] + values[context[i*3+1]*N + c] - values[context[i*3+2]*N + c];
	} else if(context) {
		for(var i = 1; i < n; i++)
			for(var c = 0; c < N; c++)
				values[i*N + c] += values[context[i*3]*N + c];
	} else {
		for(var i = N; i < nvert*N; i++)
				values[i] += values[i - N];
	}
},

postDelta: function() {},

dequantize: function(nvert) {
	var t= this;
	var n = t.n*nvert;
	switch(t.format) {
	case 0:
	case 1: break;
	case 2:
	case 3: 
	case 4: 
	case 5: 
		for(var i = 0; i < n; i++)
			buffer[i] = coords[i];
		break;
	case 6:
	case 7: 
		for(var i = 0; i < n; i++)
			buffer[i] = coords[i]*q;
		break;
	}
}

};

function ColorAttr(name, q, components, format, strategy) {
	Attribute.call(this, name, q, components, format, strategy);
	this.qc = [];
}

ColorAttr.prototype = Object.create(Attribute.prototype);
ColorAttr.prototype.decode = function(nvert, stream) {
	for(var c = 0; c < 4; c++)
		t.qc[c] = stream.readUChar();
	Attribute.prototype.decode.call(this, nvert, stream);
}

ColorAttr.prototype.dequantize = function(nvert) {
	var t = this;
	for(var i = 0; i < nvert; i++) {
		var offset = i*4;

		var e0 = t.values[offset + 0] * t.qc[0];
		var e1 = t.values[offset + 1] * t.qc[1];
		var e2 = t.values[offset + 2] * t.qc[2];

		t.buffer[offset + 0] = (e2 + e0)&0xff;
		t.buffer[offset + 1] = e0;
		t.buffer[offset + 2] = (e1 + e0)&0xff;
		t.buffer[offset + 3] = t.values[offset + 3] * t.qc[3];
	}
}

function NormalAttr(name, q, components, format, strategy) {
	Attribute.call(this, name, q, components, format, strategy);
}

NormalAttr.prototype = Object.create(Attribute.prototype);
NormalAttr.prototype.init = function(nvert, nface) {
	var t = this;
	var n = nvert*t.components;
	t.values = new Int32Array(2);  //local workspace 

	//init output buffers
	switch(t.format) {
	case 3: t.buffer = new Int16Array(n); break;
	case 6:
	case 7: t.buffer = new Float32Array(n); break;
	default: throw "Error if reading";
	}
},

NormalAttr.prototype.toSphere = function(o, values, buffer, unit) {
	var av0 = v[o] > 0? v[o]:-v[o];
	var av1 = v[o+1] > 0? v[o+1]:-v[o+1];
	buffer[o] = values[o];
	buffer[o+1] = values[o+1];
	buffer[o+2] = unit - av0 - av1;
	if (buffer[o+2] < 0) {
		buffer[o] = ((v[o] > 0)? 1 : -1)*(unit - av1);
		buffer[o+1] = ((v[o+1] > 0)? 1 : -1)*(unit - av1[0]);
	}
	var len = 1/Math.sqrt(buffer[o]*buffer[o] + buffer[o+1]*buffer[o+1] + buffer[o+2]*buffer[o+2]);
	if(t.format == 3) {
		len *= 32767;
	}
	buffer[o] *= len;
	buffer[o+1] *= len;
	buffer[o+2] *= len;
}

NormalAttr.prototype.decode = function(nvert, stream) {
	var t = this;
	t.prediction = stream.readUChar();

	stream.decodeArray(t.values, 2);
}

function IndexAttr(nvert, nface, format) {
	var t = this;
	if(!format || format == 0) //uint32 by default
		t.faces = new Uint32Array(nface);
	else if(format == 2)
		t.faces = new Uint16Buffer(nface);
	else
		throw "Unsupported format";
}

IndexAttr.prototype = {
decode: function(stream) {
	var t = this;

	var n = stream.readInt();
	t.groups = new Array(n);
	for(var i = 0; i < n; i++)
		t.groups[i] = stream.readInt();
		var stunstall = new Tunstall;

	var tunstall = new Tunstall;
	t.clers = tunstall.decompress(stream);
	t.bitstream = stream.readBitStream();
}
};



function NxzDecoder(data) {
	var t = this;
	var stream = t.stream = new Stream(data);
	
	var magic = stream.readInt();
	if(magic != 2021286656) return;

	var version = stream.readInt();
	t.entropy = stream.readUChar();
	var n = stream.readInt();

	t.geometry = {};
	t.attributes = {};
	for(let i = 0; i < n; i++) {
		var a = {};
		var name = stream.readString();
		var q = stream.readFloat();
		var components = stream.readUChar(); //internal number of components
		var format = stream.readUChar();     //default format (same as it was in input), can be overridden
		var strategy = stream.readUChar();
		var attr = Attribute;
		switch(name) {
		case "color":  attr = ColorAttr; break;
		case "normal": attr = NormalAttr; break;
		}
		t.attributes[name] = new attr(name, q, components, format, strategy);
	}

//TODO move this vars into an array.
	t.geometry.nvert = t.nvert = t.stream.readInt();
	t.geometry.nface = t.nface = t.stream.readInt();

	console.log(t.nvert, t.nface);
}

NxzDecoder.prototype = {

decode: function() {
	var t = this;

	t.last = new Uint32Array(t.nvert*3); //for parallelogram prediction
	t.last_count = 0;

	for(let i in t.attributes)
		t.attributes[i].init(t.nvert, t.nface);

	if(t.nface == 0)
		t.decodePointCloud();
	else
		t.decodeMesh();

	return t.geometry;
},

decodePointCloud: function() {
	var t = this;
	for(var i in t.attributes) {
		var a = t.attributes[i];
		a.decode(t.nvert, t.stream);
		a.deltaDecode(t.nvert);
		a.dequantize(t.nvert);
		t.geometry[a.name] = a.buffer;
	}
},

decodeMesh: function() {
	var t = this;
	t.index = new IndexAttr(t.nvert, t.nface, 0); //USE
	t.index.decode(t.stream);

	t.vertex_count = 0;
	var start = 0;
	var cler = 0;
	for(var p = 0; p < this.patches.length; p++) {
		var end = this.patches[p];
		this.decodeConnectivity(end - start, start*3);
		start = end;
	}
	//dequantize positions
	var tot = this.node.nvert*3;
	var coords = this.coords;
	var stack = this.stack;
	for(var i = 0; i < tot; ) {
		coords[i] = (coords[i] + stack[3])*stack[6]; i++;
		coords[i] = (coords[i] + stack[4])*stack[6]; i++;
		coords[i] = (coords[i] + stack[5])*stack[6]; i++;
	}
	if(this.sig.texcoords) {
		var t_tot = this.node.nvert*2;
		var t_coords = this.texcoords;
		for(var i = 0; i < tot; ) {
			t_coords[i] = (t_coords[i] + stack[9])*stack[11]; i++;
			t_coords[i] = (t_coords[i] + stack[10])*stack[11]; i++;
		}		
	}
},

decodeNormals: function() {
	var norm_q = this.stream.readChar();

	var dtunstall = new Tunstall;
	var diffs = dtunstall.decompress(this.stream);

	var stunstall = new Tunstall;
	var signs = stunstall.decompress(this.stream);
	var bitstream = this.stream.readBitStream();

	var side = (1<<(16 - norm_q))>>>0;
	var diffcount = 0;
	var signcount = 0;

	if(!this.sig.indices) {
		for(var k = 0; k < 2; k++) {
			var on = 0;
			for(var i = 0; i < this.node.nvert; i++) {
				var d = this.decodeDiff(diffs[diffcount++], bitstream);
				on = on + d;
				this.normals[3*i + k] = on*side;
			}
		}
		for(var i = 0; i < this.node.nvert; i++) {
			var offset = i*3;
			var x = this.normals[offset + 0];
			var y = this.normals[offset + 1];
			var z = 32767.0*32767.0 - x*x - y*y;

        	if(z < 0) z = 0;
        	z = Math.sqrt(z);
			if(z > 32767) z = 32767;
			if(signs[i] == 0)
				z = -z;
			this.normals[offset + 2] = z;
		}
		return;
	}

	var boundary = this.markBoundary();
	this.computeNormals();

	if(this.sig.texcoords) //hack, fixing normals makes it worse actually
		return;

	var stat = 0;
	//get difference between original and predicted
	for(var i = 0; i < this.node.nvert; i++) {
		if(!boundary[i]) continue;
		var offset = i*3;
		var x = (this.normals[offset + 0]/side);
		var y = (this.normals[offset + 1]/side);
		var dx = this.decodeDiff(diffs[diffcount++], bitstream);
		var dy = this.decodeDiff(diffs[diffcount++], bitstream);
		x = (x + dx)*side;
		y = (y + dy)*side;

        var z = 32767.0*32767.0 - x*x - y*y;

        if(z < 0) z = 0;
        z = Math.sqrt(z);
        //sign
        if(z > 32767.0) z = 32767.0;
        var signbit = signs[signcount++];
//        if(this.normals[offset+2] < 0 != signbit)
        if((this.normals[offset+2] < 0 && signbit == 0) || (this.normals[offset+2] > 0 && signbit == 1))
        	z = -z;
		this.normals[offset + 0] = x;
		this.normals[offset + 1] = y;
        this.normals[offset + 2] = z;
	}
},

decodeColors: function() {
	var color_q = [];
	for(var k = 0; k < 4; k++)
		color_q[k] = this.stream.readChar();

	var diffs = [];
	for(var k = 0; k < 4; k++) {
		var tunstall = new Tunstall;;
		diffs[k] = tunstall.decompress(this.stream);
	}
	var bitstream = this.stream.readBitStream();

	var count = 0;
	if(this.sig.indices) {
		for(var i = 0; i < this.node.nvert; i++) {
			var last  = this.last[i]*4;
			var offset = i*4;

			for(var k = 0; k < 4; k++) {
				var c = this.decodeDiff(diffs[k][count], bitstream);

				if(last >= 0)
					c += this.colors[last + k];
				this.colors[offset] = c;
				offset++;
			}
			count++;
		}
	} else {
		for(var k = 0; k < 4; k++)
			this.colors[k] = this.decodeDiff(diffs[k][count], bitstream); 
		count++;

		var offset = 4;
		for(var i = 1; i < this.node.nvert; i++) {
			for(var k = 0; k < 4; k++) {
				var d = this.decodeDiff(diffs[k][count], bitstream); 
				this.colors[offset] = this.colors[offset-4] + d;
				offset ++;
			}
			count++;
		}
	}

	var steps = [];
	for(var k = 0; k < 4; k++)
		steps[k] = (1<<(8 - color_q[k]));

	//convert to rgb
	for(var i = 0; i < this.node.nvert; i++) {
		var offset = i*4;

		var e0 = this.colors[offset + 0] * steps[0];
		var e1 = this.colors[offset + 1] * steps[1];
		var e2 = this.colors[offset + 2] * steps[2];

		this.colors[offset + 0] = (e2 + e0)&0xff;
		this.colors[offset + 1] = e0;
		this.colors[offset + 2] = (e1 + e0)&0xff;
	}
},

//how to determine if a vertex is a boundary without topology:
//for each edge a vertex is in, add or subtract the id of the other vertex depending on order
//for internal vertices sum is zero.
//unless we have strange configurations and a lot of sfiga, zero wont happen. //TODO think about this
markBoundary: function() {
//	var boundary = new Uint8Array(this.node.nvert);
	var count = new Uint32Array(this.node.nvert);

	var offset = 0;
	for(var i = 0; i < this.node.nface; i++) {
		count[this.faces[offset + 0]] += this.faces[offset + 1] - this.faces[offset + 2];
		count[this.faces[offset + 1]] += this.faces[offset + 2] - this.faces[offset + 0];
		count[this.faces[offset + 2]] += this.faces[offset + 0] - this.faces[offset + 1];
		offset += 3;
	}
	return count;
//	for(var i = 0; i < this.node.nvert; i++)
//		if(count[i] != 0)
//			boundary[i] = true;
//	return boundary;
},

norm: function(buffer, a, b, c) { //a b c offsets in the buffer
	var ba0 = buffer[b+0] - buffer[a+0];
	var ba1 = buffer[b+1] - buffer[a+1];
	var ba2 = buffer[b+2] - buffer[a+2];

	var ca0 = buffer[c+0] - buffer[a+0];
	var ca1 = buffer[c+1] - buffer[a+1];
	var ca2 = buffer[c+2] - buffer[a+2];

	var p = [];
	p[0] = ba1*ca2 - ba2*ca1;
	p[1] = ba2*ca0 - ba0*ca2;
	p[2] = ba0*ca1 - ba1*ca0;
	return p;
},

normalize: function(buffer, offset) {
	var x = buffer[offset + 0];
	var y = buffer[offset + 1];
	var z = buffer[offset + 2];
	var n = Math.sqrt(x*x + y*y + z*z);
	if(n > 0) {
		buffer[offset + 0] = x/n;
		buffer[offset + 1] = y/n;
		buffer[offset + 2] = z/n;
	}
},

computeNormals:function() {
	var tmp_normals = new Float32Array(this.node.nvert*3);

	var offset = 0;
	for(var i = 0; i < this.node.nface; i++) {
		var a = 3*this.faces[offset + 0];
		var b = 3*this.faces[offset + 1];
		var c = 3*this.faces[offset + 2];

		var buffer = this.coords;
		var ba0 = buffer[b+0] - buffer[a+0];
		var ba1 = buffer[b+1] - buffer[a+1];
		var ba2 = buffer[b+2] - buffer[a+2];

		var ca0 = buffer[c+0] - buffer[a+0];
		var ca1 = buffer[c+1] - buffer[a+1];
		var ca2 = buffer[c+2] - buffer[a+2];

		var n0 = ba1*ca2 - ba2*ca1;
		var n1 = ba2*ca0 - ba0*ca2;
		var n2 = ba0*ca1 - ba1*ca0;

		tmp_normals[a + 0] += n0;
		tmp_normals[a + 1] += n1;
		tmp_normals[a + 2] += n2;
		tmp_normals[b + 0] += n0;
		tmp_normals[b + 1] += n1;
		tmp_normals[b + 2] += n2;
		tmp_normals[c + 0] += n0;
		tmp_normals[c + 1] += n1;
		tmp_normals[c + 2] += n2;
		offset += 3;
	}

	//normalize
	var offset = 0;
	for(var i = 0; i < this.node.nvert; i++) {
		var x = tmp_normals[offset + 0];
		var y = tmp_normals[offset + 1];
		var z = tmp_normals[offset + 2];
		var n = Math.sqrt(x*x + y*y + z*z);
		if(n > 0) {
			tmp_normals[offset + 0] = x/n;
			tmp_normals[offset + 1] = y/n;
			tmp_normals[offset + 2] = z/n;
		}
		this.normals[offset + 0] = tmp_normals[offset + 0]*32767;
		this.normals[offset + 1] = tmp_normals[offset + 1]*32767;
		this.normals[offset + 2] = tmp_normals[offset + 2]*32767;
		offset += 3;
	}
},

decodeDiff: function(diff, bitstream) {
	var val;
	if(diff == 0) {
		val = 1;
	} else {
		val = 1<<(diff);
		val |= bitstream.read(diff);
	};
	val--; //vall is always >= 1
	if(val & 0x1)
		val = -((val+1)>>1);
	else
		val = val>>1;
	return val;
},

/* an edge is:   uint16_t face, uint16_t side, uint32_t prev, next, bool deleted
I do not want to create millions of small objects, I will use aUint32Array.
Problem is how long, sqrt(nface) we will over blow using nface.
*/

decodeConnectivity: function(length, start) {

	var t = this;
	var ctunstall = new Tunstall;
	var clers = ctunstall.decompress(this.stream);
	var cler_count = 0;

	var dtunstall = new Tunstall;
	var diffs = dtunstall.decompress(this.stream);
	var diff_count = 0;

	var tdiffs;
	var tdiff_count = 0;
	if(t.sig.texcoords) {
		var ttunstall = new Tunstall;
		tdiffs = ttunstall.decompress(this.stream);	
	}

	var bitstream = this.stream.readBitStream(bitstream);

	var current_face = 0;          //keep track of connected component start
	//t.vertex_count = 0;
	var front = new Uint32Array(this.node.nface*18);
	var front_count = 0; //count each integer so it's front_back*5
	function addFront(_v0, _v1, _v2, _prev, _next) {
		front[front_count++] = _v0;
		front[front_count++] = _v1;
		front[front_count++] = _v2;
		front[front_count++] = _prev;
		front[front_count++] = _next;
		front[front_count++] = 0; //deleted
	}
	function _next(t) {
		t++;
		if(t == 3) t = 0;
		return t;
	}
	function _prev(t) {
		t--;
		if(t == -1) t = 2;
		return t;
	}

	var delayed = [];
	var faceorder = [];

	var faces_count = start; //count indices in this.faces array
	var totfaces = length;
//	var estimated = [0, 0, 0]; //no! use stack.
	var stack = this.stack;
	var coords = this.coords;
	var texcoords = this.texcoords;
	var hasTexCoords = t.sig.texcoords;

	while(totfaces > 0) {
		if(!faceorder.length && !delayed.length) {
			if(current_face == this.node.nface) break; //no more faces to encode exiting

			stack[0] = stack[1] = stack[2] = 0;
			stack[7] = stack[8] = 0; //texcoords
			var last_index = -1;
			var index = [];
			for(var k = 0; k < 3; k++) {
				this.last[this.last_count++] = last_index;
				var diff = diffs[diff_count++];
				var tdiff = diff && hasTexCoords? tdiffs[tdiff_count++] : 0;
				var v = this.decodeVertex(bitstream, diff, tdiff);
				index[k] = v; 
				this.faces[faces_count++] = v;
				stack[0] = coords[v*3];
				stack[1] = coords[v*3+1];
				stack[2] = coords[v*3+2]; 
				if(t.sig.texcoords) {
					stack[7] = texcoords[v*2];
					stack[8] = texcoords[v*2+1];
				}
				last_index = v;
			}
			var current_edge = front_count;
			for(var k = 0; k < 3; k++) {
				faceorder.push(front_count);
				front[front_count++] = index[_next(k)];
				front[front_count++] = index[_prev(k)]; 
				front[front_count++] = index[k];
				front[front_count++] = current_edge + _prev(k)*6;
				front[front_count++] = current_edge + _next(k)*6;
				front_count++;
//				addFront(index[_next(k)], index[_prev(k)], index[k], current_edge + _prev(k)*6, current_edge + _next(k)*6);
			}
			current_face++;
			totfaces--;
			continue;
		}
		var f;
		if(faceorder.length) 
			f = faceorder.shift();
		else 
			f = delayed.pop();
		
		var edge_start = f;

		if(front[edge_start + 5]) continue; //deleted
		front[edge_start + 5] = 1; //set edge as deleted anyway

		var c = clers[cler_count++];
		if(c == 4) continue; //BOUNDARY

		var v0   = front[edge_start + 0];
		var v1   = front[edge_start + 1];
		var v2   = front[edge_start + 2];
		var prev = front[edge_start + 3];
		var next = front[edge_start + 4];

		var first_edge = front_count; //points to new edge to be inserted
		var opposite = -1;
		if(c == 0) { //VERTEX
			//predict position based on v0, v1 and v2
			for(var k = 0; k < 3; k++) 
				stack[k] = coords[v0*3 + k] + coords[v1*3 + k] - coords[v2*3 + k];

			if(hasTexCoords)
				for(var k = 0; k < 2; k++)
					stack[7+k] = texcoords[v0*2 + k]  + texcoords[v1*2 + k] - texcoords[v2*2 + k];
			
			var diff = diffs[diff_count++];
			var tdiff = diff && hasTexCoords? tdiffs[tdiff_count++] : 0;
			opposite = this.decodeVertex(bitstream, diff, tdiff);
			if(diff != 0)
				this.last[this.last_count++] = v1;

			front[prev + 4] = first_edge;
			front[next + 3] = first_edge + 6;
			faceorder.unshift(front_count);

			front[front_count++] = v0;
			front[front_count++] = opposite;
			front[front_count++] = v1;
			front[front_count++] = prev;
			front[front_count++] = first_edge+6;
			front_count++; 
//			addFront(v0, opposite, v1, prev, first_edge + 6);

			faceorder.push(front_count);

			front[front_count++] = opposite;
			front[front_count++] = v1;
			front[front_count++] = v0;
			front[front_count++] = first_edge; 
			front[front_count++] = next;
			front_count++; 
//			addFront(opposite, v1, v0, first_edge, next);

		} else if(c == 3) { //END
			front[prev + 5] = 1;
			front[next + 5] = 1;
			front[front[prev + 3] + 4] = front[next + 4];
			front[front[next + 4] + 3] = front[prev + 3];
			opposite = front[prev + 0];

		} else if(c == 1) { //LEFT
			front[prev + 5] = 1; //deleted
			front[front[prev + 3] + 4] = first_edge;
			front[next + 3] = first_edge;
			opposite = front[prev + 0];

			faceorder.unshift(front_count);

			front[front_count++] = opposite;
			front[front_count++] = v1;
			front[front_count++] = v0;
			front[front_count++] = front[prev +3];
			front[front_count++] = next;
			front_count++; 
//			addFront(opposite, v1, v0, front[prev + 3], next);

		} else if(c == 2) { //RIGHT
			front[next + 5] = 1;
			front[front[next + 4] + 3] = first_edge;
			front[prev + 4] = first_edge;
			opposite = front[next + 1];


			faceorder.unshift(front_count);

			front[front_count++] = v0;
			front[front_count++] = opposite;
			front[front_count++] = v1;
			front[front_count++] = prev;
			front[front_count++] = front[next+4];
			front_count++; 
//			addFront(v0, opposite, v1, prev, front[next + 4]);

		} else if(c == 5) { //DELAY
			front[edge_start + 5] = 0;
			delayed.push(edge_start);
			continue;
		}
		this.faces[faces_count++] = v1;
		this.faces[faces_count++] = v0;
		this.faces[faces_count++] = opposite;
		totfaces--;
	}
},
   
decodeVertex: function(bitstream, diff, tdiff) {
	if(diff == 0) 
		return bitstream.read(16);

	var v = this.vertex_count++;

	var max = 1<<(diff-1);

	for(var k = 0; k < 3; k++) {
		var d = bitstream.read(diff) - max;
		this.coords[v*3+k] = this.stack[k] + d; //stack 0-3 is used as extimated
	}
	if(this.sig.texcoords) {
		var tmax = 1<<(tdiff-1);
		for(var k = 0; k < 2; k++) {
			var d = bitstream.read(tdiff) - tmax;
			this.texcoords[v*2+k] = this.stack[7+k] + d; //stack 7-9 is used as extimated
		}
	}
	return v;
},

decodeDiff: function(diff, bitstream) {
	var val;
	if(diff == 0) {
		return 0;
	} 
	val = 1<<diff;
	val += bitstream.read(diff);


	if(val & 0x1) 
		val >>>= 1;
	else 
		val = -(val>>>1);

	return val;
}

};

var tot = 0;


