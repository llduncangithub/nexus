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




function Attribute(name, q, components, type, strategy) {
	var t = this;
	t.name = name;
	t.q = q; //float
	t.components = components; //integer
	t.type = type;
	t.strategy = strategy;
}

Attribute.prototype = {

Type: { UINT32:0, INT32:1, UINT16:2, INT16:3, UINT8:4, INT8:5, FLOAT:6, DOUBLE:7 }, 

Strategy: { PARALLEL:1, CORRELATED:2 },

init: function(nvert, nface) {
	var t = this;
	var n = nvert*t.components;
	t.values = new Int32Array(n);  //local workspace 

	//init output buffers
	switch(t.type) {
	case t.Type.UINT32:
	case t.Type.INT32: t.values = t.buffer = new Int32Array(n); break; //no point replicating.
	case t.Type.UINT16:
	case t.Type.INT16: t.buffer = new Int16Array(n); break;
	case t.Type.UINT8: t.buffer = new Uint8Array(n); break;
	case t.Type.INT8: t.buffer  = new Int8Array(n); break;
	case t.Type.FLOAT:
	case t.Type.DOUBLE: t.buffer = new Float32Array(n); break;
	default: throw "Error if reading";
	}
},

decode: function(nvert, stream) {
	var t = this;
	if(t.strategy & t.Strategy.CORRELATED) //correlated
		stream.decodeArray(t.components, t.values);
	else
		stream.decodeValues(t.components, t.values);
},

deltaDecode: function(nvert, context) {
	var t = this;
	var values = t.values;
	var N = t.components;

	if(t.strategy & t.Strategy.PARALLEL) { //parallel
		var n = context.length/3;
		for(var i = 1; i < n; i++) {
			for(var c = 0; c < N; c++) {
				values[i*N + c] += values[context[i*3]*N + c] + values[context[i*3+1]*N + c] - values[context[i*3+2]*N + c];
			}
		}
	} else if(context) {
		var n = context.length/3;
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
	var n = t.components*nvert;
	switch(t.type) {
	case t.Type.UINT32:
	case t.Type.INT32: break;
	case t.Type.UINT16:
	case t.Type.INT16: 
	case t.Type.UINT8: 
	case t.Type.INT8: 
		for(var i = 0; i < n; i++)
			t.buffer[i] = t.values[i]*t.q;
		break;
	case t.Type.FLOAT:
	case t.Type.DOUBLE: 
		for(var i = 0; i < n; i++)
			t.buffer[i] = t.values[i]*t.q;
		break;
	}
}

};


/* COLOR ATTRIBUTE */

function ColorAttr(name, q, components, type, strategy) {
	Attribute.call(this, name, q, components, type, strategy);
	this.qc = [];
}

ColorAttr.prototype = Object.create(Attribute.prototype);
ColorAttr.prototype.decode = function(nvert, stream) {
	for(var c = 0; c < 4; c++)
		this.qc[c] = stream.readUChar();
	Attribute.prototype.decode.call(this, nvert, stream);
}

ColorAttr.prototype.dequantize = function(nvert) {
	var t = this;
	for(var i = 0; i < nvert; i++) {
		var offset = i*4;
		var rgboff = i*3;

		var e0 = t.values[offset + 0];
		var e1 = t.values[offset + 1];
		var e2 = t.values[offset + 2];

		t.buffer[rgboff + 0] = ((e2 + e0)* t.qc[0])&0xff;
		t.buffer[rgboff + 1] = e0* t.qc[1];
		t.buffer[rgboff + 2] = ((e1 + e0)* t.qc[2])&0xff;
//		t.buffer[offset + 0] = t.values[offset + 3] * t.qc[3];
	}
}

/* NORMAL ATTRIBUTE */

function NormalAttr(name, q, components, type, strategy) {
	Attribute.call(this, name, q, components, type, strategy);
}

NormalAttr.prototype = Object.create(Attribute.prototype);

NormalAttr.prototype.Prediction = { DIFF: 0, ESTIMATED: 1, BORDER: 2 };

NormalAttr.prototype.init = function(nvert, nface) {
	var t = this;
	var n = nvert*t.components;
	t.values = new Int32Array(2*nvert);  //local workspace 

	//init output buffers
	switch(t.type) {
	case t.Type.INT16: t.buffer = new Int16Array(n); break;
	case t.Type.FLOAT:
	case t.Type.DOUBLE: t.buffer = new Float32Array(n); break;
	default: throw "Error if reading";
	}
};

NormalAttr.prototype.decode = function(nvert, stream) {
	var t = this;
	t.prediction = stream.readUChar();

	stream.decodeArray(2, t.values);
};

NormalAttr.prototype.deltaDecode = function(nvert, context) {
	var t = this;
	if(t.prediction != t.Prediction.DIFF)
		return;

	if(context) {
		for(var i = 1; i < nvert; i++) {
			for(var c = 0; c < 2; c++) {
				var d = t.values[i*2 + c];
				t.values[i*2 + c] += t.values[context[i*3]*2 + c];
			}
		}
	} else { //point clouds assuming values are already sorted by proximity.
		for(var i = 2; i < nvert*2; i++) {
			var d = t.values[i];
			t.values[i] += t.values[i-2];
		}
	}
};

NormalAttr.prototype.postDelta = function(nvert, nface, attrs, index) {
	var t = this;
	//for border and estimate we need the position already deltadecoded but before dequantized
	if(t.prediction == t.Prediction.DIFF)
		return;

	if(!attrs.position)
		throw "No position attribute found. Use DIFF normal strategy instead.";

	var coord = attrs.position;

	var estimated = new Float32Array(nvert*3);
	t.estimateNormals(nvert, coord.values, nface, index.values, estimated);

	if(t.prediction == t.Precition.BORDER) {
		t.boundary = new Uint8Array(nvert);
		t.markBoundary(nvert, nface, index, boundary);
	}

	t.computeNormals(nvert, estimated);
/*	switch(format) {
	case FLOAT:
		computeNormals((Point3f *)buffer, estimated);
		break;
	case INT16:
		computeNormals((Point3s *)buffer, estimated);
		break;
	default: throw "Format not supported for normal attribute (float, int16 only)";
	}*/
} 



NormalAttr.prototype.dequantize = function(nvert) {
	var t = this;
	if(t.prediction != t.Prediction.DIFF)
		return;

	for(var i = 0; i < nvert; i++)
		t.toSphere(i, t.values, t.buffer, t.q)
}


NormalAttr.prototype.computeNormals(nvert, estimated) {

	var count = 0; //here for the border.
	if(t.prediction == t.Prediction.ESTIMATED) {
		for(var i = 0; i < nvert; i++) {
			t.toOcta(i, estimated, values, t.q);
			t.toSphere(i, t.buffer, t.values, (int)q);
		}
	} else { //BORDER
		for(var i = 0; i < nvert; i++) {
			Point3f &e = estimated[i];
			int32_t *d = &diffs[count*2];
			Point3s &n = normals[i];

			if(boundary[i]) {
				Point2i qn = toOcta(e, (int)q);
				n = toSphere(Point2s(qn[0] + d[0], qn[0] + d[1]), (int)q);
				count++;
			} else {//no correction
				for(int k = 0; k < 3; k++)
					n[k] = (int16_t)(e[k]*32767.0f);
			}
		}
	}
}

NormalAttr.prototype.toSphere = function(i, input, out, unit) {

	var t = this;
	var j = i*2;
	var k = i*3;
	var av0 = input[j] > 0? input[j]:-input[j];
	var av1 = input[j+1] > 0? input[j+1]:-input[j+1];
	out[k] = input[j];
	out[k+1] = input[j+1];
	out[k+2] = unit - av0 - av1;
	if (out[k+2] < 0) {
		out[k] = (input[j] > 0)? unit - av1 : av1 - unit;
		out[k+1] = (input[j+1] > 0)? unit - av0: av0 - unit;
	}
	var len = 1/Math.sqrt(out[k]*out[k] + out[k+1]*out[k+1] + out[k+2]*out[k+2]);
	if(t.type == t.Type.INT16) {
		len *= 32767;
	}
	out[k] *= len;
	out[k+1] *= len;
	out[k+2] *= len;
}

NormalAttr.prototype.toOcta = function(i, input, output, unit) {
	var k = i*2;
	var j = i*3; //input

	var av0 = input[j] > 0? input[j]:-input[j];
	var av1 = input[j+1] > 0? input[j+1]:-input[j+1];
	var av2 = input[j+2] > 0? input[j+2]:-input[j+2];
	var len = av0 + av1 + av2;
	var p0 = av0/len;
	var p1 = av1/len;

	var ap0 = p0 > 0? p0: -p0;
	var ap1 = p1 > 0? p1: -p1;

	if(input[j+2] < 0) {
		p0 = (input[j] < 0)? 1.0 - ap1 : ap1 - 1;
		p1 = (input[j] < 0)? 1.0 - ap0 : ap0 - 1;
	}
	output[k] = p0*unit;
	output[k+1] = p1*unit;

/*
		Point2f p(v[0], v[1]);
		p /= (fabs(v[0]) + fabs(v[1]) + fabs(v[2]));

		if(v[2] < 0) {
			p = Point2f(1.0f - fabs(p[1]), 1.0f - fabs(p[0]));
			if(v[0] < 0) p[0] = -p[0];
			if(v[1] < 0) p[1] = -p[1];
		}
		return Point2i(p[0]*unit, p[1]*unit); */
}


/* INDEX ATTRIBUTE */

function IndexAttr(nvert, nface, type) {
	var t = this;
	if(!type || type == Type.UINT32) //uint32 by default
		t.faces = new Uint32Array(nface*3);
	else if(type == Type.UINT16)
		t.faces = new Uint16Buffer(nface*3);
	else
		throw "Unsupported type";

	t.prediction = new Uint32Array(nvert*3);
}

IndexAttr.prototype = {
decode: function(stream) {
	var t = this;

	var max_front = stream.readInt();
	t.front = new Int32Array(max_front*5);
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

