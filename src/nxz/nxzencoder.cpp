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

#include <assert.h>
#include <deque>
#include <algorithm>

//#include <QTime>

#include "tunstall.h"
#include "../nxszip/fpu_precision.h"
#include "nxzencoder.h"

using namespace nx;
using namespace std;

static int ilog2(uint64_t p) {
	int k = 0;
	while ( p>>=1 ) { ++k; }
	return k;
}


NxzEncoder::NxzEncoder(uint32_t _nvert, uint32_t _nface):
	flags(0), nvert(_nvert), nface(_nface),
	coord_o(2147483647), uv_o(2147483647), coord_q(0), uv_q(0), coord_bits(12), uv_bits(12), norm_bits(10),
	coord_size(0), normal_size(0), color_size(0), face_size(0), uv_size(0) {

	color_bits[0] = color_bits[1] = color_bits[2] = color_bits[3] = 6;
	qcoords.resize(nvert);
	index.resize(nface*3);
}


/* TODO: euristic for point clouds should be: about 1/10 of nearest neighbor.
  for now just use volume and number of points
  could also quantize to find a side where the points halve and use 1/10 */

void NxzEncoder::addCoords(float *buffer, float q, Point3f o) {
	flags |= COORD;
	Point3f *coords = (Point3f *)buffer;

	if(q == 0) {
		Point3f max(-FLT_MAX), min(FLT_MAX);
		for(uint32_t i = 0; i < nvert; i++) {
			min.setMin(coords[i]);
			max.setMax(coords[i]);
		}
		max -= min;
		q = pow(max[0]*max[1]*max[2], 2.0/3.0)/nvert;
		if(o == Point3f(FLT_MAX))
			o = min;
	}
	if(o == Point3f(FLT_MAX)) {
		for(uint32_t i = 0; i < nvert; i++)
			o.setMin(coords[i]);
	}

	coord_q = q;

	Point3i cmax(-2147483647);
	for(uint32_t i = 0; i < nvert; i++) {
		Point3i &q = qcoords[i];
		Point3f &p = coords[i];
		for(int k = 0; k < 3; k++)
			q[k] = floor((p[k] - o[k] + coord_q/2.0f)/coord_q);
	}
	setCoordBits();
}

/* if not q specified use 1/10th of average leght of edge  */
void NxzEncoder::addCoords(float *buffer, uint32_t *_index, float q, Point3f o) {
	flags |= COORD;
	Point3f *coords = (Point3f *)buffer;
	if(q == 0) {
		double average = 0; //for precision when using large number of faces
		for(uint32_t f = 0; f < nface*3; f += 3)
			average += (coords[index[f]] - coords[index[f+1]]).norm();
		q = (float)(average/nface)/10.0f;
	}
	addCoords(buffer, q, o);

	flags |= INDEX;
	index.resize(nface*3);
	memcpy(_index, &*index.begin(), nface*12);
}

void NxzEncoder::addCoords(float *buffer, uint16_t *_index, float q, Point3f o) {
	vector<uint32_t> tmp(nface*3);
	for(uint32_t i = 0; i < nvert*3; i++)
		tmp[i] = _index[i];
	addCoords(buffer, &*tmp.begin(), q, o);
}


//you do the quantization step
void NxzEncoder::addCoords(int *buffer) {
	flags |= COORD;
	memcpy(buffer, &*qcoords.begin(), nvert*12);
	setCoordBits();
}

void NxzEncoder::addCoords(int *buffer, uint32_t *_index) {
	addCoords(buffer);

	flags |= INDEX;
	memcpy(_index, &*index.begin(), nface*12);
}

void NxzEncoder::addCoords(int *buffer, uint16_t *_index) {
	vector<uint32_t> tmp(nface*3);
	for(uint32_t i = 0; i < nface*3; i++)
		tmp[i] = _index[i];
	addCoords(buffer, &*tmp.begin());
}

void NxzEncoder::setCoordBits() {
	for(auto &q: qcoords)
		coord_o.setMin(q);

	Point3i cmax(-2147483647);
	for(auto &q: qcoords) {
		q -= coord_o;
		cmax.setMax(q);
	}

	coord_bits = 1+std::max(std::max(ilog2(cmax[0]), ilog2(cmax[1])), ilog2(cmax[2]));
	assert(coord_bits < 22);
}


/* NORMALS */


void NxzEncoder::addNormals(float *buffer, int bits) {
	flags |= NORMAL;
	qnormals.resize(nvert);
	norm_bits = bits;
	norm_q = pow(2, -bits+1);
	int qmax = 1<<(bits-1);

	Point3f *normals = (Point3f *)buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		Point3i &qn = qnormals[i];
		Point3f &n = normals[i];
		for(int k = 0; k < 3; k++) {
			qn[k] = (int)std::round(n[k]/norm_q);
			if(qn[k] > qmax) qn[k] = qmax;
			if(qn[k] < -qmax) qn[k] = -qmax;
		}
	}
}

void NxzEncoder::addNormals(int16_t *buffer, int bits) {
	assert(bits <= 17);
	vector<float> tmp(nvert*3);
	for(uint32_t i = 0; i < nvert*3; i++)
		tmp[i] = buffer[i]/32767.0f;
	addNormals(&*tmp.begin(), bits);
}

/* COLORS */

void NxzEncoder::addColors(unsigned char *buffer, int lumabits, int chromabits, int alphabits) {
	flags |= COLOR;
	color_bits[0] = lumabits;
	color_bits[1] = color_bits[2] = chromabits;
	color_bits[3] = alphabits;
	for(int k = 0; k < 4; k++)
		color_q[k] = pow(1, -color_bits[k]);

	qcolors.resize(nvert);
	memcpy(buffer, &*qcolors.begin(), nvert*4);
	for(auto &c: qcolors) {
		c = c.toYCC();
		for(int k = 0; k < 4; k++)
			c[k] = (int)floor((c[k] + color_q[k]/2.0f)/color_q[k]);
	}
}

void NxzEncoder::addUV(float *buffer, float q) {
	flags |= UV;
	if(q == 0) q = pow(2, -12);
	uv_q = q;
	quvs.resize(nvert);
	Point3f *uvs = (Point3f *)buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		quvs[i][0] = (int)round(uvs[i][0]/uv_q);
		quvs[i][1] = (int)round(uvs[i][1]/uv_q);
	}

	for(auto &q: quvs)
		uv_o.setMin(q);

	Point2i cmax(-2147483647);
	for(auto &q: quvs) {
		q -= uv_o;
		cmax.setMax(q);
	}
	uv_bits = 1 + std::max(ilog2(cmax[0]), ilog2(cmax[1]));
}

void NxzEncoder::addData(float *buffer, float q, float offset) {
	flags |= DATA + qdatas.size();
	qdatas.resize(qdatas.size()+1);
	auto &qdata = qdatas.back();
	qdata.resize(nvert);
	for(uint32_t i = 0; i < nvert; i++)
		qdata[i] = (int)round((buffer[i] - offset)/q);
	data_q.push_back(q);

	setDataBits();
}

void NxzEncoder::addData(int *buffer) {
	flags |= DATA + qdatas.size();
	qdatas.resize(qdatas.size()+1);
	auto &qdata = qdatas.back();
	qdata.resize(nvert);
	memcpy(buffer, &*qdata.begin(), nvert*4);
	data_q.push_back(1.0f);

	setDataBits();
}

void NxzEncoder::setDataBits() {
	int min = 2147483647;
	for(auto &q: qdatas.back())
		if(q < min) min = q;
	data_o.push_back(min);

	int cmax(-2147483647);
	for(auto &q: qdatas.back()) {
		q -= min;
		if(q > cmax) cmax = q;
	}
	data_bits.push_back(1 + ilog2(cmax));
}



void NxzEncoder::encode() {
	FpuPrecision::store();
	FpuPrecision::setFloat();

	stream.reserve(nvert);

	stream.write<int>(flags);

	if(flags & NORMAL)
		stream.write<char>(norm_bits);

	if(flags & COLOR)
		for(int k = 0; k < 4; k++)
			stream.write<char>(color_bits[k]);


	if(flags & UV) {
		stream.write<int>(uv_o[0]);
		stream.write<int>(uv_o[1]);
		stream.write<float>(uv_q);
		stream.write<char>(uv_bits);
	}

	for(uint32_t i = 0; i < qdatas.size(); i++) {
		stream.write<int>(data_o[i]);
		stream.write<float>(data_q[i]);
		stream.write<char>(data_bits[i]);
	}

	if(nface)
		encodeMesh();
	else
		encodePointCloud();

	/*
	cout << "Coord bpv: " << coord_size*8/(float)node.nvert << endl;
	cout << "Face bpv: " << face_size*8/(float)node.nvert << endl;
	cout << "Norm bpv: " << normal_size*8/(float)node.nvert << endl;
	cout << "Color bpv: " << color_size*8/(float)node.nvert << endl; */

	FpuPrecision::restore();
}

void NxzEncoder::encodePointCloud() {
	zpoints.resize(nvert);
	for(uint32_t i = 0; i < nvert; i++) {
		Point3i &q = qcoords[i];
		zpoints[i] = ZPoint(q[0], q[1], q[2], coord_bits, i);
	}
	sort(zpoints.rbegin(), zpoints.rend());//, greater<ZPoint>());

	//reorder attributes
	if(flags & NORMAL) {
		vector<Point3i> tmp(nvert);
		for(uint32_t i = 0; i < nvert; i++)
			tmp[i] = qnormals[zpoints[i].pos];
		swap(tmp, qnormals);
	}
	if(flags & COLOR) {
		vector<Color4b> tmp(nvert);
		for(uint32_t i = 0; i < nvert; i++)
			tmp[i] = qcolors[zpoints[i].pos];
		swap(tmp, qcolors);
	}
	if(flags & UV) {
		vector<Point2i> tmp(nvert);
		for(uint32_t i = 0; i < nvert; i++)
			tmp[i] = quvs[zpoints[i].pos];
		swap(tmp, quvs);
	}
	for(auto &qdata: qdatas) {
		vector<int> tmp(nvert);
		for(uint32_t i = 0; i < nvert; i++)
			tmp[i] = qdata[zpoints[i].pos];
		swap(tmp, qdata);
	}

	//we need to remove duplicated vertices
	int count = 0;
	for(unsigned int i = 1; i < nvert; i++) {
		if(zpoints[i] == zpoints[count])
			continue;
		count++;
		zpoints[count] = zpoints[i];

		if(qnormals.size())      qnormals[count] = qnormals[i];
		if(qcolors.size())       qcolors[count] = qcolors[i];
		if(quvs.size())          quvs[count] = quvs[i];
		for(auto &qdata: qdatas) qdata[count] = qdata[i];
	}
	count++;

	zpoints.resize(count);
	if(qnormals.size())      qnormals.resize(count);
	if(qcolors.size())       qcolors.resize(count);
	if(quvs.size())          quvs.resize(count);
	for(auto &qdata: qdatas) qdata.resize(count);

	nvert = count;
	stream.write<int>(nvert);

	encodeCoords();
	if(flags & NORMAL) encodeNormals();
	if(flags & COLOR)  encodeColors();
	if(flags & UV)     encodeUvs();
	if(flags & DATA)   encodeDatas();
}

void NxzEncoder::encodeCoords() {
	vector<uchar> diffs;
	BitStream bitstream(nvert/2);
	bitstream.write(zpoints[0].bits, coord_bits*3);

	for(uint32_t pos = 1; pos < nvert; pos++) {
		ZPoint &p = zpoints[pos-1]; //previous point
		ZPoint &q = zpoints[pos]; //current point
		uchar d = p.difference(q);
		diffs.push_back(d);
		//we can get away with d diff bits (should be d+1) because the first bit will always be 1 (sorted q> p)
		bitstream.write(q.bits, d); //rmember to add a 1, since
	}

	int start = stream.size();

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	stream.write(bitstream);

	coord_size = stream.size() - start;
}

void NxzEncoder::encodeNormals() {
	//TODO is it worth having diffs for component?
	vector<uchar> diffs;
	vector<uchar> signs;
	BitStream bitstream(nvert/64);

	for(int k = 0; k < 2; k++) {
		int old = qnormals[0][k];
		bitstream.write(Tunstall::toUint(old), norm_bits);
		for(unsigned int i = 1; i < qnormals.size(); i++) {
			int d = qnormals[i][k] -old;
			encodeDiff(diffs, bitstream, d);
			old = qnormals[i][k];
		}
	}
	for(uint32_t i = 0; i < qnormals.size(); i++) {
		bool signbit = qnormals[i][2] > 0;
		signs.push_back(signbit);
	}
	int start = stream.size();

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	//TODO worth compressing, dont thinks so.?
	Tunstall::compress(stream, &*signs.begin(), signs.size());
	stream.write(bitstream);

	normal_size = stream.size() - start;
}

void NxzEncoder::encodeColors() {
	BitStream bitstream(nvert/2);

	int start = stream.size();
	for(int k = 0; k < 4; k++) {
		vector<uchar> diffs;

		int old = qcolors[0][k];
		bitstream.write(old, color_bits[k]);
		for(uint32_t i = 1; i < nvert; i++) {
			int d = qcolors[i][k] - old;
			encodeDiff(diffs, bitstream, d);
			old = qcolors[i][k];
		}
		Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	}
	stream.write(bitstream);
	color_size = stream.size() - start;
}


void NxzEncoder::encodeUvs() {
	BitStream bitstream(nvert/2);

	int start = stream.size();
	for(int k = 0; k < 2; k++) {
		vector<uchar> diffs;
		int old = quvs[0][k];
		bitstream.write(old, uv_bits);
		for(uint32_t i = 1; i < nvert; i++) {
			int d = quvs[i][k] - old;
			encodeDiff(diffs, bitstream, d);
			old = quvs[i][k];
		}
		Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	}
	stream.write(bitstream);
	uv_size = stream.size() - start;
}

void NxzEncoder::encodeDatas() {
	for(uint32_t k = 0; k < qdatas.size(); k++) {
		auto &qdata = qdatas[k];
		BitStream bitstream(nvert/2);
		vector<uchar> diffs;

		int start = stream.size();
		int old = qdata[0];
		bitstream.write(old, data_bits[k]);
		for(uint32_t i = 1; i < nvert; i++) {
			int d = qdata[i] - old;
			encodeDiff(diffs, bitstream, d);
			old = qdata[i];
		}
		Tunstall::compress(stream, &*diffs.begin(), diffs.size());

		stream.write(bitstream);
		data_size.push_back(stream.size() - start);
	}
}



//compact in place faces in data, update patches information, compute topology and encode each patch.
void NxzEncoder::encodeMesh() {
	//TODO what if nface == 0?

	encoded.resize(nvert, -1);
	last.reserve(nvert);
	order.reserve(nvert);

	if(!groups.size()) groups.push_back(nface);
	//remove degenerate faces
	uint32_t start =  0;
	uint32_t count = 0;
	for(uint32_t &end: groups) {
		for(uint32_t i = start; i < end; i++) {
			uint32_t *face = &index[i*3];

			if(face[0] == face[1] || face[0] == face[2] || face[1] == face[2])
				continue;

			if(count != i) {
				uint32_t *dest = &index[count*3];
				dest[0] = face[0];
				dest[1] = face[1];
				dest[2] = face[2];
			}
			count++;
		}
		start = end;
		end = count;
	}
	index.resize(count*3);
	nface = count;

	stream.write<int>(nvert);
	stream.write<int>(nface);

	index_bits = ilog2(nface);

	if(flags & NORMAL) {
		markBoundary(); //mark boundary points on original vertices.
		vector<Point3i> estimated(nvert);
		computeNormals(estimated);
		for(uint32_t i = 0; i < nvert; i++) {
			qnormals[i][0] -= estimated[i][0];
			qnormals[i][1] -= estimated[i][1];
			qnormals[i][2] = (estimated[i][2]*qnormals[i][2] < 0)?1:0;
		}
	}

	BitStream bitstream(nface);
	start =  0;
	for(uint32_t &end: groups) {
		encodeFaces(bitstream, start, end);
		start = end;
	}
	Tunstall::compress(stream, &*clers.begin(), clers.size());
	Tunstall::compress(stream, &*dcoords.begin(), nvert);
	if(flags & NORMAL) Tunstall::compress(stream, &*dnormals.begin(), nvert);
	if(flags & COLOR)
		for(int k = 0; k < 4; k++)
			Tunstall::compress(stream, &*dcolors[k].begin(), nvert);

	if(flags & UV)     Tunstall::compress(stream, &*duvs.begin(), nvert);
	for(auto &ddata: ddatas)
		Tunstall::compress(stream, &*ddata.begin(), nvert);
	stream.write(bitstream);
}

/*
void NxzEncoder::encodeMeshNormals() {

	vector<uchar> diffs;
	vector<uchar> signs;
	BitStream bitstream(nvert/64);

	markBoundary(); //mark boundary points on original vertices.
	vector<Point3i> estimated_normals(nvert);
	computeNormals(estimated_normals);

	for(unsigned int i = 0; i < nvert; i++) {
		int pos = order[i];
		if(!boundary[pos]) continue;
		Point3i &computed = estimated_normals[pos];
		Point3i &original = qnormals[pos];
		for(int comp = 0; comp < 2; comp++) {
			int d = (int)(original[comp]/side - computed[comp]/side); //act1ual value - predicted
			encodeDiff(diffs, bitstream, d);
		}
		bool signbit = (computed[2]*original[2] < 0);
		signs.push_back(signbit);
	}

	int start = stream.size();

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	Tunstall::compress(stream, &*signs.begin(), signs.size());
	stream.write(bitstream);

	normal_size = stream.size() - start;
}
 */

//how to determine if a vertex is a boundary without topology:
//for each edge a vertex is in, add or subtract the id of the other vertex depending on order
//for internal vertices sum is zero.
//unless we have strange configurations and a lot of sfiga, zero wont happen. //TODO think about this
void NxzEncoder::markBoundary() {
	boundary.resize(nvert, false);

	vector<int> count(nvert, 0);
	for(uint32_t i = 0; i < nface; i++) {
		uint32_t *f = &index[i*3];
		count[f[0]] += (int)f[1] - (int)f[2];
		count[f[1]] += (int)f[2] - (int)f[0];
		count[f[2]] += (int)f[0] - (int)f[1];
	}
	for(uint32_t i = 0; i < nvert; i++)
		if(count[i] != 0)
			boundary[i] = true;
}


void NxzEncoder::computeNormals(vector<Point3i> &estimated) {
	estimated.resize(nvert, Point3i(0, 0, 0));

	for(uint32_t i = 0; i < nface; i++) {
		uint32_t *face = &index[i*3];
		Point3i &p0 = qcoords[face[0]];
		Point3i &p1 = qcoords[face[1]];
		Point3i &p2 = qcoords[face[2]];
		Point3i n = (( p1 - p0) ^ (p2 - p0));
		estimated[face[0]] += n;
		estimated[face[1]] += n;
		estimated[face[2]] += n;
	}
	//normalize
	for(Point3i &n: estimated) {
		float norm = sqrt(float((float)n[0]*n[0] + (float)n[1]*n[1] + (float)n[2]*n[2]));
		for(int k = 0; k < 2; k++)
			n[k] = (int)std::round((n[k]/norm)/norm_q);
	}
}

/*
void NxzEncoder::encodeColors() {

	BitStream bitstream(node.nvert/2);
	vector<uchar> diffs[4];

	int steps[4];
	for(int k = 0; k < 4; k++)
		steps[k] = (1<<(8 - color_bits[k]));


		vector<Color4b> qcolors(nvert);
		for(unsigned int i = 0; i < nvert; i++) {
			for(int k = 0; k < 4; k++) {
				qcolors[i][k] = colors[i][k]; ///steps[k]*steps[k]; AARRRGRGGH
			}
			qcolors[i] = toYCC(qcolors[i]);
		}

		for(int i = 0; i < nvert; i++) {
			Color4b &c = qcolors[order[i]];
			int l = last[i];
			Color4b b;
			if(l < 0) //Is this 0, 0, 0 useful? Don't think so.
				b = Color4b(0, 0, 0, 0);
			else {
				b = qcolors[last[i]];
			}

			for(int k = 0; k < 4; k++) {
				int d = c[k]/steps[k] - b[k]/steps[k];
				encodeDiff(diffs[k], bitstream, d);
			}
		}



	int start = stream.size();

	for(int k = 0; k < 4; k++)
		stream.write<char>(color_bits[k]);

	for(int k = 0; k < 4; k++)
		Tunstall::compress(stream, &*diffs[k].begin(), diffs[k].size());

	stream.write(bitstream);

	color_size = stream.size() - start;
} */

class McFace {
public:
	uint32_t f[3];
	uint32_t t[3]; //topology: opposite face
	uint32_t i[3]; //index in the opposite face of this face: faces[f.t[k]].t[f.i[k]] = f;
	McFace(uint32_t v0 = 0, uint32_t v1 = 0, uint32_t v2 = 0) {
		f[0] = v0; f[1] = v1; f[2] = v2;
		t[0] = t[1] = t[2] = 0xffff;
	}
	bool operator<(const McFace &face) const {
		if(f[0] < face.f[0]) return true;
		if(f[0] > face.f[0]) return false;
		if(f[1] < face.f[1]) return true;
		if(f[1] > face.f[1]) return false;
		return f[2] < face.f[2];
	}
	bool operator>(const McFace &face) const {
		if(f[0] > face.f[0]) return true;
		if(f[0] < face.f[0]) return false;
		if(f[1] > face.f[1]) return true;
		if(f[1] < face.f[1]) return false;
		return f[2] > face.f[2];
	}
};

class CEdge { //compression edges
public:
	uint32_t face;
	uint32_t side; //opposited to side vertex of face
	uint32_t prev, next;
	bool deleted;
	CEdge(uint32_t f = 0, uint32_t s = 0, uint32_t p = 0, uint32_t n = 0):
		face(f), side(s), prev(p), next(n), deleted(false) {}
};


class McEdge { //topology edges
public:
	uint32_t face, side;
	uint32_t v0, v1;
	bool inverted;
	//McEdge(): inverted(false) {}
	McEdge(uint32_t _face, uint32_t _side, uint32_t _v0, uint32_t _v1): face(_face), side(_side), inverted(false) {

		if(_v0 < _v1) {
			v0 = _v0; v1 = _v1;
			inverted = false;
		} else {
			v1 = _v0; v0 = _v1;
			inverted = true;
		}
	}
	bool operator<(const McEdge &edge) const {
		if(v0 < edge.v0) return true;
		if(v0 > edge.v0) return false;
		return v1 < edge.v1;
	}
	bool match(const McEdge &edge) {
		if(inverted && edge.inverted) return false;
		if(!inverted && !edge.inverted) return false;
		return v0 == edge.v0 && v1 == edge.v1;
	}
};

static void buildTopology(vector<McFace> &faces) {
	//create topology;
	vector<McEdge> edges;
	for(unsigned int i = 0; i < faces.size(); i++) {
		McFace &face = faces[i];
		for(int k = 0; k < 3; k++) {
			int kk = k+1;
			if(kk == 3) kk = 0;
			int kkk = kk+1;
			if(kkk == 3) kkk = 0;
			edges.push_back(McEdge(i, k, face.f[kk], face.f[kkk]));
		}
	}
	sort(edges.begin(), edges.end());

	McEdge previous(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff);
	for(unsigned int i = 0; i < edges.size(); i++) {
		McEdge &edge = edges[i];
		if(edge.match(previous)) {
			uint32_t &edge_side_face = faces[edge.face].t[edge.side];
			uint32_t &previous_side_face = faces[previous.face].t[previous.side];
			if(edge_side_face == 0xffffffff && previous_side_face == 0xffffffff) {
				edge_side_face = previous.face;
				faces[edge.face].i[edge.side] = previous.side;
				previous_side_face = edge.face;
				faces[previous.face].i[previous.side] = edge.side;
			}
		} else
			previous = edge;
	}
}
static int next_(int t) {
	t++;
	if(t == 3) t = 0;
	return t;
}
static int prev_(int t) {
	t--;
	if(t == -1) t = 2;
	return t;
}

void NxzEncoder::encodeFaces(BitStream &bitstream, int start, int end) {
	int useless = 0;

	vector<McFace> faces(end - start);
	for(int i = start; i < end; i++) {
		uint32_t * f = &index[i*3];
		faces[i - start] = McFace(f[0], f[1], f[2]);
		assert(f[0] != f[1] && f[1] != f[2] && f[2] != f[0]);
	}

	buildTopology(faces);

	unsigned int current = 0;          //keep track of connected component start

	vector<int> delayed;
	deque<int> faceorder;
	vector<CEdge> front;

	vector<bool> visited(faces.size(), false);
	unsigned int totfaces = faces.size();

	//	vector<int> test_faces;
	int counting = 0;
	while(totfaces > 0) {
		if(!faceorder.size() && !delayed.size()) {

			while(current != faces.size()) {   //find first triangle non visited
				if(!visited[current]) break;
				current++;
			}
			if(current == faces.size()) break; //no more faces to encode exiting

			//encode first face: 3 vertices indexes, and add edges
			unsigned int current_edge = front.size();
			McFace &face = faces[current];
			Point3i coord_estimated(0, 0, 0);
			Point2i uv_estimated(0, 0);
			int last_index = -1;
			for(int k = 0; k < 3; k++) {
				int index = face.f[k];
				last.push_back(last_index);

				encodeVertex(index, coord_estimated, uv_estimated, bitstream, last_index);

				last_index = index;
				coord_estimated = qcoords[index];
				if(flags & UV)
					uv_estimated = qtexcoords[index];

				faceorder.push_back(front.size());
				front.push_back(CEdge(current, k, current_edge + prev_(k), current_edge + next_(k)));
			}
			counting++;
			visited[current] = true;
			current++;
			totfaces--;
			continue;
		}
		int c;
		if(faceorder.size()) {
			c = faceorder.front();
			faceorder.pop_front();
		} else {
			c = delayed.back();
			delayed.pop_back();
		}
		CEdge &e = front[c];
		if(e.deleted) continue;
		e.deleted = true;

		//opposite face is the triangle we are encoding
		uint32_t opposite_face = faces[e.face].t[e.side];
		int opposite_side = faces[e.face].i[e.side];

		if(opposite_face == 0xffffffff || visited[opposite_face]) { //boundary edge or glue
			clers.push_back(BOUNDARY);
			continue;
		}

		assert(opposite_face < faces.size());
		McFace &face = faces[opposite_face];

		int k2 = opposite_side;
		int k0 = next_(k2);
		int k1 = next_(k0);

		//check for closure on previous or next edge
		int eprev = e.prev;
		int enext = e.next;
		assert(enext < front.size());
		assert(eprev < front.size());
		CEdge &previous_edge = front[eprev];
		CEdge &next_edge = front[enext];

		int first_edge = front.size();
		bool close_left = (faces[previous_edge.face].t[previous_edge.side] == opposite_face);
		bool close_right = (faces[next_edge.face].t[next_edge.side] == opposite_face);

		if(close_left && close_right) {
			clers.push_back(END);
			previous_edge.deleted = true;
			next_edge.deleted = true;
			front[previous_edge.prev].next = next_edge.next;
			front[next_edge.next].prev = previous_edge.prev;

		} else if(close_left) {
			clers.push_back(LEFT);
			previous_edge.deleted = true;
			front[previous_edge.prev].next = first_edge;
			front[enext].prev = first_edge;
			faceorder.push_front(front.size());
			front.push_back(CEdge(opposite_face, k1, previous_edge.prev, enext));

		} else if(close_right) {
			clers.push_back(RIGHT);
			next_edge.deleted = true;
			front[next_edge.next].prev = first_edge;
			front[eprev].next = first_edge;
			faceorder.push_front(front.size());
			front.push_back(CEdge(opposite_face, k0, eprev, next_edge.next));

		} else {
			int v0 = face.f[k0];
			int v1 = face.f[k1];
			int opposite = face.f[k2];

			if(encoded[opposite] != -1 && faceorder.size()) { //split, but we can still delay it.
				e.deleted = false; //undelete it.
				delayed.push_back(c);
				clers.push_back(DELAY);
				continue;
			}

			clers.push_back(VERTEX);
			//compute how much would it take to save vertex information:
			//we need to estimate opposite direction using v0 + v1 -
			int v2 = faces[e.face].f[e.side];

			Point3i coord_predicted = qcoords[v0] + qcoords[v1] - qcoords[v2];
			Point2i uv_predicted(0, 0);
			if(flags & UV)
				uv_predicted = qtexcoords[v0] + qtexcoords[v1] - qtexcoords[v2];


			if(encoded[opposite] == -1) {//only first time
				last.push_back(v0);
			} else
				useless++;  //we encoutered it already.

			encodeVertex(opposite, coord_predicted, uv_predicted, bitstream, v0);


			previous_edge.next = first_edge;
			next_edge.prev = first_edge + 1;
			faceorder.push_front(front.size());
			front.push_back(CEdge(opposite_face, k0, eprev, first_edge+1));
			faceorder.push_back(front.size());
			front.push_back(CEdge(opposite_face, k1, first_edge, enext));
		}

		counting++;
		assert(!visited[opposite_face]);
		visited[opposite_face] = true;
		totfaces--;
	}
	/*	vector<int> freq(6, 0);
	for(int i: clers)
		freq[i]++;
	cout << "Clers.size: " << clers.size() << " face.size: " << faces.size() << endl;

	for(int i: freq)
		cout << "I: " << 100.0*i/clers.size() << endl;

	//save obj
	QFile file("test.obj");
	file.open(QFile::WriteOnly);
	QTextStream tstream(&file);
	for(Point3i p: qpoints)
		tstream << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";

	for(int i = 0; i < test_faces.size(); i += 3) {
		int *f = &test_faces[i];
		tstream << "f " << (f[0]+1) << " " << (f[1] +1) << " " << (f[2]+1) << "\n";
	}
*/

/*

	int stream_start = stream.size();
	Tunstall::compress(stream, &*clers.begin(), clers.size());

	face_size += stream.size() - stream_start;
	stream_start = stream.size();

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	if(sig.vertex.hasTextures()) {
		int bits =0;
		for(auto &n: tdiffs) {
			bits += 2*n;
			//			cout << "N: " << (int)n << endl;
		}
		Tunstall::compress(stream, &*tdiffs.begin(), tdiffs.size());



	}
	stream.write(bitstream);

	coord_size += stream.size() - stream_start; */
}

int needed(int v) {
	int n = 1;
	while(1) {
		if(v >= 0) {
			if(v < (1<<(n-1)))
				break;
		} else
			if(v + (1<<(n-1)) >= 0)
				break;
		n++;
	}
	return n;
}

void NxzEncoder::encodeVertex(int target, const Point3i &predicted, const Point2i &texpredicted, BitStream &bitstream, int last) {

	//compute how much would it take to save vertex information:
	//we need to estimate opposite direction using v0 + v1 -

	assert(target < nvert);
	if(encoded[target] != -1) { //write index of target.
		dcoords.push_back(0);
		uint64_t bits = encoded[target];
		bitstream.write(bits, index_bits); //this should be related to nvert, actually.
		return;
	}
	assert(order.size() < nvert);
	encoded[target] = order.size();
	order.push_back(target);

	{
		Point3i d = qcoords[target] - predicted; //this is the estimated point.
		int diff = 0;
		for(int k = 0; k < 3; k++) {
			int n = needed(d[k]);
			if(n > diff)
				diff = n;
		}
		assert(diff > 0);
		for(int k = 0; k < 3; k++)
			d[k] += 1<<(diff-1); //bring it into the positive.

		dcoords.push_back(diff);
		bitstream.write(d[0], diff);
		bitstream.write(d[1], diff);
		bitstream.write(d[2], diff);
	}

	if(flags & NORMAL) {
		Point3i &d = qnormals[target];
		if(boundary[target]) {
			int diff = 0;
			for(int k = 0; k < 2; k++) {
				int n = needed(d[k]);
				if(n > diff)
					diff = n;
			}
			assert(diff > 0);
			for(int k = 0; k < 2; k++)
				d[k] += 1<<(diff-1); //bring it into the positive.
			dnormals.push_back(diff);
			bitstream.write(d[0], diff);
			bitstream.write(d[1], diff);
			bitstream.write(d[2], 1);
		}
	}
	if(flags & COLOR) {
		Color4b &d = qcolors[target];
		if(last >= 0)
			d -= qcolors[last];
		for(int k = 0; k < 4; k++)
			encodeDiff(dcolors[k], bitstream, d[k]);
	}

	if(flags & UV) {
		Point2i dt = qtexcoords[target] - texpredicted; //this is the estimated point

		int tdiff = 0;
		for(int k = 0; k < 2; k++) {
			int n = needed(dt[k]);
			if(n > tdiff)
				tdiff = n;
			if(tdiff >= 22) {
				cerr << "Texture precision required cannot be bigger than 2^-21. \n";
			}
		}
		for(int k = 0; k < 2; k++)
			dt[k] += 1<<(tdiff-1);

		duvs.push_back(tdiff);
		bitstream.write(dt[0], tdiff);
		bitstream.write(dt[1], tdiff);
	}

	for(int k = 0; k < qdatas.size(); k++) {
		auto &qdata = qdatas[k];
		int d = qdata[target];
		if(last >= 0)
			d -= qdata[last];
		encodeDiff(ddatas[k], bitstream, d);
	}
}

//val can be zero.
void NxzEncoder::encodeDiff(vector<uchar> &diffs, BitStream &stream, int val) {
	val = Tunstall::toUint(val)+1;
	int ret = ilog2(val);
	diffs.push_back(ret);
	if(ret > 0)
		stream.write(val, ret);
}

