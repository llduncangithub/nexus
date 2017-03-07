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
#include "nxzencoder.h"

using namespace nx;
using namespace std;

static int ilog2(uint64_t p) {
	int k = 0;
	while ( p>>=1 ) { ++k; }
	return k;
}


NxzEncoder::NxzEncoder(uint32_t _nvert, uint32_t _nface, Entropy en):
	flags(0), nvert(_nvert), nface(_nface),
	coord(0.0f, Point3i(0)),
	norm(0.0f, Point2i(0)),
	uv(0.0f, Point2i(0)),
	header_size(0), entropy(en), normals_prediction(ESTIMATED), current_vertex(0) {
	/*	coord_o(2147483647), uv_o(2147483647), coord_q(0), uv_q(0), coord_bits(12), uv_bits(12), norm_bits(10),
	coord_size(0), normal_size(0), color_size(0), face_size(0), uv_size(0) {
	color_bits[0] = color_bits[1] = color_bits[2] = color_bits[3] = 6; */
	color[0] = Attribute<uchar>(0.0f, 0);
	color[1] = Attribute<uchar>(0.0f, 0);
	color[2] = Attribute<uchar>(0.0f, 0);
	color[3] = Attribute<uchar>(0.0f, 0);

	coord.values.resize(nvert);
	face.values.resize(nface*3);
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

	coord.q = q;

	//	Point3i cmax(-2147483647);
	for(uint32_t i = 0; i < nvert; i++) {
		Point3i &q = coord.values[i];
		Point3f &p = coords[i];
		for(int k = 0; k < 3; k++)
			q[k] = floor((p[k] - o[k] + coord.q/2.0f)/coord.q);
	}
	setCoordBits();
}

/* if not q specified use 1/10th of average leght of edge  */
void NxzEncoder::addCoords(float *buffer, uint32_t *_index, float q, Point3f o) {
	flags |= INDEX;
	face.values.resize(nface*3);
	memcpy(&*face.values.begin(),_index,  nface*12);

	flags |= COORD;
	Point3f *coords = (Point3f *)buffer;
	if(q == 0) {
		double average = 0; //for precision when using large number of faces
		for(uint32_t f = 0; f < nface*3; f += 3)
			average += (coords[_index[f]] - coords[_index[f+1]]).norm();
		q = (float)(average/nface)/20.0f;
	}
	coord.diffs.resize(nvert);
	addCoords(buffer, q, o);
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
	memcpy(&*coord.values.begin(), buffer, nvert*12);
	setCoordBits();
}

void NxzEncoder::addCoords(int *buffer, uint32_t *_index) {
	coord.diffs.resize(nvert);
	addCoords(buffer);

	flags |= INDEX;
	memcpy(&*face.values.begin(),_index,  nface*12);
}

void NxzEncoder::addCoords(int *buffer, uint16_t *_index) {
	vector<uint32_t> tmp(nface*3);
	for(uint32_t i = 0; i < nface*3; i++)
		tmp[i] = _index[i];
	addCoords(buffer, &*tmp.begin());
}

void NxzEncoder::setCoordBits() {
	for(auto &q: coord.values)
		coord.o.setMin(q);

	if(1) {
		Point3i cmax(-2147483647);
		for(auto &q: coord.values) {
			q -= coord.o;
			cmax.setMax(q);
		}

		int bits = 1+std::max(std::max(ilog2(cmax[0]), ilog2(cmax[1])), ilog2(cmax[2]));
		if(!(flags & INDEX) && bits >= 22)
			cerr << "Quantiziation above 22 bits for point clouds is not supported..." << endl;
		else
			cout << "Quantization coords in " << bits << " bits\n";
	}
}


/* NORMALS */


void NxzEncoder::addNormals(float *buffer, int bits, Normals no) {
	normals_prediction = no;
	flags |= NORMAL;
	norm.values.resize(nvert);
	norm.diffs.resize(nvert);
	norm.q = pow(2, bits-1);

	Point3f *normals = (Point3f *)buffer;
	for(uint32_t i = 0; i < nvert; i++)
		norm.values[i] = encodeNormal(normals[i], (int)norm.q);
}

void NxzEncoder::addNormals(int16_t *buffer, int bits, Normals no) {
	normals_prediction = no;
	assert(bits <= 17);
	vector<float> tmp(nvert*3);
	for(uint32_t i = 0; i < nvert*3; i++)
		tmp[i] = buffer[i]/32767.0f;
	addNormals(&*tmp.begin(), bits, no);
}

/* COLORS */

void NxzEncoder::addColors(unsigned char *buffer, int lumabits, int chromabits, int alphabits) {
	flags |= COLOR;
	int q[4];
	q[0] = color[0].q = 1<<(8 - lumabits);
	q[1] = q[2] = color[1].q = color[2].q = 1<<(8 - chromabits);
	q[3] = color[3].q = 1<<(8-alphabits);
	for(int k = 0; k < 4; k++) {
		color[k].values.resize(nvert);
		color[k].diffs.resize(nvert);
	}

	Color4b *colors = (Color4b *)buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		Color4b c = colors[i].toYCC();
		for(int k = 0; k < 4; k++)
			color[k].values[i] = (int)floor((c[k] + q[k]/2.0f)/q[k]);
	}
}

void NxzEncoder::addUV(float *buffer, float q) {
	flags |= UV;
	if(q == 0) q = pow(2, -12);
	uv.q = q;
	uv.values.resize(nvert);
	uv.diffs.resize(nvert);
	Point2f *uvs = (Point2f *)buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		uv.values[i][0] = (int)round(uvs[i][0]/uv.q);
		uv.values[i][1] = (int)round(uvs[i][1]/uv.q);
	}

	for(auto &q: uv.values)
		uv.o.setMin(q);

	if(1) {
		Point2i cmax(-2147483647);
		for(auto &q: uv.values) {
			q -= uv.o;
			cmax.setMax(q);
		}
		int bits = 1 + std::max(ilog2(cmax[0]), ilog2(cmax[1]));
		cout << "Texture cooridnates quantization in " << bits << " bits\n";
	}
}

void NxzEncoder::addData(float *buffer, float q, float offset) {
	flags |= DATA*(1<<data.size());
	data.resize(data.size()+1);
	auto &da = data.back();
	da.q = q;
	da.values.resize(nvert);
	da.diffs.resize(nvert);
	for(uint32_t i = 0; i < nvert; i++)
		da.values[i] = (int)round((buffer[i] - offset)/q);
	setDataBits();
}

void NxzEncoder::addData(int *buffer) {
	flags |= DATA*(1<<data.size());
	data.resize(data.size()+1);
	auto &d = data.back();
	d.q = 1.0f;
	d.values.resize(nvert);
	memcpy(&*d.values.begin(), buffer, nvert*4);
	setDataBits();
}

void NxzEncoder::setDataBits() {
	auto &d = data.back();
	d.o = 2147483647;
	for(auto &q: d.values)
		if(q < d.o) d.o = q;
}

void NxzEncoder::encode() {
	stream.reserve(nvert);

	stream.write<int>(flags);
	stream.write<int>(nvert);
	stream.write<char>(entropy);

	stream.write<float>(coord.q);
	stream.write<int>(coord.o[0]);
	stream.write<int>(coord.o[1]);
	stream.write<int>(coord.o[2]);

	if(flags & NORMAL) {
		stream.write<uchar>(normals_prediction);
		stream.write<float>(norm.q);
	}

	if(flags & COLOR)
		for(int k = 0; k < 4; k++)
			stream.write<float>(color[k].q);


	if(flags & UV) {
		stream.write<float>(uv.q);
		stream.write<int>(uv.o[0]);
		stream.write<int>(uv.o[1]);
	}

	for(auto &d: data) {
		stream.write<float>(d.q);
		stream.write<int>(d.o);
	}

	if(flags & INDEX)
		encodeMesh();
	else
		encodePointCloud();
}

void NxzEncoder::encodePointCloud() {
	std::vector<ZPoint> zpoints(nvert);

	for(uint32_t i = 0; i < nvert; i++) {
		Point3i &q = coord.values[i];
		zpoints[i] = ZPoint(q[0], q[1], q[2], 21, i);
	}
	sort(zpoints.rbegin(), zpoints.rend());//, greater<ZPoint>());

	int count = 0;
	for(unsigned int i = 1; i < nvert; i++) {
		if(zpoints[i] == zpoints[count])
			continue;
		count++;
		zpoints[count] = zpoints[i];
	}
	count++;
	nvert = count;
	zpoints.resize(nvert);

	header_size = stream.elapsed();

	encodeZPoints(zpoints);

	//reorder attributes, compute differences end encode
	if(flags & NORMAL) {
		norm.diffs[0] = norm.values[zpoints[0].pos];
		for(uint32_t i = 1; i < nvert; i++)
			norm.diffs[i] = norm.values[zpoints[i].pos] - norm.values[zpoints[i-1].pos];

		encodeNormals();
	}

	if(flags & COLOR) {
		for(int k = 0; k < 4; k++) {
			color[k].diffs[0] = color[k].values[zpoints[0].pos];
			for(uint32_t i = 1; i < nvert; i++)
				color[k].diffs[i] = color[k].values[zpoints[i].pos] - color[k].values[zpoints[i-1].pos];
		}
		encodeColors();
	}

	if(flags & UV) {
		uv.diffs[0] = uv.values[zpoints[0].pos];
		for(uint32_t i = 1; i < nvert; i++)
			uv.diffs[i] = uv.values[zpoints[i].pos] - uv.values[zpoints[i-1].pos];
		encodeUvs();
	}

	for(auto &da: data) {
		vector<int> tmp(nvert);
		da.diffs[0] = da.values[zpoints[0].pos];
		for(uint32_t i = 1; i < nvert; i++)
			da.diffs[i] = da.values[zpoints[i].pos] - da.values[zpoints[i-1].pos];
		encodeDatas();
	}
}

void NxzEncoder::encodeZPoints(std::vector<ZPoint> &zpoints) {
	vector<uchar> diffs;
	BitStream bitstream(nvert/2);
	bitstream.write(zpoints[0].bits, 63);

	for(uint32_t pos = 1; pos < nvert; pos++) {
		ZPoint &p = zpoints[pos-1]; //previous point
		ZPoint &q = zpoints[pos]; //current point
		uchar d = p.difference(q);
		diffs.push_back(d);
		//we can get away with d diff bits (should be d+1) because the first bit will always be 1 (sorted q> p)
		bitstream.write(q.bits, d); //rmember to add a 1, since
	}

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	stream.write(bitstream);
	coord.size = stream.elapsed();
}

void NxzEncoder::encodeCoords() {
	vector<uchar> diffs;
	BitStream bitstream(nvert/2);

	for(Point3i &d: coord.diffs)
		encodeDiff(diffs, bitstream, d);

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	stream.write(bitstream);
	coord.size = stream.elapsed();
}

void NxzEncoder::encodeNormals() {
	vector<uchar> diffs;
	vector<uchar> signs;
	BitStream bitstream(nvert/64);

	//TOO: VS using Point2i diff
	for(uint32_t i = 0; i < nvert; i++)
		if(normals_prediction != BORDER || boundary[i])
			encodeDiff(diffs, bitstream, norm.diffs[i]);

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	stream.write(bitstream);

	norm.size = stream.elapsed();
}

void NxzEncoder::encodeColors() {

	//TODO: Slightly more efficient to use only one bitstream,
	//TODO again: worth separating components?

	for(int k = 0; k < 4; k++) {
		BitStream bitstream(nvert/2);
		vector<uchar> diffs;

		//use complement to 255 (hence char conversion)
		for(uchar &c: color[k].diffs)
			encodeDiff(diffs, bitstream, (char)c);

		Tunstall::compress(stream, &*diffs.begin(), diffs.size());
		stream.write(bitstream);
		color[k].size = stream.elapsed();
	}
}


void NxzEncoder::encodeUvs() {
	vector<uchar> diffs;
	BitStream bitstream(nvert/2);

	for(Point2i &u: uv.diffs)
		encodeDiff(diffs, bitstream, u);

	Tunstall::compress(stream, &*diffs.begin(), diffs.size());
	stream.write(bitstream);
	uv.size = stream.elapsed();
}

void NxzEncoder::encodeDatas() {
	for(auto &da: data) {
		BitStream bitstream(nvert/2);
		vector<uchar> diffs;

		for(int &v: da.diffs)
			encodeDiff(diffs, bitstream, v);

		Tunstall::compress(stream, &*diffs.begin(), diffs.size());
		stream.write(bitstream);
		da.size += stream.elapsed();
	}
}



//compact in place faces in data, update patches information, compute topology and encode each patch.
void NxzEncoder::encodeMesh() {
	encoded.resize(nvert, -1);

	if(!groups.size()) groups.push_back(nface);
	//remove degenerate faces
	uint32_t start =  0;
	uint32_t count = 0;
	for(uint32_t &end: groups) {
		for(uint32_t i = start; i < end; i++) {
			int *f = &face.values[i*3];

			if(f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
				continue;

			if(count != i) {
				int *dest = &face.values[count*3];
				dest[0] = f[0];
				dest[1] = f[1];
				dest[2] = f[2];
			}
			count++;
		}
		start = end;
		end = count;
	}
	face.values.resize(count*3);
	nface = count;

	stream.write<int>(nface);
	stream.write<int>(groups.size());
	for(uint32_t &end: groups)
		stream.write<int>(end);

	header_size = stream.elapsed();

	if((flags & NORMAL) && normals_prediction != DIFF) {
		computeNormals(norm.diffs);
		if(normals_prediction == BORDER)
			markBoundary(); //mark boundary points on original vertices.
		for(uint32_t i = 0; i < nvert; i++)
			norm.values[i] -= norm.diffs[i];
	}


	BitStream bitstream(nface);
	start =  0;
	for(uint32_t &end: groups) {
		encodeFaces(bitstream, start, end);
		start = end;
	}

	Tunstall::compress(stream, &*clers.begin(), clers.size());
	stream.write(bitstream);
	face.size = stream.elapsed();

	encodeCoords();
	if(flags & NORMAL)
		encodeNormals();

	if(flags & COLOR)
		encodeColors();

	if(flags & UV)
		encodeUvs();

	encodeDatas();

}

//How to determine if a vertex is a boundary without topology:
//for each edge a vertex is in, add or subtract the id of the other vertex depending on order
//for internal vertices sum is zero.
//unless we have non manifold strange configurations zero wont happen.

void NxzEncoder::markBoundary() {
	boundary.resize(nvert, false);

	vector<int> count(nvert, 0);
	for(uint32_t i = 0; i < nface; i++) {
		int *f = &face.values[i*3];
		count[f[0]] += (int)f[1] - (int)f[2];
		count[f[1]] += (int)f[2] - (int)f[0];
		count[f[2]] += (int)f[0] - (int)f[1];
	}
	for(uint32_t i = 0; i < nvert; i++)
		if(count[i] != 0)
			boundary[i] = true;
}


void NxzEncoder::computeNormals(vector<Point2i> &estimated) {
	vector<Point3f> tmp(nvert, Point3f(0, 0, 0));
	estimated.resize(nvert);
	for(uint32_t i = 0; i < nface; i++) {
		int *f = &face.values[i*3];
		Point3i &p0 = coord.values[f[0]];
		Point3i &p1 = coord.values[f[1]];
		Point3i &p2 = coord.values[f[2]];
		Point3i qn = (( p1 - p0) ^ (p2 - p0));
		Point3f n(qn[0], qn[1], qn[2]);
		tmp[f[0]] += n;
		tmp[f[1]] += n;
		tmp[f[2]] += n;
	}
	//normalize
	for(uint32_t i = 0; i < nvert; i++) {
		Point3f &n = tmp[i];
		n /= n.norm();
		estimated[i] = encodeNormal(n, norm.q);
	}
}

class McFace {
public:
	uint32_t f[3];
	uint32_t t[3]; //topology: opposite face
	uint32_t i[3]; //index in the opposite face of this face: faces[f.t[k]].t[f.i[k]] = f;
	McFace(uint32_t v0 = 0, uint32_t v1 = 0, uint32_t v2 = 0) {
		f[0] = v0; f[1] = v1; f[2] = v2;
		t[0] = t[1] = t[2] = 0xffffffff;
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
		int * f = &face.values[i*3];
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
				//				last.push_back(last_index);

				encodeVertex(index, coord_estimated, uv_estimated, bitstream, last_index);

				last_index = index;
				coord_estimated = coord.values[index];
				if(flags & UV)
					uv_estimated = uv.values[index];

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
		assert(enext < (int)front.size());
		assert(eprev < (int)front.size());
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
			if(encoded[opposite] == -1)
				clers.push_back(VERTEX);
			else
				clers.push_back(SPLIT);

			//compute how much would it take to save vertex information:
			//we need to estimate opposite direction using v0 + v1 -
			int v2 = faces[e.face].f[e.side];

			Point3i coord_predicted = coord.values[v0] + coord.values[v1] - coord.values[v2];
			Point2i uv_predicted(0, 0);
			if(flags & UV)
				uv_predicted = uv.values[v0] + uv.values[v1] - uv.values[v2];


			if(encoded[opposite] == -1) {//only first time
				//				last.push_back(v0);
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
}
//TODO faster version?
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

	assert(target < (int)nvert);
	if(encoded[target] != -1) { //write index of target.
		uint64_t bits = encoded[target];
		bitstream.write(bits, 32); //this should be related to nvert, actually.
		return;
	}

	//notice how vertex needs to be reordered
	coord.diffs[current_vertex] = coord.values[target] - predicted;

	if((flags & NORMAL)) {
		Point2i &dt = norm.diffs[target];
		dt = norm.values[target];
		if(normals_prediction == DIFF) {
			if(last >= 0)
				dt -= norm.values[last];
			if(dt[0] < -coord.q) dt[0] += coord.q;
			else if(dt[0] > +coord.q) dt[0] -= coord.q;
			if(dt[1] < -coord.q) dt[1] += coord.q;
			else if(dt[1] > +coord.q) dt[1] -= coord.q;
		}
	}

	if(flags & COLOR)
		for(int k = 0; k < 4; k++) {
			uchar tar = color[k].values[target];
			uchar las = ((last < 0)? 0: color[k].values[last]);
//			if(k == 2)
//				cout << "k: " << k << " last: " << last << " t: " << (int)tar << " l: " << (int)las << endl;
			color[k].diffs[current_vertex] = color[k].values[target] - ((last < 0)? 0: color[k].values[last]);
			//			if(k == 1 && last >= 0)
			//				cout << (int)(char)color[k].diffs[target] << " ";
		}

	if(flags & UV)
		uv.diffs[current_vertex] = uv.values[target] - texpredicted; //this is the estimated point

	for(auto &da: data)
		da.diffs[current_vertex] = da.values[target] - (last < 0)? 0: da.values[last];

	encoded[target] = current_vertex++;
}


/*
 0    -> 0 []
-1, 1 -> 1, [0,1]
-2, 2
-3, 3 -> 2 [0, 3]   (1<<2) max;
-4, 4
..
-7, 7     3 [0,7] if(val>>1 < (1<<(diff-1)) val = -(1<<diff) + 1 + val;

var middle = (1<<(diff-1));
if(val < middle) val = -val -middle;

//encode:
0 - > 0
			   000  001 010 011
-4   4         100  101 110 111
-7   7;        111

if(d == 0) d = 0;
if(d < 0) n = -d; else n = d;
diff = ilog2(n)
var middle = 1<<(diff-1);
if(d < 0)
	d = -middle - val;

*/

//val can be zero.
void NxzEncoder::encodeDiff(vector<uchar> &diffs, BitStream &stream, int val) {
	/*
	val = Tunstall::toUint(val)+1;
	int ret = ilog2(val);
	diffs.push_back(ret);
	if(ret > 0)
		stream.write(val, ret); */

	int ref = val;

	if(val == 0) {
		diffs.push_back(0);
//		cout << "Encode diff: " << ref << " to: " << 0 << " bits and " << 0 << " as bits " << endl;

		return;
	}
	int ret = ilog2(abs(val)) + 1;  //0 -> 0, [1,-1] -> 1 [-2,-3,2,3] -> 2
	diffs.push_back(ret);
	int middle = (1<<ret)>>1;
	if(val < 0) val = -val -middle;
	stream.writeUint(val, ret);

	//cout << "Encode diff: " << ref << " to: " << ret << " bits and " << val << " as bits " << endl;
}

void NxzEncoder::encodeDiff(vector<uchar> &diffs, BitStream &bitstream, const Point2i &val) {
	int diff = std::max(needed(val[0]), needed(val[1]));
	assert(diff > 0);
	diffs.push_back(diff);

	int max = 1<<(diff-1);
	bitstream.writeUint(val[0] + max, diff);
	bitstream.writeUint(val[1] + max, diff);
}

void NxzEncoder::encodeDiff(vector<uchar> &diffs, BitStream &bitstream, const Point3i &val) {
	int diff = std::max(needed(val[0]), needed(val[1]));
	diff = std::max(diff, needed(val[2]));
	assert(diff > 0);
	diffs.push_back(diff);

	int max = 1<<(diff-1);
	bitstream.writeUint(val[0] + max, diff);
	bitstream.writeUint(val[1] + max, diff);
	bitstream.writeUint(val[2] + max, diff);
}


//encode between -unit and +unit;
Point2i NxzEncoder::encodeNormal(Point3f v, int unit) {
	Point2f p(v[0], v[1]);
	p /= (abs(v[0]) + abs(v[1]) + abs(v[2]));

	if(v[2] < 0) {
		p = Point2f(1.0f - abs(p[1]), 1.0f - abs(p[0]));
		if(v[0] < 0) p[0] = -p[0];
		if(v[1] < 0) p[1] = -p[1];
	}
	return Point2i(p[0]*unit, p[1]*unit);
}
