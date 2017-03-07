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

#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <deque>

#include "tunstall.h"
#include "nxzdecoder.h"

using namespace std;
using namespace nx;

class DEdge2 { //decompression edges
public:
	uint32_t v0, v1, v2; //used for parallelogram prediction
	uint32_t prev, next;
	bool deleted;
	DEdge2(uint32_t a = 0, uint32_t b = 0, uint32_t c = 0, uint32_t p = 0, uint32_t n = 0):
		v0(a), v1(b), v2(c), prev(p), next(n), deleted(false) {}
};

NxzDecoder::NxzDecoder(int len, uchar *input):
	short_normals(false), short_index(false), vertex_count(0) {

	stream.init(len, input);

	flags = stream.read<int>();
	nvert = stream.read<int>();
	entropy = (Entropy)stream.read<uchar>();

	coord.q = stream.read<char>();
	coord.o[0] = stream.read<int>();
	coord.o[1] = stream.read<int>();
	coord.o[2] = stream.read<int>();

	if(flags & NORMAL) {
		normals_prediction = (Normals)stream.read<uchar>();
		norm.q = stream.read<float>();
	}

	if(flags & COLOR)
		for(int k = 0; k < 4; k++)
			color[k].q = stream.read<float>();

	if(flags & UV) {
		uv.q = stream.read<float>();
		uv.o[0] = stream.read<int>();
		uv.o[1] = stream.read<int>();
	}
	int d = DATA;
	while(flags & d) {
		float q = stream.read<float>();
		int o = stream.read<int>();
		data.push_back(Attribute<int>(q, o));
		d *= 2;
	}
}

void NxzDecoder::setCoords(float *buffer) { coord.buffer = buffer; }

void NxzDecoder::setNormals(float *buffer) { norm.buffer = buffer; }


void NxzDecoder::setNormals(int16_t *buffer) {
	norm.buffer = buffer;
	short_normals = true;
}

void NxzDecoder::setColors(uchar *buffer) { color[0].buffer = buffer; }

void NxzDecoder::setUvs(float *buffer) { uv.buffer = buffer; }

void NxzDecoder::setData(int pos, float *buffer) { data[pos].buffer = buffer; }

void NxzDecoder::setIndex(uint32_t *buffer) { face.buffer = buffer; }

void NxzDecoder::setIndex(int16_t *buffer) {
	face.buffer = buffer;
	short_index = true;
}

void NxzDecoder::decode() {
	if(flags & INDEX)
		decodeMesh();
	else
		decodePointCloud();

}

void NxzDecoder::decodePointCloud() {
	decodeZPoints();

	if(flags & NORMAL)
		decodeNormals();

	if(flags & COLOR)
		decodeColors();

	if(flags & UV)
		decodeUvs();

	if(flags & DATA)
		decodeDatas();
}

void NxzDecoder::decodeZPoints() {
	std::vector<uchar> diffs;
	Tunstall::decompress(stream, diffs);

	BitStream bitstream;
	stream.read(bitstream);

	vector<ZPoint> zpoints(nvert);
	bitstream.read(63, zpoints[0].bits);
	//TODO Optimize: we could do it without zpoints array
	for(size_t i = 1; i < zpoints.size(); i++) {
		ZPoint &p = zpoints[i];
		p = zpoints[i-1];
		uchar d = diffs[i-1];
		p.setBit(d, 1);
		uint64_t e = 0; //1<<d;
		bitstream.read(d, e);
		p.bits &= ~((1LL<<d) -1);
		p.bits |= e;
	}

	Point3f *coords = (Point3f *)coord.buffer;
	for(size_t i = 0; i < zpoints.size(); i++)
		coords[i] = zpoints[i].toPoint(coord.o, coord.q);
}

void NxzDecoder::decodeCoords() {
	std::vector<uchar> diffs;
	Tunstall::decompress(stream, diffs);

	BitStream bitstream;
	stream.read(bitstream);

	Point3i *coords = (Point3i *)coord.buffer;
	for(uint32_t i = 0; i < nvert; i++)
		decodeDiff(diffs[i], bitstream, coords[i]);
}

void NxzDecoder::decodeNormals() {
	std::vector<uchar> diffs;
	Tunstall::decompress(stream, diffs);

	BitStream bitstream;
	stream.read(bitstream);

	int count =0 ;
	if(short_normals) {
		Point3s *normals = (Point3s *)norm.buffer;
		for(uint32_t i = 0; i < nvert; i++)
			if(normals_prediction != BORDER || boundary[i])
				decodeDiff(diffs[count++], bitstream, normals[i]);
	} else {
		Point3i *normals = (Point3i *)norm.buffer;
		for(uint32_t i = 0; i < nvert; i++)
			if(normals_prediction != BORDER || boundary[i])
				decodeDiff(diffs[count++], bitstream, normals[i]);
	}
}

void NxzDecoder::decodeColors() {

	Color4b *colors = (Color4b *)color[0].buffer;

	vector<uchar> diffs;
	BitStream bitstream;

	for(int k = 0; k < 4; k++) {
		Tunstall::decompress(stream, diffs);
		stream.read(bitstream);

		for(uint32_t i = 0; i < nvert; i++)
			colors[i][k] = (char)decodeDiff(diffs[i], bitstream);
	}
}

void NxzDecoder::decodeUvs() {
	std::vector<uchar> diffs;
	Tunstall::decompress(stream, diffs);

	BitStream bitstream;
	stream.read(bitstream);

	Point2i *uvs = (Point2i *)uv.buffer;
	for(uint32_t i = 0; i < nvert; i++)
		decodeDiff(diffs[i], bitstream, uvs[i]);
}

void NxzDecoder::decodeDatas() {
	for(auto &da: data) {
		std::vector<uchar> diffs;
		Tunstall::decompress(stream, diffs);

		BitStream bitstream;
		stream.read(bitstream);

		int *d = (int *)da.buffer;
		for(uint32_t i = 0; i < nvert; i++)
			d[i] = decodeDiff(diffs[i], bitstream);
	}
}


void NxzDecoder::dequantize() {
	if(flags & INDEX) {
		Point3i *coords = (Point3i *)coord.buffer;
		Point3f *points = (Point3f *)coords;
		for(uint32_t i = 0; i < nvert; i++) {
			Point3i &q = coords[i];
			Point3f &p = points[i];
			p[0] = (q[0] + coord.o[0])*coord.q;
			p[1] = (q[1] + coord.o[1])*coord.q;
			p[2] = (q[2] + coord.o[2])*coord.q;
		}
	}
	if(flags & NORMAL) {
	}
}

template <class F>
void integrateNormals(int nvert, int nface, F *index, Point3f *coords, Point3f *normals) {
	for(int i = 0; i < nface; i++) {
		Point3f &p0 = coords[index[0]];
		Point3f &p1 = coords[index[1]];
		Point3f &p2 = coords[index[2]];
		Point3f n = (p1 - p0) ^ (p2 - p0);
		normals[index[0]] += n;
		normals[index[1]] += n;
		normals[index[2]] += n;
		index += 3;
	}
}

void NxzDecoder::computeNormals(Point3f *normals3f) {

	Point3f *coords = (Point3f *)coord.buffer;

	if(short_index)
		integrateNormals<uint16_t>(nvert, nface, (uint16_t *)face.buffer, coords, normals3f);
	else
		integrateNormals<uint32_t>(nvert, nface, (uint32_t *)face.buffer, coords, normals3f);

	for(unsigned int i = 0; i < nvert; i++) {
		Point3f &n = normals3f[i];
		n /= n.norm();
	}
}

void NxzDecoder::computeNormals(Point3s *normals3s) {
	vector<Point3f> normals3f(nvert, Point3f(0, 0, 0));
	computeNormals(&*normals3f.begin());

	for(unsigned int i = 0; i < nvert; i++) {
		Point3f &nf = normals3f[i];
		Point3s &ns = normals3s[i];
		for(int k = 0; k < 3; k++)
			ns[k] = (int16_t)(nf[k]*32767.0f);
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

void NxzDecoder::decodeMesh() {
	last.reserve(nvert);

	nface = stream.read<int>();
	groups.resize(stream.read<int>());
	for(uint32_t &g: groups)
		g = stream.read<int>();
	Tunstall tunstall;
	tunstall.decompress(stream, clers);
	BitStream bitstream;
	stream.read(bitstream);

	decodeCoords();

	if(flags & NORMAL)
		decodeNormals();

	if(flags & COLOR)
		decodeColors();

	if(flags & UV)
		decodeUvs();

	decodeDatas();

	uint32_t start = 0;
	uint32_t cler = 0; //keeps track of current cler
	for(uint32_t &end: groups) {
		decodeFaces(bitstream, start*3, end*3, cler);
		start = end;
	}

	dequantize();
}

void NxzDecoder::decodeFaces(BitStream &bitstream, uint32_t start, uint32_t end, uint32_t &cler) {

	uint16_t *faces16 = ((uint16_t *)face.buffer);
	uint32_t *faces32 = ((uint32_t *)face.buffer);

	bool hasUv = (flags & UV) != 0;

	Point3i *coords = (Point3i *)coord.buffer;
	Point2i *uvs = (Point2i *)uv.buffer;

	vector<uchar> clers;

	//unsigned int current = 0; //keep track of connected component start

	vector<int> delayed;
	deque<int> faceorder;
	vector<DEdge2> front;


	while(start < end) {
		if(!faceorder.size() && !delayed.size()) {

			Point3i estimated(0, 0, 0);
			Point2i texestimated(0, 0);

			int last_index = -1;
			int index[3];
			for(int k = 0; k < 3; k++) {
				last.push_back(last_index);
				int v = decodeVertex(estimated, texestimated, last_index);
				index[k] = v;
				if(short_index)
					faces16[start++] = v;
				else
					faces32[start++] = v;

				estimated = coords[v];
				if(hasUv)
					texestimated = uvs[v];
				last_index = v;
			}
			int current_edge = front.size();
			for(int k = 0; k < 3; k++) {
				faceorder.push_back(front.size());
				front.push_back(DEdge2(index[next_(k)], index[prev_(k)], index[k], current_edge + prev_(k), current_edge + next_(k)));
			}
			continue;
		}

		int f;
		if(faceorder.size()) {
			f = faceorder.front();
			faceorder.pop_front();
		} else {
			f = delayed.back();
			delayed.pop_back();
		}
		DEdge2 &e = front[f];
		if(e.deleted) continue;
		e.deleted = true;
		int prev = e.prev;
		int next = e.next;
		int v0 = e.v0;
		int v1 = e.v1;

		DEdge2 &previous_edge = front[prev];
		DEdge2 &next_edge = front[next];

		int c = clers[cler++];
		if(c == BOUNDARY) continue;

		int first_edge = front.size();
		int opposite = -1;

		if(c == VERTEX || c == SPLIT) {

			if(c == SPLIT) {
				opposite = bitstream.readUint(32);
			} else {
				//parallelogram prediction.
				int v2 = e.v2;
				Point3i predicted = coords[v0] + coords[v1] - coords[v2];
				Point2i texpredicted(0, 0);
				if(hasUv)
					texpredicted = uvs[v0] + uvs[v1] - uvs[v2];

				//not v0!, in encode we use face, here edge which is inverted
				opposite = decodeVertex(predicted, texpredicted, v1);
//				if(diff != 0)
//					last.push_back(v1);
			}
			previous_edge.next = first_edge;
			next_edge.prev = first_edge + 1;
			faceorder.push_front(front.size());
			front.push_back(DEdge2(v0, opposite, v1, prev, first_edge + 1));
			faceorder.push_back(front.size());
			front.push_back(DEdge2(opposite, v1, v0, first_edge, next));

		} else if(c == END) {
			previous_edge.deleted = true;
			next_edge.deleted = true;
			front[previous_edge.prev].next = next_edge.next;
			front[next_edge.next].prev = previous_edge.prev;
			opposite = previous_edge.v0;

		} else if(c == LEFT) {
			previous_edge.deleted = true;
			front[previous_edge.prev].next = first_edge;
			front[next].prev = first_edge;
			opposite = previous_edge.v0;
			faceorder.push_front(front.size());
			front.push_back(DEdge2(opposite, v1, v0, previous_edge.prev, next));

		} else if(c == RIGHT) {
			next_edge.deleted = true;
			front[next_edge.next].prev = first_edge;
			front[prev].next = first_edge;
			opposite = next_edge.v1;
			faceorder.push_front(front.size());
			front.push_back(DEdge2(v0, opposite, v1, prev, next_edge.next));

		} else if(c == DELAY) {
			e.deleted = false;
			delayed.push_back(f);
			continue;
		} else {
			assert(0);
		}
		if(short_index) {
			faces16[start++] = v1;
			faces16[start++] = v0;
			faces16[start++] = opposite;
		} else {
			faces32[start++] = v1;
			faces32[start++] = v0;
			faces32[start++] = opposite;
		}
	}
}

int NxzDecoder::decodeVertex(const Point3i &predicted, const Point2i &texpredicted,
							 int last_index) {
	int v = vertex_count++;

	Point3i &p = ((Point3i *)coord.buffer)[v];
	p += predicted;

	if(flags & UV) {
		//TODO these needds to be precomputed.
		Point2i &u = (((Point2i *)&uv.buffer)[v]);
		u += texpredicted;
	}

	if(last_index <= 0)
		return v;

	if(flags & NORMAL) {
		if(normals_prediction == DIFF) {
			if(short_normals) {
				Point3s &n = ((Point3s *)norm.buffer)[v];
				Point3s &ref = ((Point3s *)norm.buffer)[last_index];
				n[0] += ref[0];
				n[1] += ref[1];
				if(n[0] < -coord.q) n[0] += coord.q;
				else if(n[0] > +coord.q) n[0] -= coord.q;
				if(n[1] < -coord.q) n[1] += coord.q;
				else if(n[1] > +coord.q) n[1] -= coord.q;
			} else {
				Point3i &n = ((Point3i *)norm.buffer)[v];
				Point3i &ref = ((Point3i *)norm.buffer)[last_index];
				n[0] += ref[0];
				n[1] += ref[1];
				if(n[0] < -coord.q) n[0] += coord.q;
				else if(n[0] > +coord.q) n[0] -= coord.q;
				if(n[1] < -coord.q) n[1] += coord.q;
				else if(n[1] > +coord.q) n[1] -= coord.q;
			}
		}
	}
	if(flags & COLOR) {
		Color4b &c = ((Color4b *)color[0].buffer)[v];
		c += ((Color4b *)color[0].buffer)[last_index];
	}

	for(auto &da: data) {
		int &d = ((int *)da.buffer)[v];
		d += ((int *)da.buffer)[last_index];
	}

	return v;
}




int NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream) {
	if(diff == 0)
		return 0;
/*
	//TODO could we use the max system to make it quicker?
	uint64_t c = 1<<(diff);
	bitstream.read(diff, c);
	int val = (int)c;

	if(val & 0x1)
		val >>= 1;
	else
		val = -(val >> 1);
	return val;
	*/

	int val = bitstream.readUint(diff);
	int middle = (1<<diff)>>1;
	if(val < middle)
		val = -val -middle;
	return val;
}


void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point3i &p) {
	//assert(diff < 22);

	//int mask = (1<<diff)-1;
	int max = (1<<(diff-1));
	p[0] = bitstream.readUint(diff) - max;
	p[1] = bitstream.readUint(diff) - max;
	p[2] = bitstream.readUint(diff) - max;
	//TODO is this really faster than just do multiple readings?

	/*	bitstream.read(diff*3, bits);
	for(int k = 2; k >= 0; k--) {
		p[k] = (bits & mask) - max;
		bits >>= diff;
	} */
}

void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point3s &p) {
	//assert(diff < 22);

	//int mask = (1<<diff)-1;
	int max = (1<<(diff-1));
	p[0] = bitstream.readUint(diff) - max;
	p[1] = bitstream.readUint(diff) - max;
	p[2] = bitstream.readUint(diff) - max;
	//TODO is this really faster than just do multiple readings?

	/*	bitstream.read(diff*3, bits);
	for(int k = 2; k >= 0; k--) {
		p[k] = (bits & mask) - max;
		bits >>= diff;
	} */
}


void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point2i &p) {
	//assert(diff < 22);

	//int mask = (1<<diff)-1;
	int max = (1<<(diff-1));
	p[0] = bitstream.readUint(diff) - max;
	p[1] = bitstream.readUint(diff) - max;
	//TODO is this really faster than just do multiple readings?

	/*	bitstream.read(diff*3, bits);
	for(int k = 2; k >= 0; k--) {
		p[k] = (bits & mask) - max;
		bits >>= diff;
	} */
}

Point2i NxzDecoder::encodeNormal(Point3f v, int unit) {
	Point2f p(v[0], v[1]);
	p /= (abs(v[0]) + abs(v[1]) + abs(v[2]));

	if(v[2] < 0) {
		p = Point2f(1.0f - abs(p[1]), 1.0f - abs(p[0]));
		if(v[0] < 0) p[0] = -p[0];
		if(v[1] < 0) p[1] = -p[1];
	}
	return Point2i(p[0]*unit, p[1]*unit);
}


Point3f NxzDecoder::decodeNormal3i(Point2i v, int unit) {
	Point3f n(v[0], v[1], unit - abs(v[0]) -abs(v[1]));
	if (n[2] < 0) {
		n[0] = ((v[0] > 0)? 1 : -1)*(unit - abs(v[1]));
		n[1] = ((v[1] > 0)? 1 : -1)*(unit - abs(v[0]));
	}
	n /= n.norm();
	return n;
}

Point3s NxzDecoder::decodeNormal3s(Point2s v, int unit) {
	Point3f n(v[0], v[1], unit - abs(v[0]) -abs(v[1]));
	if (n[2] < 0) {
		n[0] = ((v[0] > 0)? 1 : -1)*(unit - abs(v[1]));
		n[1] = ((v[1] > 0)? 1 : -1)*(unit - abs(v[0]));
	}
	n /= n.norm();
	return Point3s(n[0]*32767, n[1]*32767, n[2]*32767);
}
