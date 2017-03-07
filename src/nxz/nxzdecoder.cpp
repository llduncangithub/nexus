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
	vertex_count(0), short_normals(false), short_index(false) {

	stream.init(len, input);

	flags = stream.read<int>();
	nvert = stream.read<int>();
	entropy = (Entropy)stream.read<uchar>();

	coord.q = stream.read<char>();
	coord.o[0] = stream.read<int>();
	coord.o[1] = stream.read<int>();
	coord.o[2] = stream.read<int>();

	if(flags & NORMAL) {
		normals = (Normals)stream.read<uchar>();
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

void NxzDecoder::decodePointClout() {
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
		decodeDiff(diffs, bitstream, coords[i]);
}

void NxzDecoder::decodeNormals() {
	std::vector<uchar> diffs;
	Tunstall::decompress(stream, diffs);

	BitStream bitstream;
	stream.read(bitstream);

	for(uint32_t i = 0; i < nvert; i++)
		if(normals != BORDER || boundary[i])
			encodeDiff(diffs, bitstream, norm.diffs[i]);

	Point3s *norms_s = (Point3s *)norm.buffer;
	Point3i *norms_i = (Point3i *)norm.buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		Point2i n;
		decodeDiff(diffs, bitstream, n);
		float x = n[0];
		float y = n[1];
		float z = 1.0 * 1.0 - n[0]*n[0] - n[1]*n[1];
		if(z < 0) z = 0;
		else z = sqrt(z);
		if(z > 1.0) z = 1.0;
		if(!signs[i]) z = -z;
		if(short_normals) {
			norms[i][0] = (int16_t)x*32767.0f;
			norms[i][1] = (int16_t)y*32767.0f;
			norms[i][2] = (int16_t)z*32767.0f;
		} else
			norms[i] = Point3f(x, y, z);
	}

	/*	markBoundary();
		Point3s *estimated_normals = data.normals(sig, node.nvert);
		computeNormals(norm.values);
		if(sig.vertex.hasTextures()) //hack, actually fixing normals makes it worse if textures are enabled.
			return;

		size_t diffcount = 0;
		int signcount = 0;
		//gert difference between original and predicted
		for(unsigned int i = 0; i < node.nvert; i++) {
			vcg::Point3s &N = estimated_normals[i];
			if(!boundary[i]) continue;
			if(diffcount >= diffs.size()) break; //TODO assert here.
			for(int comp = 0; comp < 2; comp++) {
				int n = N[comp]/side;
				assert(diffcount < diffs.size());
				int diff = decodeDiff(diffs[diffcount++], bitstream);
				N[comp] = (n + diff)*side;
			}
			float x = N[0];
			float y = N[1];
			float z = 32767.0f*32767.0f - x*x - y*y;

			if(z < 0) z = 0;
			z = sqrt(z);
			//sign
			if(z > 32767.0f) z = 32767.0f;
			bool signbit = signs[signcount++];
			if((N[2] < 0) != signbit)
				z = -z;
			N[2] = (int16_t)z;
		}
	} */
}



void NxzDecoder::decodeColors() {

	Color4b *colors = data.colors(sig, node.nvert);

	for(int k = 0; k < 4; k++)
		color_q[k] = stream.read<char>();

	std::vector<uchar> diffs[4];
	for(int k = 0; k < 4; k++) {
		Tunstall tunstall;
		tunstall.decompress(stream, diffs[k]);
	}
	BitStream bitstream;
	stream.read(bitstream);

	bool indexed = sig.face.hasIndex();
	colors[0] = Color4b(0, 0, 0, 0);

	int count = 0;

	if(indexed) {
		for(int i = 0; i < node.nvert; i++) {
			int l = last[i];
			Color4b &c = colors[i];
			for(int k = 0; k < 4; k++) {
				c[k] = decodeDiff(diffs[k][count], bitstream);
				if(l >= 0)
					c[k] += colors[l][k];
			}
			count++;
		}

	} else {
		Color4b on(0, 0, 0, 0);
		for(int i = 0; i < node.nvert; i++) {
			Color4b &c = colors[i];

			for(int k = 0; k < 4; k++)
				on[k] += decodeDiff(diffs[k][count], bitstream);
			count++;
			c = on;
		}
	}

	int steps[4];
	for(int k = 0; k < 4; k++)
		steps[k] = (1<<(8 - color_q[k]));

	for(int i = 0; i < node.nvert; i++) {
		Color4b &c = colors[i];
		for(int k = 0; k < 4; k++)
			c[k] *= steps[k];
		c = toRGB(c);

	}
}

void NxzDecoder::dequantize() {
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

template <class F>
void integrateNormals(int nvert, int nface, F *index, Point3f *coords, Point3f *normals) {
	for(int i = 0; i < nface; i++) {
		Point3f &p0 = coords[index[0]];
		Point3f &p1 = coords[index[1]];
		Point3f &p2 = coords[index[2]];
		Point3f n = (p1 - p0) ^ (p2 - p0);
		normal[index[0]] += n;
		normal[index[1]] += n;
		normal[index[2]] += n;
		index += 3;
	}
}

void NxzDecoder::computeNormals(Point3f *normals3f) {

	Point3f *coords = (Point3f *)coord.buffer;

	if(short_index)
		integrateNormals(nvert, nface, (uint16_t *)face, coords, normals3f);
	else
		integrateNormals(nvert, nface, (uint32_t *)face, coords, normals3f);

	for(unsigned int i = 0; i < normal.size(); i++) {
		Point3f &n = normal3f[i];
		n /= n.norm();
	}
}

void NxzDecoder::computeNormals(Point3s *normal3s) {
	vector<Point3f> normals3f(nvert, Point3f(0, 0, 0));
	computeNormals(&*normals3f.begin());

	for(unsigned int i = 0; i < nvert; i++) {
		Point3f &nf = normals3f[i];
		Point3s &ns = normals3s[i];
		for(int k = 0; k < 3; k++)
			normals3s[i][k] = (int16_t)(n[k]*32767.0f);
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
	for(uint32_t &end: groups) {
		decodeFaces(bitstream, start, end);
		start = end;
	}

	dequantize();
}

void NxzDecoder::decodeFaces(BitStream &bitstream, int nface, int count) {

	uint16_t *faces16 = ((uint16_t *)face.buffer);
	uint32_t *faces32 = ((uint32_t *)face.buffer);

	bool hasTextures = (flags & UV) != 0;

	Point3i *coords = (Point3i *)coord.buffer;
	Point2i *texcoords = (Point2i *)uv.buffer;

	vector<uchar> clers;
	vector<uchar> diffs;
	vector<uchar> tdiffs;
	BitStream bitstream;


	Tunstall tunstall1;
	tunstall1.decompress(stream, diffs);
	if(hasTextures) {
		Tunstall tunstall2;
		tunstall2.decompress(stream, tdiffs);
	}



	unsigned int current = 0;          //keep track of connected component start

	vector<int> delayed;
	deque<int> faceorder;
	vector<DEdge2> front; // for fast mode we need to pop_front

	unsigned int totfaces = nface;

	int diff_count = 0;
	int tdiff_count = 0;
	int cler_count = 0;


	while(totfaces > 0) {

		if(!faceorder.size() && !delayed.size()) {

			if(current == nface) break; //no more faces to encode exiting

			vcg::Point3i estimated(0, 0, 0);
			vcg::Point2i texestimated(0, 0);

			int last_index = -1;
			int index[3];
			for(int k = 0; k < 3; k++) {
				last.push_back(last_index);
				int diff = diffs[diff_count++];
				int tdiff = diff && hasTextures? tdiffs[tdiff_count++] : 0;
				int v = decodeVertex(estimated, texestimated, bitstream, diff, tdiff);
				index[k] = v;
				if(short_index)
					faces16[count++] = v;
				else
					faces32[count++] = v;
				estimated = coords[v];
				if(sig.vertex.hasTextures())
					texestimated = texcoords[v];
				last_index = v;
			}
			int current_edge = front.size();
			for(int k = 0; k < 3; k++) {
				faceorder.push_back(front.size());
				front.push_back(DEdge2(index[next_(k)], index[prev_(k)], index[k], current_edge + prev_(k), current_edge + next_(k)));
			}

			current++;
			totfaces--;
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

		assert(cler_count < clers.size());
		int c = clers[cler_count++];
		if(c == BOUNDARY) continue;


		int first_edge = front.size();
		int opposite = -1;

		if(c == VERTEX) {
			//we need to estimate opposite direction using v0 + v1 -
			int v2 = e.v2;

			Point3i predicted = coords[v0] + coords[v1] - coords[v2];
			Point2i texpredicted(0, 0);
			if(sig.vertex.hasTextures())
				texpredicted = texcoords[v0] + texcoords[v1] - texcoords[v2];

			int diff = diffs[diff_count++];
			int tdiff = diff && hasTextures? tdiffs[tdiff_count++] : 0;
			opposite = decodeVertex(predicted, texpredicted, bitstream, diff, tdiff);
			if(diff != 0)
				last.push_back(v1); //not v0!, in encode we use face, here edge which is inverted



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
			faces16[count++] = v1;
			faces16[count++] = v0;
			faces16[count++] = opposite;
		} else {
			faces32[count++] = v1;
			faces32[count++] = v0;
			faces32[count++] = opposite;
		}
		totfaces--;
	}
}

int NxzDecoder::decodeVertex(const Point3i &predicted, const Point2i &texpredicted, BitStream &bitstream, int diff, int tdiff) {

	static int count = 0;
	count++;
	if(diff == 0) {
		uint64_t bits = 0;
		bitstream.read(16, bits);
		assert(bits < node.nvert);
		return int(bits);
	}

	int v = vertex_count++;


	Point3i &p = *(((Point3i *)coord.buffer)[v]);
	decodeDiff(diff, bitstream, p);
	p += predicted;


	if(flags & UV) {
		//TODO these needds to be precomputed.
		Point2i &u = *(((Point2i *)&uv.buffer)[v]);
		decodeDiff(tdiff, bitstream, u);
		u += texpredicted;
	}
	return v;
}




int NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream) {
	if(diff == 0)
		return 0;

	//TODO could we use the max system to make it quicker?
	uint64_t c = 1<<(diff);
	bitstream.read(diff, c);
	int val = (int)c;

	if(val & 0x1)
		val >>= 1;
	else
		val = -(val >> 1);
	return val;
}


void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point3i &p) {
	//assert(diff < 22);

	//int mask = (1<<diff)-1;
	int max = (1<<(diff-1));
	p[0] = bitstream.readUint32(diffs) - max;
	p[1] = bitstream.readUint32(diffs) - max;
	p[2] = bitstream.readUint32(diffs) - max;
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
	p[0] = bitstream.readUint32(diffs) - max;
	p[1] = bitstream.readUint32(diffs) - max;
	//TODO is this really faster than just do multiple readings?

	/*	bitstream.read(diff*3, bits);
	for(int k = 2; k >= 0; k--) {
		p[k] = (bits & mask) - max;
		bits >>= diff;
	} */
}



