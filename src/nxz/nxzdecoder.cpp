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

void NxzDecoder::NxzDecoder(int len, uchar *input):
	vertex_count(0), short_normals(false), short_index(false) {

	stream.init(len, input);

	flags = stream.read<int>();
	nvert = stream.read<int>();
	entropy = stream.read<char>();

	coord.q = stream.read<char>();
	coord.o[0] = stream.read<int>();
	coord.o[1] = stream.read<int>();
	coord.o[2] = stream.read<int>();

	if(flags & NORMAL) {
		normals = stream.read<uchar>();
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

void NxzDecoder::setData(int pos, float *buffer) { data[i].buffer = buffer; }

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




void NxzDecoder::decodeCoordinates() {
	std::vector<uchar> diffs;
	Tunstall::decompress(stream, diffs);

	BitStream bitstream;
	stream.read(bitstream);

	Point3i *coords = (Point3i *)coord.buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		decodeDiff(diffs, bitstream, coords[i]);
	}
}

void NxzDecoder::dequantize() {
	{
		float step = pow(2.0f, (float)coord_q);
		vcg::Point3f *coords = data.coords();
		for(int i = 0; i < node.nvert; i++) {
			Point3i p = *(Point3i *)&coords[i];
			p[0] += min[0];
			p[1] += min[1];
			p[2] += min[2];
			coords[i] = Point3f(p[0]*step, p[1]*step, p[2]*step);
		}
	}
	if(sig.vertex.hasTextures()) {
		float step = pow(2.0f, (float)uv_q);
		vcg::Point2f *texcoords = data.texCoords(sig, node.nvert);
		for(int i = 0; i < node.nvert; i++) {
			Point2i p = *(Point2i *)&texcoords[i];
			p[0] += tmin[0];
			p[1] += tmin[1];
			texcoords[i] = Point2f(p[0]*step, p[1]*step);
		}
	}
}

void NxzDecoder::decodeNormals() {

	Point3s *norms = data.normals(sig, node.nvert);
	norm_q = stream.read<char>();

	std::vector<uchar> diffs;
	Tunstall tunstall;
	tunstall.decompress(stream, diffs);

	std::vector<uchar> signs;
	Tunstall tunstall1;
	tunstall1.decompress(stream, signs);

	BitStream bitstream;
	stream.read(bitstream);

	int side = 1<<(16 - norm_q);

	if(!sig.face.hasIndex()) {
		int count = 0;
		for(int k = 0; k < 2; k++) {

			int on = 0; //old normal;
			for(unsigned int i = 0; i < node.nvert; i++) {
				Point3s &c = norms[i];
				int n = c[k];
				//int16_t &n = c[k];

				int d = decodeDiff(diffs[count++], bitstream);
				n = on + d;
				on = n;
				n *= side;
				if(n >=32768) n = 32767;
				c[k] = n;
			}
		}

		for(unsigned int i = 0; i < node.nvert; i++) {
			Point3s &c = norms[i];
			double z = 32767.0 * 32767.0 - c[0]*c[0] - c[1]*c[1];
			if(z < 0) z = 0;
			z = sqrt(z);
			if(z > 32767.0)
				z = 32767.0;
			c[2] = (int16_t)z; //sqrt(32767.0 * 32767.0 - c[0]*c[0] - c[1]*c[1]);
			assert(c[2] >= 0);
			if(!signs[i]) {
				c[2] = -c[2];
			}
		}

	} else {

		markBoundary();
		Point3s *estimated_normals = data.normals(sig, node.nvert);
		computeNormals(estimated_normals);
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
	}
}

void NxzDecoder::markBoundary() {
	if(!sig.face.hasIndex()) {
		boundary.resize(node.nvert, true);
		return;
	}
	boundary.resize(node.nvert, false);

	uint16_t *faces = data.faces(sig, node.nvert);
	vector<int> count(node.nvert, 0);
	for(int i = 0; i < node.nface; i++) {
		uint16_t *f = faces + i*3;
		count[f[0]] += (int)f[1] - (int)f[2];
		count[f[1]] += (int)f[2] - (int)f[0];
		count[f[2]] += (int)f[0] - (int)f[1];
	}
	for(int i = 0; i < node.nvert; i++)
		if(count[i] != 0)
			boundary[i] = true;
}


void NxzDecoder::computeNormals(Point3i *estimated_normals) {

	Point3f *coords = data.coords();
	uint16_t *faces = data.faces(sig, node.nvert);
	vector<Point3f> normals(node.nvert, Point3f(0, 0, 0));

	for(int i = 0; i < node.nface; i++) {
		uint16_t *face = faces + i*3;
		Point3f &p0 = coords[face[0]];
		Point3f &p1 = coords[face[1]];
		Point3f &p2 = coords[face[2]];
		Point3f n = (( p1 - p0) ^ (p2 - p0));
		normals[face[0]] += n;
		normals[face[1]] += n;
		normals[face[2]] += n;
	}
	//normalize
	for(unsigned int i = 0; i < normals.size(); i++) {
		Point3f &n = normals[i];
		float norm = sqrt(float(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]));
		for(int k = 0; k < 3; k++)
			estimated_normals[i][k] = (int16_t)(n[k]*32767.0f/norm);
	}
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
	 last.reserve(node.nvert);

	 uint16_t *faces = data.faces(sig, node.nvert);

	 //decodeFaces(node.nface, faces);
	 int start = 0;
	 for(uint32_t p = node.first_patch; p < node.last_patch(); p++) {
		 Patch &patch = patches[p];
		 uint end = patch.triangle_offset;
		 decodeFaces(end - start, faces + start*3);
		 start = end;
	 }
	 dequantize();
 }

void NxzDecoder::decodeFaces(int nface, uint16_t *faces) {

	bool hasTextures = sig.vertex.hasTextures();

	Point3i *coords = (Point3i *)data.coords();
	Point2i *texcoords = (Point2i *)data.texCoords(sig, node.nvert);

	vector<uchar> clers;
	vector<uchar> diffs;
	vector<uchar> tdiffs;
	BitStream bitstream;

	Tunstall tunstall;
	tunstall.decompress(stream, clers);
	Tunstall tunstall1;
	tunstall1.decompress(stream, diffs);
	if(hasTextures) {
		Tunstall tunstall2;
		tunstall2.decompress(stream, tdiffs);
	}


	stream.read(bitstream);

	unsigned int current = 0;          //keep track of connected component start

	vector<int> delayed;
	deque<int> faceorder;
	vector<DEdge2> front; // for fast mode we need to pop_front

	int count = 0;            //keep track of decoded triangles
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
				faces[count++] = v;
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
		faces[count++] = v1;
		faces[count++] = v0;
		faces[count++] = opposite;
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

	uint64_t c = 1<<(diff);
	bitstream.read(diff, c);
	int val = (int)c;

	if(val & 0x1)
		val >>= 1;
	else
		val = -(val >> 1);
	return val;
}

void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point2i &p) {

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



