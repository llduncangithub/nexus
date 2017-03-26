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
	stream.entropy = (Stream::Entropy)stream.read<uchar>();

	int nattr = stream.read<int>();

	for(int i = 0; i < nattr; i++) {

		std::string name =  stream.readString();
		float q = stream.read<float>();
		uint32_t components = stream.read<uchar>();
		uint32_t strategy = stream.read<uchar>();

		Attribute23 *attr = nullptr;
		if(name == "position")
			attr = new GenericAttr<int>(components);
		else if(name == "normal")
			attr = new NormalAttr();
		else if(name == "color")
			attr = new ColorAttr();
		else if(name == " uv")
			attr = new GenericAttr<int>(components);

		attr->q = q;
		attr->strategy = strategy;
		data[name] = attr;
	}
}

bool NxzDecoder::setAttribute(const char *name, char *buffer, Attribute23::Format format) {
	if(data.find(name) == data.end()) return false;
	Attribute23 *attr = data[name];
	attr->format = format;
	attr->buffer = buffer;
	return true;
}

bool NxzDecoder::setAttribute(const char *name, char *buffer, Attribute23 *attr) {
	if(data.find(name) == data.end()) return false;
	Attribute23 *found = data[name];
	attr->q = found->q;
	attr->strategy = found->strategy;
	attr->N = found->N;
	attr->buffer = buffer;
	delete data[name];
	data[name] = attr;
	return true;
}


void NxzDecoder::decode() {
	nvert = stream.read<uint32_t>();
	nface = stream.read<uint32_t>();

	if(nface > 0)
		decodeMesh();
	else
		decodePointCloud();
}

void NxzDecoder::decodePointCloud() {

	std::vector<nx::Face> dummy;

	for(auto it: data)
		it.second->decode(nvert, stream);
	for(auto it: data)
		it.second->deltaDecode(nvert, dummy);
/*	for(auto it: data)
		it.second->postDelta(nvert, dummy); */
	for(auto it: data)
		it.second->dequantize(nvert);

}
/*
	Zpoint encoding gain 1 bit (because we know it's sorted: diffs are positive!, but it's 3/2 slower and limited to 22 bits precision.

	Point3f *coords = (Point3f *)coord.buffer;

	ZPoint p;
	bitstream.read(63, p.bits);
	coords[0] = p.toPoint(coord.q);
	for(size_t i = 1; i < nvert; i++) {
		uchar d = diffs[i-1];
		p.setBit(d, 1);
		bitstream.read(d, p.bits);
		coords[i] = coords[0] + p.toPoint(coord.q);
	}
	*/

void NxzDecoder::decodeMesh() {
	groups.resize(stream.read<int>());
	for(uint32_t &g: groups)
		g = stream.read<int>();

	stream.decompress(clers);
	BitStream bitstream;
	stream.read(bitstream);

	for(auto it: data)
		it.second->decode(nvert, stream);

	prediction.resize(nvert);

	uint32_t start = 0;
	uint32_t cler = 0; //keeps track of current cler
	for(uint32_t &end: groups) {
		decodeFaces(start*3, end*3, cler, bitstream);
		start = end;
	}
#ifdef PRESERVED_UNREFERENCED
	//decode unreferenced vertices
	while(vertex_count < nvert) {
		int last = vertex_count-1;
		prediction[vertex_count++] = Face(last, last, last);
	}
#endif

	for(auto it: data)
		it.second->deltaDecode(nvert, prediction);

	//TODO cleanup this mess.
	std::vector<uint32_t> index(nface*3);
	uint16_t *faces16 = ((uint16_t *)face.buffer);
	uint32_t *faces32 = ((uint32_t *)face.buffer);
	if(short_index)
		for(uint32_t i = 0; i < nface*3; i++)
			index[i] = faces16[i];
	else
		for(uint32_t i = 0; i < nface*3; i++)
			index[i] = faces32[i];

	for(auto it: data)
		it.second->postDelta(nvert, data, index);

	for(auto it: data)
		it.second->dequantize(nvert);
}

static int ilog2(uint64_t p) {
	int k = 0;
	while ( p>>=1 ) { ++k; }
	return k;
}

void NxzDecoder::decodeFaces(uint32_t start, uint32_t end, uint32_t &cler, BitStream &bitstream) {

	uint16_t *faces16 = ((uint16_t *)face.buffer);
	uint32_t *faces32 = ((uint32_t *)face.buffer);

	//edges of the mesh to be processed
	vector<DEdge2> front;
	front.reserve((end - start)*3);

	//faceorder is used to minimize split occourrence positioning in front and in back the new edges to be processed.
	vector<int>faceorder;
	faceorder.reserve((end - start)/2);
	uint32_t order = 0;

	//delayed again minimize split by further delay problematic splits
	vector<int> delayed;

	//TODO test if recording number of bits needed for splits improves anything. (very small but cost is zero.
	int splitbits = ilog2(nvert) + 1;

	int new_edge = -1; //last edge added which sohuld be the first to be processed, no need to store it in faceorder.

	while(start < end) {
		if(new_edge == -1 && order >= faceorder.size() && !delayed.size()) {

			int last_index = vertex_count-1;
			int index[3];

			int split =  0; //bitmask for vertex already decoded/
			if(clers[cler] == SPLIT) { //lookahead
				cler++;
				split = bitstream.readUint(3);
			}

			for(int k = 0; k < 3; k++) {
				int v; //TODO just use last_index.
				if(split & (1<<k))
					v = bitstream.readUint(splitbits);
				else {
					assert(vertex_count < (int)prediction.size());
					prediction[vertex_count] = Face(last_index, last_index, last_index);
					last_index = v = vertex_count++;
				}
				index[k] = v;
				if(short_index)
					faces16[start++] = v;
				else
					faces32[start++] = v;
			}
			int current_edge = front.size();
			faceorder.push_back(front.size());
			front.emplace_back(index[1], index[2], index[0], current_edge + 2, current_edge + 1);
			faceorder.push_back(front.size());
			front.emplace_back(index[2], index[0], index[1], current_edge + 0, current_edge + 2);
			faceorder.push_back(front.size());
			front.emplace_back(index[0], index[1], index[2], current_edge + 1, current_edge + 0);
			continue;
		}

		int f;
		if(new_edge != -1) {
			f = new_edge;
			new_edge = -1;

		} else if(order < faceorder.size()) {
			f = faceorder[order++];
		} else if (delayed.size()){
			f = delayed.back();
			delayed.pop_back(); //or popfront?

		} else {
			throw "Decoding topology failed";
		}

		DEdge2 &e = front[f];
		if(e.deleted) continue;
		//e.deleted = true; //each edge is processed once at most.

		int c = clers[cler++];
		if(c == BOUNDARY) continue;

		int v0 = e.v0;
		int v1 = e.v1;

		DEdge2 &previous_edge = front[e.prev];
		DEdge2 &next_edge = front[e.next];

		new_edge = front.size(); //index of the next edge to be added.
		int opposite = -1;

		if(c == VERTEX) {
			if(clers[cler] == SPLIT) { //lookahead
				cler++;
				opposite = bitstream.readUint(splitbits);
			} else {
				//Edge is inverted respect to encoding hence v1-v0 inverted.
				prediction[vertex_count] = Face(v1, v0, e.v2);
				opposite = vertex_count++;
			}

			previous_edge.next = new_edge;
			next_edge.prev = new_edge + 1;

			front.emplace_back(v0, opposite, v1, e.prev, new_edge + 1);
			faceorder.push_back(front.size());
			front.emplace_back(opposite, v1, v0, new_edge, e.next);

		} else if(c == LEFT) {
			previous_edge.deleted = true;
			front[previous_edge.prev].next = new_edge;
			front[e.next].prev = new_edge;
			opposite = previous_edge.v0;

			front.emplace_back(opposite, v1, v0, previous_edge.prev, e.next);

		} else if(c == RIGHT) {
			next_edge.deleted = true;
			front[next_edge.next].prev = new_edge;
			front[e.prev].next = new_edge;
			opposite = next_edge.v1;

			front.emplace_back(v0, opposite, v1, e.prev, next_edge.next);

		} else if(c == DELAY) {
			e.deleted = false;
			delayed.push_back(f);
			new_edge = -1;
			continue;

		} else if(c == END) {
			previous_edge.deleted = true;
			next_edge.deleted = true;
			front[previous_edge.prev].next = next_edge.next;
			front[next_edge.next].prev = previous_edge.prev;
			opposite = previous_edge.v0;
			new_edge = -1;

		} else {
			assert(0);
		}
		assert(v0 != v1);
		assert(v1 != opposite);
		assert(v0 != opposite);

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
