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
	stream.entropy = (Stream::Entropy)stream.read<uchar>();

	coord.q = stream.read<float>();
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
	data.resize(stream.read<int>());
	for(auto &da: data) {
		da.q = stream.read<float>();
		da.o = stream.read<int>();
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
	stream.decompress(diffs);

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
	stream.decompress(diffs);
	BitStream bitstream;
	stream.read(bitstream);

	Point3i *coords = (Point3i *)coord.buffer;
	for(uint32_t i = 0; i < nvert; i++)
		decodeDiff(diffs[i], bitstream, coords[i]);
}

void NxzDecoder::decodeNormals() {
	std::vector<uchar> diffs;
	stream.decompress(diffs);

	BitStream bitstream;
	stream.read(bitstream);
	if(!norm.buffer)
		return;

	//if normals_prediction == BORDER, diffs.size() will be small.
	//cheating to know the size in computeNormals.
	norm.size = diffs.size();

	Point3s *normals = (Point3s *)norm.buffer;
	Point3i *normali = (Point3i *)norm.buffer;
	for(uint32_t i = 0; i < diffs.size(); i++) {
		int pos = i;
		//store the difference at the back of the array, otherwise would be ovewritten
		if(normals_prediction == BORDER)
			pos += nvert - diffs.size();
		if(short_normals) {
			Point2s &n = *(Point2s *)&(normals[pos]);
			decodeDiff(diffs[i], bitstream, n);
		} else {
			Point2i &n = *(Point2i *)&(normali[pos]);
			decodeDiff(diffs[i], bitstream, n);
		}
	}
}

void NxzDecoder::decodeColors() {

	Color4b *colors = (Color4b *)color[0].buffer;

	vector<uchar> diffs;
	BitStream bitstream;

	for(int k = 0; k < 4; k++) {
		stream.decompress(diffs);
		stream.read(bitstream);
		if(colors)
			for(uint32_t i = 0; i < nvert; i++)
				colors[i][k] = decodeDiff(diffs[i], bitstream);
	}
}

void NxzDecoder::decodeUvs() {
	std::vector<uchar> diffs;
	stream.decompress(diffs);

	BitStream bitstream;
	stream.read(bitstream);

	if(!uv.buffer)
		return;

	Point2i *uvs = (Point2i *)uv.buffer;
	for(uint32_t i = 0; i < nvert; i++)
		decodeDiff(diffs[i], bitstream, uvs[i]);
}

void NxzDecoder::decodeDatas() {
	for(auto &da: data) {
		std::vector<uchar> diffs;
		stream.decompress(diffs);

		BitStream bitstream;
		stream.read(bitstream);

		if(!da.buffer) continue;
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

	if(flags & NORMAL && norm.buffer) {
		Point3s *normals = (Point3s *)norm.buffer;
		Point3f *normalf = (Point3f *)norm.buffer;

		if(normals_prediction == DIFF) {
			for(uint32_t i = 0; i < nvert; i++) {
				if(short_normals) {
					Point3s &n = normals[i];
					n = Normal::decode(Point2s(n[0], n[1]), (int)norm.q);
				} else {
					Point3f &n = normalf[i];
					Point2i ni = *(Point2i *)&n;
					n = Normal::decode(ni, (int)norm.q);
				}
			}
		} else {
			//decode now!
			vector<Point3f> estimated(nvert, Point3f(0, 0, 0));
			estimateNormals(&*estimated.begin());

			if(normals_prediction == BORDER) {
				if(short_index)
					markBoundary<uint16_t>(); //mark boundary points on original vertices.
				else
					markBoundary<uint32_t>();
			}
			if(short_normals)
				computeNormals(normals, &*estimated.begin());
			else
				computeNormals(normalf, &*estimated.begin());
		}
	}

	if(flags & COLOR && color[0].buffer) {
		Color4b *colors = (Color4b *)color[0].buffer;
		for(uint32_t i = 0; i < nvert; i++) {
			Color4b &c = colors[i];

			c = c.toRGB();
			for(int k = 0; k < 4; k++)
				c[k] *= (int)color[k].q;
		}
	}
}


template <class F>
void integrateNormals(int nface, F *index, Point3f *coords, Point3f *normals) {
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

void NxzDecoder::estimateNormals(Point3f *normals3f) {

	Point3f *coords = (Point3f *)coord.buffer;

	if(short_index)
		integrateNormals<uint16_t>(nface, (uint16_t *)face.buffer, coords, normals3f);
	else
		integrateNormals<uint32_t>(nface, (uint32_t *)face.buffer, coords, normals3f);

	for(unsigned int i = 0; i < nvert; i++) {
		Point3f &n = normals3f[i];
		n /= n.norm();
	}
}

void NxzDecoder::computeNormals(Point3s *normals, Point3f *estimated) {
	//for BORDERS we have a few diffs stored in the normals,
	//at the back of the array
	int count = norm.size; //here for the border.
	for(unsigned int i = 0; i < nvert; i++) {
		Point3f &e = estimated[i];
		Point3s &d = normals[count++];
		Point3s &n = normals[i];

		if(normals_prediction == ESTIMATED || boundary[i]) {
			Point2i qn = Normal::encode(e, (int)norm.q);
			n = Normal::decode(Point2s(qn[0] + d[0], qn[0] + d[1]), (int)norm.q);
		} else {//no correction
			for(int k = 0; k < 3; k++)
				n[k] = (int16_t)(e[k]*32767.0f);
		}
	}
}

void NxzDecoder::computeNormals(Point3f *normals, Point3f *estimated) {
	//for BORDERS we have a few diffs stored in the normals,
	//at the back of the array
	int count = norm.size; //here for the border.
	Point3i *diffs = (Point3i *)normals;
	for(unsigned int i = 0; i < nvert; i++) {
		Point3f &e = estimated[i];
		Point3i &d = diffs[count++];
		Point3f &n = normals[i];
		if(normals_prediction == ESTIMATED || boundary[i]) {
			Point2i qn = Normal::encode(e, (int)norm.q);
			n = Normal::decode(Point2i(qn[0] + d[0], qn[0] + d[1]), (int)norm.q);
		} else //no correction
			n = e;
	}
}

void NxzDecoder::decodeMesh() {
	nface = stream.read<int>();
	groups.resize(stream.read<int>());
	for(uint32_t &g: groups)
		g = stream.read<int>();

	stream.decompress(clers);
	BitStream bitstream;
	stream.read(bitstream);

	decodeCoords();

	//HERE I NEED THE NORMALS FOR THE BORDER! SIGH
	if(flags & NORMAL)
		decodeNormals();

	if(flags & COLOR)
		decodeColors();

	if(flags & UV)
		decodeUvs();

	decodeDatas();

	prediction.reserve(nvert*3);

	uint32_t start = 0;
	uint32_t cler = 0; //keeps track of current cler
	for(uint32_t &end: groups) {
		decodeFaces(start*3, end*3, cler, bitstream);
		start = end;
	}


	Point3i *coords = (Point3i *)coord.buffer;
	for(uint32_t i = 1, j = 3; i < nvert; i++, j+=3) {
		coords[i] += coords[prediction[j]] + coords[prediction[j+1]] - coords[prediction[j+2]];
	}
	//coords[v] += coords[v0] + coords[v1] - coords[v2];

	if(flags & UV) {
		Point2i *uvs = (Point2i *)uv.buffer;
		for(uint32_t i = 1, j = 3; i < nvert; i++, j+=3) {
			uvs[i] += uvs[prediction[j]] + uvs[prediction[j+1]] - uvs[prediction[j+2]];
		}
	}

	if(flags & NORMAL) {
		for(uint32_t i = 1, j = 3; i < nvert; i++, j += 3) {
			if(normals_prediction == DIFF) {
				if(short_normals) {
					Point3s &n = ((Point3s *)norm.buffer)[i];
					Point3s &ref = ((Point3s *)norm.buffer)[prediction[j]];
					//TODO slightly faster to avoid +=.
					n[0] += ref[0];
					n[1] += ref[1];
					if(n[0] < -norm.q)      n[0] += 2*norm.q;
					else if(n[0] > +norm.q) n[0] -= 2*norm.q;
					if(n[1] < -norm.q)      n[1] += 2*norm.q;
					else if(n[1] > +norm.q) n[1] -= 2*norm.q;

				} else {
					Point3i &n = ((Point3i *)norm.buffer)[i];
					Point3i &ref = ((Point3i *)norm.buffer)[prediction[j]];
					n[0] += ref[0];
					n[1] += ref[1];
					if(n[0] < -norm.q)      n[0] += 2*norm.q;
					else if(n[0] > +norm.q) n[0] -= 2*norm.q;
					if(n[1] < -norm.q)      n[1] += 2*norm.q;
					else if(n[1] > +norm.q) n[1] -= 2*norm.q;
				}
			}
		}
	}
	if(flags & COLOR && color[0].buffer) {
		for(uint32_t i = 1, j = 3; i < nvert; i++) {
			Color4b &c = ((Color4b *)color[0].buffer)[i];
			c += ((Color4b *)color[0].buffer)[prediction[j]];
		}
	}
/*
	for(auto &da: data) {
		int &d = ((int *)da.buffer)[v];
		d += ((int *)da.buffer)[v0];
	} */



	dequantize();
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
		if(new_edge == -1 && !faceorder.size() && !delayed.size()) {

			int last_index = vertex_count-1;
			int index[3];

			int split =  0; //bitmask for vertex already decoded/
			if(clers[cler] == SPLIT) { //lookahead
				cler++;
				split = bitstream.readUint(3);
			}

			for(int k = 0; k < 3; k++) {
				int v;
				if(split & (1<<k))
					v = bitstream.readUint(splitbits);
				else {
					prediction.push_back(last_index);
					prediction.push_back(last_index);
					prediction.push_back(last_index);
					v = vertex_count++;
//					v = decodeVertex(last_index, last_index, last_index);
				}
				index[k] = v;
				if(short_index)
					faces16[start++] = v;
				else
					faces32[start++] = v;

				last_index = v;
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
		} else {
			f = delayed.back();
			delayed.pop_back(); //or popfront?
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
				prediction.push_back(v1);
				prediction.push_back(v0);
				prediction.push_back(e.v2);
				opposite = vertex_count++;
				//opposite = decodeVertex(v1, v0, e.v2);
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

int NxzDecoder::decodeVertex(int v0, int v1, int v2) {

	uint32_t v = vertex_count++;
	assert(v < nvert);

	if(v0 == -1) //first vertex of component. nothing to do.
		return v;

	Point3i *coords = (Point3i *)coord.buffer;
	coords[v] += coords[v0] + coords[v1] - coords[v2];

	if(flags & UV) {
		Point2i *uvs = (Point2i *)uv.buffer;
		uvs[v] += uvs[v0] + uvs[v1] - uvs[v2];
	}

	if(flags & NORMAL) {
		if(normals_prediction == DIFF) {
			if(short_normals) {
				Point3s &n = ((Point3s *)norm.buffer)[v];
				Point3s &ref = ((Point3s *)norm.buffer)[v0];
				//TODO slightly faster to avoid +=.
				n[0] += ref[0];
				n[1] += ref[1];
				if(n[0] < -norm.q)      n[0] += 2*norm.q;
				else if(n[0] > +norm.q) n[0] -= 2*norm.q;
				if(n[1] < -norm.q)      n[1] += 2*norm.q;
				else if(n[1] > +norm.q) n[1] -= 2*norm.q;

			} else {
				Point3i &n = ((Point3i *)norm.buffer)[v];
				Point3i &ref = ((Point3i *)norm.buffer)[v0];
				n[0] += ref[0];
				n[1] += ref[1];
				if(n[0] < -norm.q)      n[0] += 2*norm.q;
				else if(n[0] > +norm.q) n[0] -= 2*norm.q;
				if(n[1] < -norm.q)      n[1] += 2*norm.q;
				else if(n[1] > +norm.q) n[1] -= 2*norm.q;
			}
		}
	}
	if(flags & COLOR && color[0].buffer) {
		Color4b &c = ((Color4b *)color[0].buffer)[v];
		c += ((Color4b *)color[0].buffer)[v0];
	}

	for(auto &da: data) {
		int &d = ((int *)da.buffer)[v];
		d += ((int *)da.buffer)[v0];
	}

	return v;
}




int NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream) {
	if(diff == 0)
		return 0;

	int val = bitstream.readUint(diff);
	int middle = 1<<(diff-1);
	if(val < middle)
		val = -val -middle;
	return val;
}

void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point3i &p) {
	if(diff == 0) {
		p[0] = p[1] = p[2] = 0;
		return;
	}
	int max = 1<<(diff-1);
	p[0] = bitstream.readUint(diff) - max;
	p[1] = bitstream.readUint(diff) - max;
	p[2] = bitstream.readUint(diff) - max;
}

void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point3s &p) {
	//assert(diff < 22);

	if(diff == 0) {
		p[0] = p[1] = p[2] = 0;
		return;
	}

	int max = 1<<(diff-1);
	p[0] = bitstream.readUint(diff) - max;
	p[1] = bitstream.readUint(diff) - max;
	p[2] = bitstream.readUint(diff) - max;
}

void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point2i &p) {
	//assert(diff < 22);
	if(diff == 0) {
		p[0] = p[1] =0;
		return;
	}

	int max = 1<<(diff-1);
	p[0] = bitstream.readUint(diff) - max;
	p[1] = bitstream.readUint(diff) - max;
}

void NxzDecoder::decodeDiff(uchar diff, BitStream &bitstream, Point2s &p) {
	//assert(diff < 22);

	if(diff == 0) {
		p[0] = p[1] = 0;
		return;
	}

	int max = 1<<(diff-1);
	p[0] = bitstream.readUint(diff) - max;
	p[1] = bitstream.readUint(diff) - max;
}
