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


NxzEncoder::NxzEncoder(uint32_t _nvert, uint32_t _nface, Stream::Entropy entropy):
	nvert(_nvert), nface(_nface),
	header_size(0), current_vertex(0) {

	stream.entropy = entropy;
	index.resize(nface*3);
}

//TODO optional checking for invalid (nan, infinity) values.

/* TODO: euristic for point clouds should be: about 1/10 of nearest neighbor.
  for now just use volume and number of points
  could also quantize to find a side where the points halve and use 1/10 */

bool NxzEncoder::addPositions(float *buffer, float q, Point3f o) {
	std::vector<Point3f> coords(nvert);
	Point3f *input = (Point3f *)buffer;
	for(uint32_t i = 0; i < nvert; i++)
		coords[i] = input[i] - o;

	if(q == 0) { //assume point cloud and estimate quantization
		Point3f max(-FLT_MAX), min(FLT_MAX);
		for(uint32_t i = 0; i < nvert; i++) {
			min.setMin(coords[i]);
			max.setMax(coords[i]);
		}
		max -= min;
		q = 0.02*pow(max[0]*max[1]*max[2], 2.0/3.0)/nvert;
	}
	uint32_t strategy = Attribute23::CORRELATED;
	if(nface > 0)
		strategy |= Attribute23::PARALLEL;

	return addAttribute("position", (char *)buffer, Attribute23::FLOAT, 3, q, strategy );
}


/* if not q specified use 1/10th of average leght of edge  */
bool NxzEncoder::addPositions(float *buffer, uint32_t *_index, float q, Point3f o) {
	index.resize(nface*3);
	memcpy(&*index.begin(), _index,  nface*12);

	Point3f *coords = (Point3f *)buffer;
	if(q == 0) { //estimate quantization on edge length average.
		double average = 0; //for precision when using large number of faces
		for(uint32_t f = 0; f < nface*3; f += 3)
			average += (coords[index[f]] - coords[index[f+1]]).norm();
		q = (float)(average/nface)/20.0f;
	}
	return addPositions(buffer, q, o);
}

bool NxzEncoder::addPositions(float *buffer, uint16_t *_index, float q, Point3f o) {
	vector<uint32_t> tmp(nface*3);
	for(uint32_t i = 0; i < nvert*3; i++)
		tmp[i] = _index[i];
	return addPositions(buffer, &*tmp.begin(), q, o);
}

/*
	//TODOThis should be stored somewhere for testing purpouses.
	if(1) {
		Point3i cmin(2147483647);
		Point3i cmax(-2147483647);
		for(auto &q: coord.values) {
			cmin.setMin(q);
			cmax.setMax(q);
		}
		cmax -= cmin;

		int bits = 1+std::max(std::max(ilog2(cmax[0]), ilog2(cmax[1])), ilog2(cmax[2]));
		if(!(flags & INDEX) && bits >= 22)
			cerr << "Quantiziation above 22 bits for point clouds is not supported..." << endl;
		else
			cout << "Quantization coords in " << bits << " bits\n";
	}
} */

bool NxzEncoder::addNormals(float *buffer, int bits, NormalAttr::Prediction no) {

	NormalAttr *normal = new NormalAttr(bits);
	normal->format = Attribute23::FLOAT;
	normal->prediction = no;
	bool ok = addAttribute("normal", (char *)buffer, normal);
	if(!ok) delete normal;
	return ok;
}

bool NxzEncoder::addNormals(int16_t *buffer, int bits, NormalAttr::Prediction no) {
	assert(bits <= 16);
	vector<Point3f> tmp(nvert*3);
	for(uint32_t i = 0; i < nvert; i++)
		for(int k = 0; k < 3; k++)
			tmp[i][k] = buffer[i*3 + k]/32767.0f;
	return addNormals((float *)&*tmp.begin(), bits, no);
}

bool NxzEncoder::addColors(unsigned char *buffer, int lumabits, int chromabits, int alphabits) {
	ColorAttr *color = new ColorAttr();
	color->setQ(lumabits, chromabits, alphabits);
	color->format = Attribute23::UINT8;
	bool ok = addAttribute("color", (char *)buffer, color);
	if(!ok) delete color;
	return ok;
}

bool NxzEncoder::addUV(float *buffer, float q) {
	GenericAttr<int> *uv = new GenericAttr<int>(2);
	uv->q = q;
	uv->format = Attribute23::FLOAT;
	bool ok = addAttribute("uv", (char *)buffer, uv);
	if(!ok) delete uv;
	return ok;
}

void NxzEncoder::encode() {
	stream.reserve(nvert);

//	stream.write<int>(flags);
	stream.write<uchar>(stream.entropy);

	stream.write<int>(data.size());
	for(auto it: data) {
		stream.writeString(it.first.c_str());                //name
		stream.write<float>(it.second->q);
		stream.write<uchar>(it.second->N);
		stream.write<uchar>(it.second->strategy);
	}

	if(nface > 0)
		encodeMesh();
	else
		encodePointCloud();
}

//TODO: test pointclouds
void NxzEncoder::encodePointCloud() {
	//look for positions
	if(data.find("position") == data.end())
		throw "No position attribute found. Use DIFF normal strategy instead.";

	GenericAttr<int> *coord = dynamic_cast<GenericAttr<int> *>(data["position"]);
	if(!coord)
		throw "Position attr has been overloaded, Use DIFF normal strategy instead.";

	Point3i *coords = (Point3i *)&*coord->values.begin();

	std::vector<ZPoint> zpoints(nvert);

	Point3i min(0, 0, 0);
	for(uint32_t i = 0; i < nvert; i++) {
		Point3i &q = coords[i];
		min.setMin(q);
	}
	for(uint32_t i = 0; i < nvert; i++) {
		Point3i q = coords[i] - min;
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

        stream.write<uint32_t>(nvert);
	stream.write<uint32_t>(0); //nface


	prediction.resize(nvert);
	prediction[0] = Quad(zpoints[0].pos, -1, -1, -1);
	for(uint32_t i = 1; i < nvert; i++)
		prediction[i] = Quad(zpoints[i].pos, zpoints[i-1].pos, zpoints[i-1].pos, zpoints[i-1].pos);

	for(auto it: data)
		it.second->preDelta(nvert, data, index);

	for(auto it: data)
		it.second->deltaEncode(prediction);

	for(auto it: data)
		it.second->encode(nvert, stream);
}

/*	Zpoint encoding gain 1 bit (because we know it's sorted, but it's 3/2 slower and limited to 22 bits precision.
 *
 * bitstream.write(zpoints[0].bits, 63);
	for(uint32_t pos = 1; pos < nvert; pos++) {
		ZPoint &p = zpoints[pos-1]; //previous point
		ZPoint &q = zpoints[pos]; //current point
		uchar d = p.difference(q);
		diffs.push_back(d);
		//we can get away with d diff bits (should be d+1) because the first bit will always be 1 (sorted q> p)
		bitstream.write(q.bits, d); //rmember to add a 1, since
	} */

//compact in place faces in data, update patches information, compute topology and encode each patch.
void NxzEncoder::encodeMesh() {
	encoded.resize(nvert, -1);

	if(!groups.size()) groups.push_back(nface);
	//remove degenerate faces
	uint32_t start =  0;
	uint32_t count = 0;
	for(uint32_t &end: groups) {
		for(uint32_t i = start; i < end; i++) {
			uint32_t *f = &index[i*3];

			if(f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
				continue;

			if(count != i) {
				uint32_t *dest = &index[count*3];
				dest[0] = f[0];
				dest[1] = f[1];
				dest[2] = f[2];
			}
			count++;
		}
		start = end;
		end = count;
	}
	index.resize(count*3);
	nface = count;



	BitStream bitstream(nvert/4);
	prediction.resize(nvert);

	start =  0;
	for(uint32_t &end: groups) {
		encodeFaces(start, end, bitstream);
		start = end;
	}
#ifdef PRESERVED_UNREFERENCED
	//encoding unreferenced vertices
	for(uint32_t i = 0; i < nvert; i++) {
		if(encoded[i] != -1)
			continue;
		int last = cur`rent_vertex-1;
		prediction.emplace_back(Quad(i, last, last, last));
		current_vertex++;
	}
#endif
	cout << "Unreference vertices: " << nvert - current_vertex << " remaining: " << current_vertex << endl;
	nvert = current_vertex;
	prediction.resize(nvert);

	stream.write<int>(nvert);
	stream.write<int>(nface);

	stream.write<int>(groups.size());
	for(uint32_t &end: groups)
		stream.write<int>(end);

	header_size = stream.elapsed();

	for(auto it: data)
		it.second->preDelta(nvert, data, index);

	for(auto it: data)
		it.second->deltaEncode(prediction);


	stream.compress(clers.size(), &*clers.begin());
	stream.write(bitstream);
	index_size = stream.elapsed();

	for(auto it: data)
		it.second->encode(nvert, stream);

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
	uint32_t side; //opposited to side vertex of face (edge 2 opposite to vertex 2)
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

void NxzEncoder::encodeFaces(int start, int end, BitStream &bitstream) {

	vector<McFace> faces(end - start);
	for(int i = start; i < end; i++) {
		uint32_t * f = &index[i*3];
		faces[i - start] = McFace(f[0], f[1], f[2]);
		assert(f[0] != f[1] && f[1] != f[2] && f[2] != f[0]);
	}

	buildTopology(faces);

	unsigned int current = 0;          //keep track of connected component start

	vector<int> delayed;
	//TODO move to vector + order
	deque<int> faceorder;
	vector<CEdge> front;

	vector<bool> visited(faces.size(), false);
	unsigned int totfaces = faces.size();

	int splitbits = ilog2(nvert) + 1;

	int last_index = -1;
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

			//int last_index = current_vertex -1; //ecccalla' la cazzata. mi servirebbe l'ultimo NON sorted.
			int split = 0;
			for(int k = 0; k < 3; k++) {
				int vindex = face.f[k];
				if(encoded[vindex] != -1)
					split |= 1<<k;
			}

			if(split) {
				clers.push_back(SPLIT);
				bitstream.writeUint(split, 3);
			}

			for(int k = 0; k < 3; k++) {
				uint32_t vindex = face.f[k];
				int &enc = encoded[vindex];

				if(enc != -1) {
					bitstream.writeUint(enc, splitbits);
				} else {
					//quad uses presorting indexing. (diff in attribute are sorted, values are not).
					prediction[current_vertex] = Quad(vindex, last_index, last_index, last_index);
					enc = current_vertex++;
					last_index = vindex;
				}
			}

			faceorder.push_back(front.size());
			front.emplace_back(current, 0, current_edge + 2, current_edge + 1);
			faceorder.push_back(front.size());
			front.emplace_back(current, 1, current_edge + 0, current_edge + 2);
			faceorder.push_back(front.size());
			front.emplace_back(current, 2, current_edge + 1, current_edge + 0);


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
			front.emplace_back(opposite_face, k1, previous_edge.prev, enext);

		} else if(close_right) {
			clers.push_back(RIGHT);
			next_edge.deleted = true;
			front[next_edge.next].prev = first_edge;
			front[eprev].next = first_edge;
			faceorder.push_front(front.size());
			front.emplace_back(opposite_face, k0, eprev, next_edge.next);

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
			if(encoded[opposite] != -1) {
				clers.push_back(SPLIT);
				bitstream.writeUint(encoded[opposite], splitbits);

			} else {
				//vertex needed for parallelogram prediction
				int v2 = faces[e.face].f[e.side];
				prediction[current_vertex] = Quad(opposite, v0, v1, v2);
				encoded[opposite] = current_vertex++;
				last_index = opposite;
			}

			previous_edge.next = first_edge;
			next_edge.prev = first_edge + 1;
			faceorder.push_front(front.size());
			front.emplace_back(opposite_face, k0, eprev, first_edge+1);
			faceorder.push_back(front.size());
			front.emplace_back(opposite_face, k1, first_edge, enext);
		}

		counting++;
		assert(!visited[opposite_face]);
		visited[opposite_face] = true;
		totfaces--;
	}
}
