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
/*
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
		q = 0.01*pow(max[0]*max[1]*max[2], 2.0/3.0)/nvert;
		if(o == Point3f(FLT_MAX))
			o = min;
	}
	if(o == Point3f(FLT_MAX)) {
		for(uint32_t i = 0; i < nvert; i++)
			o.setMin(coords[i]);
	}

	coord.q = q;

	for(uint32_t i = 0; i < nvert; i++) {
		Point3i &q = coord.values[i];
		Point3f &p = coords[i];
		for(int k = 0; k < 3; k++)
			q[k] = floor(p[k]/coord.q);
	}
	setCoordBits();
} */

bool NxzEncoder::addPositions(Point3f *buffer, float q, Point3f o) {
	std::vector<Point3f> coords(nvert);
	for(uint32_t i = 0; i < nvert; i++)
		coords[i] = buffer[i] - o;

	if(q == 0) { //assume point cloud and estimate quantization
		Point3f max(-FLT_MAX), min(FLT_MAX);
		for(uint32_t i = 0; i < nvert; i++) {
			min.setMin(coords[i]);
			max.setMax(coords[i]);
		}
		max -= min;
		q = 0.02*pow(max[0]*max[1]*max[2], 2.0/3.0)/nvert;
	}

	return addAttribute("position", (char *)buffer, Attribute23::FLOAT, 3, q, Attribute23::CORRELATED | Attribute23::PARALLEL);
}


/* if not q specified use 1/10th of average leght of edge  */
bool NxzEncoder::addPositions(Point3f *buffer, uint32_t *_index, float q, Point3f o) {
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

bool NxzEncoder::addPositions(Point3f *buffer, uint16_t *_index, float q, Point3f o) {
	vector<uint32_t> tmp(nface*3);
	for(uint32_t i = 0; i < nvert*3; i++)
		tmp[i] = _index[i];
	return addPositions(buffer, &*tmp.begin(), q, o);
}

/*void NxzEncoder::setCoordBits() {
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


/* NORMALS */

bool NxzEncoder::addNormals(Point3f *buffer, int bits, NormalAttr::Prediction no) {

	NormalAttr *normal = new NormalAttr(bits);
	normal->format = Attribute23::FLOAT;
	bool ok = addAttribute("normal", (char *)buffer, normal);
	if(!ok) delete normal;
	return ok;

/*	normals_prediction = no;
	flags |= NORMAL;
	norm.values.resize(nvert);
	norm.diffs.resize(nvert);
	norm.q = pow(2, bits-1);

	Point3f *normals = (Point3f *)buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		Point2i dt = norm.values[i] = Normal::encode(normals[i], (int)norm.q);
		assert(abs(dt[0]) < (1<<25));
		assert(abs(dt[1]) < (1<<25));
	} */
}

bool NxzEncoder::addNormals(Point3s *buffer, int bits, NormalAttr::Prediction no) {
	assert(bits <= 16);
	vector<Point3f> tmp(nvert*3);
	for(uint32_t i = 0; i < nvert; i++)
		for(int k = 0; k < 3; k++)
			tmp[i][k] = buffer[i][k]/32767.0f;
	return addNormals(&*tmp.begin(), bits, no);
}

/* COLORS */

bool NxzEncoder::addColors(unsigned char *buffer, int lumabits, int chromabits, int alphabits) {
	ColorAttr *color = new ColorAttr();
	color->setQ(lumabits, chromabits, alphabits);
	color->format = Attribute23::UINT8;
	bool ok = addAttribute("color", (char *)buffer, color);
	if(!ok) delete color;
	return ok;
/*	flags |= COLOR;
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
		Color4b c = colors[i];
		Color4b y;
		for(int k = 0; k < 4; k++)
			y[k] = c[k]/q[k];
		y = y.toYCC();
		for(int k = 0; k < 4; k++)
			color[k].values[i] = y[k];
	}*/
}

bool NxzEncoder::addUV(float *buffer, float q) {
	GenericAttr<int> *uv = new GenericAttr<int>(2);
	uv->q = q;
	uv->format = Attribute23::FLOAT;
	bool ok = addAttribute("uv", (char *)buffer, uv);
	if(!ok) delete uv;
	return ok;

/*	flags |= UV;
	if(q == 0) q = pow(2, -12);
	uv.q = q;
	uv.values.resize(nvert);
	uv.diffs.resize(nvert);
	Point2f *uvs = (Point2f *)buffer;
	for(uint32_t i = 0; i < nvert; i++) {
		uv.values[i][0] = (int)round(uvs[i][0]/uv.q);
		uv.values[i][1] = (int)round(uvs[i][1]/uv.q);
	}

	if(1) {
		Point2i cmin( 2147483647);
		Point2i cmax(-2147483647);
		for(auto &q: uv.values) {
			cmax.setMax(q);
			cmin.setMin(q);
		}
		cmax -= cmin;
		int bits = 1 + std::max(ilog2(cmax[0]), ilog2(cmax[1]));
		cout << "Texture cooridnates quantization in " << bits << " bits\n";
	} */
}
/*
void NxzEncoder::addData(Attribute23 *attr, char *buffer, float q) {
	flags |= DATA*(1<<data.size());
	attr->q = q;
	attr->quantize(nvert, buffer);
	data.push_back(attr);

//	setDataBits();
}*/
/*
void NxzEncoder::addData(int *buffer) {
	flags |= DATA*(1<<data.size());
	Data<int, int> *attr = new Data<int, int>(1.0f);
	attr->quantize(nvert, (char *)buffer);
	data.push_back(attr);
} */

//void NxzEncoder::setDataBits() {
/*	auto &d = data.back();
	int cmin =  2147483647;
	int cmax = -2147483647;

	for(auto &q: d.values) {
		if(q < cmin) cmin = q;
		if(q > cmax) cmax = q;
	}
	cmax -= cmin;
	cout << "Data quantized in: " << 1 + ilog2(cmax) << " bits." << endl; */
//}

void NxzEncoder::encode() {
	stream.reserve(nvert);

//	stream.write<int>(flags);
	stream.write<uchar>(stream.entropy);

	stream.write<int>(data.size());
	for(auto it: data) {
		stream.writeString(it.first.c_str());                //name
		stream.write<uchar>(it.second->q);
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

	stream.write<uint32_t>(0); //nface
	stream.write<uint32_t>(nvert);

	prediction.resize(nvert);
	prediction[0] = Quad(zpoints[0].pos, -1, -1, -1);
	for(uint32_t i = 1; i < nvert; i++)
		prediction[i] = Quad(zpoints[i].pos, zpoints[i-1].pos, zpoints[i-1].pos, zpoints[i-1].pos);

/*
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
	}*/

	for(auto it: data) {
		it.second->deltaEncode(prediction);
		it.second->encode(nvert, stream);
	}
}

/*void NxzEncoder::encodeZPoints(std::vector<ZPoint> &zpoints) {
	vector<uchar> diffs;
	BitStream bitstream(nvert/2);

	coord.diffs.resize(nvert);
	coord.diffs[0] = coord.values[zpoints[0].pos];
	for(uint32_t i = 1; i < nvert; i++)
		coord.diffs[i] = coord.values[zpoints[i].pos] - coord.values[zpoints[i-1].pos];

	for(Point3i &d: coord.diffs)
		encodeDiff(diffs, bitstream, d);

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

/*	stream.compress(diffs.size(), &*diffs.begin());
	stream.write(bitstream);
	coord.size = stream.elapsed();
} */
/*
void NxzEncoder::encodeCoords() {
	//TODO: diffs + bitstream is an logentropy compressor make it a class.
	vector<uchar> diffs;
	BitStream bitstream(nvert/2);

	for(Point3i &d: coord.diffs)
		encodeDiff(diffs, bitstream, d);

	stream.compress(diffs.size(), &*diffs.begin());
	stream.write(bitstream);
	coord.size = stream.elapsed();
}

void NxzEncoder::encodeNormals() {
	vector<uchar> diffs;
	BitStream bitstream(nvert/64);

	//TOO: VS using Point2i diff
	for(uint32_t i = 0; i < nvert; i++)
		if(normals_prediction != BORDER || boundary[i])
			encodeDiff(diffs, bitstream, norm.diffs[i]);

	stream.compress(diffs.size(), &*diffs.begin());
	stream.write(bitstream);

	norm.size = stream.elapsed();
}

void NxzEncoder::encodeColors() {

	//TODO: Slightly more efficient to use only one bitstream,

	for(int k = 0; k < 4; k++) {
		BitStream bitstream(nvert/2);
		vector<uchar> diffs;

		//use complement to 255 (hence char conversion)

		for(uchar &c: color[k].diffs)
			encodeDiff(diffs, bitstream, (char)c);

		stream.compress(diffs.size(), &*diffs.begin());
		stream.write(bitstream);
		color[k].size = stream.elapsed();
	}
}


void NxzEncoder::encodeUvs() {
	vector<uchar> diffs;
	BitStream bitstream(nvert/2);

	for(Point2i &u: uv.diffs)
		encodeDiff(diffs, bitstream, u);

	stream.compress(diffs.size(), &*diffs.begin());
	stream.write(bitstream);
	uv.size = stream.elapsed();
}
*/

//void NxzEncoder::encodeDatas() {
/*	for(auto &da: data) {
		BitStream bitstream(nvert/2);
		vector<uchar> diffs;

		for(int &v: da.diffs)
			encodeDiff(diffs, bitstream, v);

		stream.compress(diffs.size(), &*diffs.begin());
		stream.write(bitstream);
		da.size += stream.elapsed();
	} */
//}



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



	BitStream bitstream;
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
		int last = current_vertex-1;
		prediction.emplace_back(Quad(i, last, last, last));
		current_vertex++;
	}
#endif
	cout << "Unreference vertices: " << nvert - current_vertex << " remaining: " << current_vertex << endl;
	nvert = current_vertex;

	stream.write<int>(nvert);
	stream.write<int>(nface);
	stream.write<int>(groups.size());
	for(uint32_t &end: groups)
		stream.write<int>(end);

	header_size = stream.elapsed();


/*	coord.diffs[0] = coord.values[prediction[0].t];
	for(uint32_t i = 1; i < nvert; i++) {
		Quad &q = prediction[i];
		coord.diffs[i] = coord.values[q.t] - (coord.values[q.a] + coord.values[q.b] - coord.values[q.c]);
	}

	if((flags & NORMAL)) {


		if((flags & NORMAL) && normals_prediction != DIFF) {
			computeNormals(norm.diffs);
			if(normals_prediction == BORDER)
				markBoundary(); //mark boundary points on original vertices.
			for(uint32_t i = 0; i < nvert; i++) {
				norm.values[i] -= norm.diffs[i];

				Point2i dt = norm.values[i];
				assert(abs(dt[0]) < (1<<25));
				assert(abs(dt[1]) < (1<<25));
			}
		}

		norm.diffs[0] = norm.values[prediction[0].t];
		for(uint32_t i = 1; i < nvert; i++) {
			Quad &q = prediction[i];

			Point2i &dt = norm.diffs[i];
			dt = norm.values[q.t];

			if(normals_prediction == DIFF) {
				dt -= norm.values[q.a];
				if(dt[0] < -norm.q)      dt[0] += 2*norm.q;
				else if(dt[0] > +norm.q) dt[0] -= 2*norm.q;
				if(dt[1] < -norm.q)      dt[1] += 2*norm.q;
				else if(dt[1] > +norm.q) dt[1] -= 2*norm.q;
			}
		}
	}

	if(flags & COLOR) {
		for(int k = 0; k < 4; k++)
			color[k].diffs[0] = color[k].values[prediction[0].t];
		for(uint32_t i = 1; i < nvert; i++) {
			Quad &q = prediction[i];
			for(int k = 0; k < 4; k++)
				color[k].diffs[i] = color[k].values[q.t] - color[k].values[q.a];
		}
	}

	if(flags & UV) {
		uv.diffs[0] = uv.values[prediction[0].t];
		for(uint32_t i = 1; i < nvert; i++) {
			Quad &q = prediction[i];
			uv.diffs[i] = uv.values[q.t] - (uv.values[q.a] + uv.values[q.b] - uv.values[q.c]);
		}
	} */

	for(auto it: data)
		it.second->deltaEncode(prediction);





	stream.compress(clers.size(), &*clers.begin());
	stream.write(bitstream);
	index_size = stream.elapsed();

/*(	encodeCoords();

	if(flags & NORMAL)
		encodeNormals();

	if(flags & COLOR)
		encodeColors();

	if(flags & UV)
		encodeUvs(); */

	for(auto it: data)
		it.second->encode(nvert, stream);

}

//How to determine if a vertex is a boundary without topology:
//for each edge a vertex is in, add or subtract the id of the other vertex depending on order
//for internal vertices sum is zero.
//unless we have non manifold strange configurations zero wont happen.
/*
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
} */
/*

void NxzEncoder::computeNormals(vector<Point2i> &estimated) {
	vector<Point3f> tmp(nvert, Point3f(0, 0, 0));
	estimated.resize(nvert);
	for(uint32_t i = 0; i < nface; i++) {
		uint32_t *f = &index[i*3];
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
		float len = n.norm();
		if(len < 0.00001) {
			estimated[i] = Point2i(0, 0);
			continue;
		}
		n /= len;
		estimated[i] = Normal::encode(n, norm.q);
		assert(abs(estimated[i][0]) < (1<<25));
		assert(abs(estimated[i][1]) < (1<<25));
	}
}
*/
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
			//Point3i coord_estimated(0, 0, 0);
			//Point2i uv_estimated(0, 0);
			int last_index = current_vertex -1;
			int split = 0;
			for(int k = 0; k < 3; k++) {
				int vindex = face.f[k];
				if(encoded[vindex] != -1)
					split |= (1<<k);
			}

			if(split) {
				clers.push_back(SPLIT);
				bitstream.writeUint(split, 3);
			}

			for(int k = 0; k < 3; k++) {
				int vindex = face.f[k];

				if(split & (1<<k))
					bitstream.writeUint(encoded[vindex], splitbits);
				else {
					prediction[current_vertex] = Quad(vindex, last_index, last_index, last_index);
					encoded[vindex] = current_vertex++;
				}

				last_index = vindex;
				//coord_estimated = coord.values[vindex];
				//if(flags & UV)
				//	uv_estimated = uv.values[vindex];

				faceorder.push_back(front.size());
				front.emplace_back(current, k, current_edge + prev_(k), current_edge + next_(k));
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
/*
void NxzEncoder::encodeDiff(vector<uchar> &diffs, BitStream &stream, int val) {
	if(val == 0) {
		diffs.push_back(0);
		return;
	}
	int ret = ilog2(abs(val)) + 1;  //0 -> 0, [1,-1] -> 1 [-2,-3,2,3] -> 2
	diffs.push_back(ret);
	int middle = (1<<ret)>>1;
	if(val < 0) val = -val -middle;
	stream.writeUint(val, ret);
}


int needed(int a) {
	if(a == 0) return 0;
	if(a == -1) return 1;
	if(a < 0) a = -a - 1;
	int n = 2;
	while(a >>= 1) n++;
	return n;
}

//0 -> 0
//1 -> [-1, 0]
//2 -> [-2, -1, 0, 1]
//3 -> [-4 -3 -2 0 1 2 3]


void NxzEncoder::encodeDiff(vector<uchar> &diffs, BitStream &bitstream, const Point2i &val) {

	assert(abs(val[0]) < (1<<25));
	assert(abs(val[1]) < (1<<25));

	int diff = std::max(needed(val[0]), needed(val[1]));
	diffs.push_back(diff);
	if(diff == 0) return;

	int max = 1<<(diff-1);
	bitstream.writeUint(val[0] + max, diff);
	bitstream.writeUint(val[1] + max, diff);
}

void NxzEncoder::encodeDiff(vector<uchar> &diffs, BitStream &bitstream, const Point3i &val) {
	assert(abs(val[0]) < (1<<25));
	assert(abs(val[1]) < (1<<25));
	assert(abs(val[2]) < (1<<25));

	int diff = std::max(needed(val[0]), needed(val[1]));
	diff = std::max(diff, needed(val[2]));
	diffs.push_back(diff);
	if(diff == 0) return;

	int max = 1<<(diff-1);
	bitstream.writeUint(val[0] + max, diff);
	bitstream.writeUint(val[1] + max, diff);
	bitstream.writeUint(val[2] + max, diff);
}

*/
