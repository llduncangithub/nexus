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
#ifndef NXZ_DECODER_H
#define NXZ_DECODER_H

#include <vector>
#include <map>
#include <algorithm>

#include "cstream.h"
#include "tunstall.h"
#include "../nxszip/fpu_precision.h"
#include "zpoint.h"
#include "nxz.h"
#include "attribute.h"
#include "color_attribute.h"
#include "normal_attribute.h"

namespace nx {

class NxzDecoder {
public:
	uint32_t nvert, nface;
//	uint32_t flags;

/*	Attribute<Point3i> coord;
	Attribute<Point3i> norm;
	Attribute<int> color[4];
	Attribute<Point2i> uv; */

	std::map<std::string, Attribute23 *> data;
	Attribute<int> face; //turn this into a pointer: face32 and face16

	Normals normals_prediction;
	Stream stream;

	NxzDecoder(int len, uchar *input);

	bool hasAttr(const char *name) { return data.count(name); }

	bool setPositions(float *buffer) { return setAttribute("position", buffer, Attribute23::FLOAT); }
	bool setNormals(float *buffer)   { return setAttribute("normal", buffer, Attribute23::FLOAT); }
	bool setNormals(int16_t *buffer) { return setAttribute("normal", buffer, Attribute23::INT16); }
	bool setColors(uchar *buffer)    { return setAttribute("color", buffer, Attribute23::INT8); }
	bool setUvs(float *buffer)       { return setAttribute("uv", buffer, Attribute23::FLOAT); }

	bool setAttr(const char *name, char *buffer, Attribute23::Format format) {
		if(data.find(name) == data.end()) return false;
		Attribute23 *attr = data[name];
		attr->format = format;
		attr->buffer = buffer;
		return true;
	}

	bool setAttr(const char *name, char *buffer, Attribute23 *attr) {
		if(data.find(name) == data.end()) return false;
		Attribute23 *found = data[name];
		if(found->id != attr->id)
			return false;
		attr->buffer = buffer;
		delete data[name];
		data[name] = attr;
		return true;
	}

	void setIndex(uint32_t *buffer);
	void setIndex(int16_t *buffer);

/*	template <class T> void setAttribute(int n, T *buffer);
	template <class T> void setAttribute(int n, Attribute *a, T *buffer); */

	void decode();

private:
	bool short_normals;
	bool short_index;
	std::vector<bool> boundary;
	std::vector<uint32_t> groups;
	std::vector<uchar> clers;

	std::vector<Face> prediction;

	int vertex_count; //keep tracks of current decoding vertex

	void decodePointCloud();
	void decodeMesh();
	void decodeZPoints();
	void decodeCoords();
	void decodeNormals();
	void decodeColors();
	void decodeUvs();
	void decodeDatas();

	void shuffle(); //shuffle vertices for point clouds

	void decodeFaces(uint32_t start, uint32_t end, uint32_t &cler, BitStream &bitstream);
	//TODO these are in common with MeshCoder, we should make a NxzEncoder class and move the common parts
	void estimateNormals(Point3f *normals3f);
	void computeNormals(Point3f *normals, Point3f *estimated);
	void computeNormals(Point3s *normals, Point3f *estimated);

	template <class F> void markBoundary();


	//we assume integer computations and float computations are equivalent for numbers < 1<<23 ? we shouldnt
//	int decodeVertex(const Point3i &predicted, const Point2i &texpredicted, int last_index);
	int decodeVertex(int v0, int v1, int v2);


	int decodeDiff(uchar diff, BitStream &stream);
	void decodeDiff(uchar diff, BitStream &stream, Point3i &p);
	void decodeDiff(uchar diff, BitStream &stream, Point3s &p);
	void decodeDiff(uchar diff, BitStream &stream, Point2i &p);
	void decodeDiff(uchar diff, BitStream &stream, Point2s &p);

	void dequantize();
};

template <class F> void NxzDecoder::markBoundary() {
	boundary.resize(nvert, false);

	std::vector<int> count(nvert, 0);
	F *f = (F *)face.buffer;
	for(uint32_t i = 0; i < nface; i++) {
		count[f[0]] += (int)f[1] - (int)f[2];
		count[f[1]] += (int)f[2] - (int)f[0];
		count[f[2]] += (int)f[0] - (int)f[1];
		f += 3;
	}
	for(uint32_t i = 0; i < nvert; i++)
		if(count[i] != 0)
			boundary[i] = true;
}



} //namespace
#endif // NXZ_DECODER_H
