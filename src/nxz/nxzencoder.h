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
#ifndef NXZ_ENCODER_H
#define NXZ_ENCODER_H

#include <vector>
#include <map>

#include <limits.h>
#include <float.h>

#include "cstream.h"
#include "zpoint.h"
#include "nxz.h"
#include "attribute.h"
#include "color_attribute.h"
#include "normal_attribute.h"

namespace nx {

class NxzEncoder {
public:
	//uint32_t flags; //keeps track of what is inside
	uint32_t nvert, nface;

/*	Attribute<Point3i> coord;
	Attribute<Point2i> norm;
	Attribute<unsigned char> color[4];
	Attribute<Point2i> uv; */
	std::map<std::string, Attribute23 *> data;

	std::vector<uint32_t> index;
	int index_size;
	int header_size;

//	Normals normals_prediction;

	Stream stream;

	NxzEncoder(uint32_t _nvert, uint32_t _nface = 0, Stream::Entropy entropy = Stream::TUNSTALL);

	bool addPositions(float *buffer, float q = 0.0f, Point3f o = Point3f(0.0f));
	bool addPositions(float *buffer, uint32_t *index, float q = 0.0f, Point3f o = Point3f(0.0f));
	bool addPositions(float *buffer, uint16_t *index, float q = 0.0f, Point3f o = Point3f(0.0f));

	bool addNormals(float *buffer, int bits, NormalAttr::Prediction no = NormalAttr::ESTIMATED);
	bool addNormals(int16_t *buffer, int bits, NormalAttr::Prediction no = NormalAttr::ESTIMATED);

	bool addColors(unsigned char *buffer, int lumabits = 6, int chromabits = 6, int alphabits = 5);

	bool addUV(float *buffer, float q = 0);

	bool addAttribute(const char *name, char *buffer, Attribute23::Format format, int components, float q, uint32_t strategy) {
		if(data.count(name)) return false;
		GenericAttr<int> *attr = new GenericAttr<int>(components);

		attr->q = q;
		attr->strategy = strategy;
		attr->format = format;
		attr->quantize(nvert, (char *)buffer);
		data[name] = attr;
		return true;
	}
	//whatever is inside is your job to fill attr variables.
	bool addAttribute(const char *name, char *buffer, Attribute23 *attr) {
		if(data.count(name)) return true;
		attr->quantize(nvert, buffer);
		data[name] = attr;
		return true;
	}

	void addGroup(int end_triangle) { groups.push_back(end_triangle); }

	void encode();

private:
	std::vector<uint32_t> groups;
	std::vector<uchar> clers;
	int current_vertex;

	std::vector<bool> boundary;
	std::vector<int> encoded;    //encoded vertex number
	std::vector<Quad> prediction;

//	void setCoordBits();
//	void setDataBits();

	void encodePointCloud();
//	void encodeZPoints(std::vector<ZPoint> &zpoints);

	void encodeMesh();
	void encodeFaces(int start, int end, BitStream &bitstream);

/*	void encodeCoords();
	void encodeNormals();
	void encodeColors();
	void encodeUvs();
	void encodeDatas(); */

/*	void computeNormals(std::vector<Point2i> &estimated_normals);
	void markBoundary();

	void encodeDiff(std::vector<uchar> &diffs, BitStream &stream, int val);
	void encodeDiff(std::vector<uchar> &diffs, BitStream &stream, const Point2i &val);
	void encodeDiff(std::vector<uchar> &diffs, BitStream &stream, const Point3i &val); */
};

} //namespace

#endif // NXZ_ENCODER_H
