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
#include "index_attribute.h"
#include "vertex_attribute.h"
#include "color_attribute.h"
#include "normal_attribute.h"

namespace nx {

class NxzEncoder {
public:

	uint32_t nvert, nface;


	IndexAttr index;
	std::map<std::string, VertexAttribute *> data;
	int header_size;

	Stream stream;

	NxzEncoder(uint32_t _nvert, uint32_t _nface = 0, Stream::Entropy entropy = Stream::TUNSTALL);

	bool addPositions(float *buffer, float q = 0.0f, Point3f o = Point3f(0.0f));
	bool addPositions(float *buffer, uint32_t *index, float q = 0.0f, Point3f o = Point3f(0.0f));
	bool addPositions(float *buffer, uint16_t *index, float q = 0.0f, Point3f o = Point3f(0.0f));

	bool addNormals(float *buffer, int bits, NormalAttr::Prediction no = NormalAttr::ESTIMATED);
	bool addNormals(int16_t *buffer, int bits, NormalAttr::Prediction no = NormalAttr::ESTIMATED);

	bool addColors(unsigned char *buffer, int lumabits = 6, int chromabits = 6, int alphabits = 5);

	bool addUvs(float *buffer, float q = 0);

	bool addAttribute(const char *name, char *buffer, VertexAttribute::Format format, int components, float q, uint32_t strategy = 0);
	//whatever is inside is your job to fill attr variables.
	bool addAttribute(const char *name, char *buffer, VertexAttribute *attr);

	void addGroup(int end_triangle) { index.groups.push_back(end_triangle); }

	void encode();

private:
	int current_vertex;

	std::vector<bool> boundary;
	std::vector<int> encoded;    //encoded vertex number
	std::vector<Quad> prediction;

	void encodePointCloud();

	void encodeMesh();
	void encodeFaces(int start, int end);
};

} //namespace

#endif // NXZ_ENCODER_H
