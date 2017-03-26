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

	std::map<std::string, Attribute23 *> data;
	Attribute<int> face; //turn this into a different class: with face32 and face16 pointers.


	NxzDecoder(int len, uchar *input);

	bool hasAttr(const char *name) { return data.count(name); }

	bool setPositions(float *buffer) { return setAttribute("position", (char *)buffer, Attribute23::FLOAT); }
	bool setNormals(float *buffer)   { return setAttribute("normal", (char *)buffer, Attribute23::FLOAT); }
	bool setNormals(int16_t *buffer) { return setAttribute("normal", (char *)buffer, Attribute23::INT16); }
	bool setColors(uchar *buffer)    { return setAttribute("color", (char *)buffer, Attribute23::UINT8); }
	bool setUvs(float *buffer)       { return setAttribute("uv", (char *)buffer, Attribute23::FLOAT); }

	bool setAttribute(const char *name, char *buffer, Attribute23::Format format);
	bool setAttribute(const char *name, char *buffer, Attribute23 *attr);

	void setIndex(uint32_t *buffer) {face.buffer = buffer; }
	void setIndex(int16_t *buffer) {
		face.buffer = buffer;
		short_index = true;
	}

	void decode();

private:
	Stream stream;

	bool short_normals;
	bool short_index;

	std::vector<uint32_t> groups;
	std::vector<uchar> clers;
	std::vector<Face> prediction;

	int vertex_count; //keep tracks of current decoding vertex

	void decodePointCloud();
	void decodeMesh();
	void decodeFaces(uint32_t start, uint32_t end, uint32_t &cler, BitStream &bitstream);
};


} //namespace
#endif // NXZ_DECODER_H
