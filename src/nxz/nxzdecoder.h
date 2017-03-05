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
#include <algorithm>

#include "cstream.h"
#include "tunstall.h"
#include "../nxszip/fpu_precision.h"
#include "zpoint.h"


namespace nx {

class NxzDecoder {
public:
	enum Attributes { INDEX = 0, COORD = 1, NORMAL = 2, COLOR = 4, UV = 8, DATA = 16 };
	enum Clers { VERTEX = 0, LEFT = 1, RIGHT = 2, END = 3, BOUNDARY = 4, DELAY = 5 };

	template <typename S> struct Attribute {
		float q; //quantization
		S o;     //origin
		int bits;
		uint32_t size; //for stats
		Attribute(): q(0.0f), bits(0) {}
		Attribute(float _q, S _o, int _b): q(_q), o(_o), bits(_b), size(0) {}
	};


	Attribute<Point3i> coord;
	Attribute<Point3i> norm;
	Attribute<int> color[4];
	Attribute<Point2i> uv;
	std::vector<Attribute<int>> data;
	Attribute<int> face;

	CStream stream;

	NxzDecoder(): vertex_count(0) {}

	void decode(int len, uchar *input);

private:
	std::vector<bool> boundary;
	std::vector<int> last;

	void decodeFaces();
	void decodeCoordinates();
	void decodeNormals();
	void decodeColors();

	void shuffle(); //shuffle vertices for point clouds

	void decodeFaces(int start, uint16_t *faces);
	//TODO these are in common with MeshCoder, we should make a NxzEncoder class and move the common parts
	void computeNormals(vcg::Point3s *estimated_normals);
	void markBoundary();

	int decodeDiff(uchar diff, BitStream &stream);
	//we assume integer computations and float computations are equivalent for numbers < 1<<23
	int decodeVertex(const vcg::Point3i &predicted, const vcg::Point2i &texpredicted, BitStream &bitstream, int diff, int tdiff);

	void dequantize();
	int vertex_count; //keep tracks of current decoding vertex
};

} //namespace
#endif // NXZ_DECODER_H
