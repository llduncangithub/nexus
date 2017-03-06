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
#include "nxz.h"


namespace nx {

class NxzDecoder {
public:
	uint32_t flags; //keeps track of what is inside
	uint32_t nvert, nface;

	Attribute<Point3i> coord;
	Attribute<Point3i> norm;
	Attribute<int> color[4];
	Attribute<Point2i> uv;
	std::vector<Attribute<int>> data;
	Attribute<int> face;

	Entropy entropy;
	Normals normals;

	CStream stream;

	void NxzDecoder(int len, uchar *input);
	bool hasNormals() { return flags & NORMAL; }
	bool hasColors() { return flags & COLOR; }
	bool hasUvs() { return flags & UV; }
	bool dataCount() { return data.size(); }
	bool hasIndex() { return flags & INDEX; }
	void setCoords(float *buffer);
	void setNormals(float *buffer);
	void setNormals(int16_t *buffer);
	void setColors(uchar *buffer);
	void setUvs(float *buffer);
	void setData(int pos, float *buffer);
	void setIndex(uint32_t *buffer);
	void setIndex(int16_t *buffer);

	void decode();

private:
	bool short_normals;
	bool short_index;
	std::vector<bool> boundary;
	std::vector<int> last;

	void decodeFaces();
	void decodeCoordinates();
	void decodeNormals();
	void decodeColors();

	void shuffle(); //shuffle vertices for point clouds

	void decodeFaces(int start, uint16_t *faces);
	//TODO these are in common with MeshCoder, we should make a NxzEncoder class and move the common parts
	void computeNormals(Point3i *estimated_normals);
	void markBoundary();


	//we assume integer computations and float computations are equivalent for numbers < 1<<23 ? we shouldnt
	int decodeVertex(const vcg::Point3i &predicted, const vcg::Point2i &texpredicted, BitStream &bitstream, int diff, int tdiff);

	int decodeDiff(uchar diff, BitStream &stream);
	void decodeDiff(uchar diff, BitStream &stream, Point3i &p);
	void decodeDiff(uchar diff, BitStream &stream, Point2i &p);

	void dequantize();
	int vertex_count; //keep tracks of current decoding vertex
};

} //namespace
#endif // NXZ_DECODER_H
