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

#include <limits.h>
#include <float.h>
//#include "point.h"

#include "cstream.h"
#include "zpoint.h"

namespace nx {

class NxzEncoder {
public:
	enum Attributes { INDEX = 0, COORD, NORMAL, COLOR, UV, DATA };
	enum Clers { VERTEX = 0, LEFT = 1, RIGHT = 2, END = 3, BOUNDARY = 4, DELAY = 5 };

	uint32_t flags; //keeps track of what is inside
	uint32_t nvert, nface;

	CStream stream;
	//compression stats
	int coord_size;
	int normal_size;
	int color_size;
	int face_size;
	int uv_size;
	std::vector<int> data_size;

	int coord_bits;    //number of bits for coordinates.
	int uv_bits;       //bumber of bits for texture coordinates
	int norm_bits;     //normal bits
	int color_bits[4]; //color bits
	std::vector<int> data_bits;
	int index_bits;


	NxzEncoder(uint32_t _nvert, uint32_t _nface = 0);
	void addCoords(float *buffer, float q = 0, Point3f o = Point3f(FLT_MAX));
	void addCoords(float *buffer, uint32_t *index, float q = 0, Point3f o = Point3f(FLT_MAX));
	void addCoords(float *buffer, uint16_t *index, float q = 0, Point3f o = Point3f(FLT_MAX));
	//you do the quantization step
	void addCoords(int *buffer);
	void addCoords(int *buffer, uint32_t *index);
	void addCoords(int *buffer, uint16_t *index);

	void addNormals(float *buffer, int bits);
	void addNormals(int16_t *buffer, int bits);

	void addColors(unsigned char *buffer, int lumabits = 6, int chromabits = 6, int alphabits = 5);

	void addUV(float *buffer, float q = 0);
	void addUV(int *buffer, float q = 0);

	void addData(float *buffer, float q, float offset = 0);
	void addData(int *buffer);

	void addGroup(int end_triangle) { groups.push_back(end_triangle); }

	void encode();

private:

	Point3i coord_o;
	Point2i uv_o;
	std::vector<float> data_o;

	float coord_q;    //coordinates quantization.
	float norm_q;     //expecting a power of 2.
	float color_q[4]; //expecting a power of 2.
	float uv_q;       //texture quantization
	std::vector<float> data_q;

	std::vector<Point3i> qcoords;
	std::vector<Point2i> qtexcoords;
	std::vector<Point3i> qnormals;
	std::vector<Color4b> qcolors;
	std::vector<Point2i> quvs;
	std::vector<std::vector<int>> qdatas;

	std::vector<uint32_t> index;
	std::vector<uint32_t> groups;
	vector<uchar> clers;

	std::vector<uchar> dcoords;
	std::vector<uchar> dnormals;
	std::vector<uchar> dcolors[4];
	std::vector<uchar> duvs;
	std::vector<std::vector<uchar>> ddatas;

	std::vector<ZPoint> zpoints; //used by point cloud only
	std::vector<int> order;      //for mesh we store for each new vertex the original vertex index.

	std::vector<int> last;       //used with order to make diffs in colors (refers to the original indexes too.
	std::vector<bool> boundary;
	std::vector<int> encoded;    //encoded vertices, use index instead of diff for those.

	void setCoordBits();
	void setDataBits();

	void encodePointCloud(); //used only for point clouds
	void encodeCoords();
	void encodeNormals();
	void encodeColors();
	void encodeUvs();
	void encodeDatas();

	void encodeMesh();
	void encodeFaces(BitStream &bitstream, int start, int end);

	void computeNormals(vector<Point3i> &estimated_normals);
	void markBoundary();

	void encodeVertex(int target, const Point3i &predicted, const Point2i &texpredicted, BitStream &bitstream, int last);
	void encodeDiff(std::vector<uchar> &diffs, BitStream &stream, int val);
};

} //namespace

#endif // NXZ_ENCODER_H
