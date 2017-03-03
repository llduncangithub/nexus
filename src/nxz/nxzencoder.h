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

#include <vcg/space/point3.h>
#include <vcg/space/point2.h>

#include "../nxszip/cstream.h"
#include "../nxszip/zpoint.h"


namespace nx {

//face topology utility structures.


class MeshEncoder {
public:
	enum Clers { VERTEX = 0, LEFT = 1, RIGHT = 2, END = 3, BOUNDARY = 4, DELAY = 5 };

	int coord_q; //global power of 2 quantization.
	int tex_q;   //texture bits

	int coord_bits; //number of bits for coordinates.
	int tex_bits; //bumber of bits for texture coordinates
	int norm_bits;  //normal bits
	int color_bits[4]; //color bits

	CStream stream;
	//compression stats
	int coord_size;
	int normal_size;
	int color_size;
	int face_size;

	MeshEncoder(uint32_t nvert, uint32_t nface):
		coord_q(0), tex_q(0), coord_bits(12), tex_bits(12), norm_bits(10),
		coord_size(0), normal_size(0), color_size(0), face_size(0) {
		color_bits[0] = color_bits[1] = color_bits[2] = color_bits[3] = 6;
	}
	void addCoords(float *buffer, float q = 0); //if not specified, use coord_bits instead
	void addNormals(float *buffer, int bits);
	void addColors(unsigned char *buffer, int lumabits, int chromabits, int alphabits);
	void addUV(float *buffer, float q =  0);
	void addData(float *buffer, float q);
	void encode();

private:
	float *coords, *normals, *uv;
	float *data[8];

	vcg::Point3i min, max; //minimum position of vertices
	vcg::Point2i tmin, tmax; //minimum position texture coordinates


	std::vector<vcg::Point3i> qpoints;
	std::vector<vcg::Point2i> qtexcoords;

	std::vector<ZPoint> zpoints; //used by point cloud only
	std::vector<int> order; //not used for point cloud for mesh we store for each new vertex the original vertex index.
	std::vector<int> reorder; //for each OLD vertex its new vertex

	std::vector<int> last; //used with order to make diffs in colors (refers to the original indexes too.
	std::vector<bool> boundary;
	std::vector<int> encoded; //encoded vertices, use index instead of diff for those.

	static vcg::Point3i quantize(vcg::Point3f &p, float side); //used in encode faces
	void quantize(); //used for point clouds
	void quantizeCoords();
	void quantizeTexCoords();

	void encodeFaces();
	void encodeCoordinates(); //used only for point clouds
	void encodeNormals();
	void encodeColors();
	void encodeTexCoords();

	void encodeFaces(int start, int end);

	void computeNormals(std::vector<vcg::Point3s> &estimated_normals);
	void markBoundary();

	void encodeVertex(int target, const vcg::Point3i &predicted, const vcg::Point2i &texpredicted, BitStream &bitstream, vector<uchar> &diffs, vector<uchar> &tdiffs);
	void encodeDiff(std::vector<uchar> &diffs, BitStream &stream, int val);
};

} //namespace

#endif // NXZ_ENCODER_H
