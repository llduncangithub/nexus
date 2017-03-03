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

#include <vcg/space/color4.h>
#include <vcg/space/point3.h>
#include <vcg/space/point2.h>

#include "../nxszip/cstream.h"
#include "../nxszip/zpoint.h"

namespace nx {

class NxzEncoder {
public:
	enum Clers { VERTEX = 0, LEFT = 1, RIGHT = 2, END = 3, BOUNDARY = 4, DELAY = 5 };


	CStream stream;
	//compression stats
	int coord_size;
	int normal_size;
	int color_size;
	int face_size;

	NxzEncoder(uint32_t _nvert, uint32_t _nface):
		nvert(_nvert), nface(_nface),
		coord_q(0), uv_q(0), coord_bits(12), uv_bits(12), norm_bits(10),
		coord_size(0), normal_size(0), color_size(0), face_size(0),
		coords(NULL), normals(NULL), uv(NULL), colors(NULL) {
		color_bits[0] = color_bits[1] = color_bits[2] = color_bits[3] = 6;
		for(int i = 0; i < 8; i++) data[i] = NULL;
	}
	void addCoords(float *buffer, float q = 0, vcg::Point3f o = vcg::Point3f(0, 0, 0)) {
		coords = (vcg::Point3f *)buffer;
		coord_q = q;
		offset = o;
	}

	void addNormals(float *buffer, int bits) {
		normals = (vcg::Point3f *)buffer;
		norm_bits = bits;
	}

	void addNormals(int16_t *buffer, int bits) {
		normals16 = (vcg::Point2s *)buffer;
		norm_bits = bits;
	}

	void addColors(unsigned char *buffer, int lumabits, int chromabits, int alphabits) {
		colors = (vcg::Color4b *)buffer;
		color_bits[0] = lumabits;
		color_bits[1] = color_bits[2] = chromabits;
		color_bits[3] = alphabits;
	}
	void addUV(float *buffer, float q =  0) {
		uv = (vcg::Point2f *)buffer;
		uv_q = q;
	}

	void addData(float *buffer, float q) {
		data[0] = buffer;
		data_q[0] = q;
	}

	void encode();

private:
	int nvert, nface;
	int coord_q; //coordinates quantization.
	vcg::Point3f offset;
	int uv_q;    //texture quantization
	int data_q[8];

	int coord_bits;    //number of bits for coordinates.
	int uv_bits;       //bumber of bits for texture coordinates
	int norm_bits;     //normal bits
	int color_bits[4]; //color bits

	vcg::Point3f *coords, *normals;
	vcg::Point2f *uv;
	vcg::Point2s *normals16;
	vcg::Color4b *colors;
	float *data[8];

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

	void encodeFaces(int start, int end);

	void computeNormals(std::vector<vcg::Point3s> &estimated_normals);
	void markBoundary();

	void encodeVertex(int target, const vcg::Point3i &predicted, const vcg::Point2i &texpredicted, BitStream &bitstream, vector<uchar> &diffs, vector<uchar> &tdiffs);
	void encodeDiff(std::vector<uchar> &diffs, BitStream &stream, int val);
};

} //namespace

#endif // NXZ_ENCODER_H
