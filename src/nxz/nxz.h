#ifndef NXZ_H
#define NXZ_H

#include <vector>
#include "point.h"
#include "cstream.h"

namespace nx {

enum Attributes { INDEX = 1, COORD = 2, NORMAL = 4, COLOR = 8, UV = 16, DATA = 32 };
enum Clers { VERTEX = 0, LEFT = 1, RIGHT = 2, END = 3, BOUNDARY = 4, DELAY = 5, SPLIT = 6};
enum Normals { DIFF = 0,      //do not estimate normals, use diffs to previous
			   ESTIMATED = 1, //estimate normals then encode differences
			   BORDER = 2 };  //encode differences only on the boundary


struct Face {
	uint32_t a, b, c;
	Face() {}
	Face(uint32_t v0, uint32_t v1, uint32_t v2): a(v0), b(v1), c(v2) {}
};

struct Quad {
	uint32_t t, a, b, c;
	Quad() {}
	Quad(uint32_t _t, uint32_t v0, uint32_t v1, uint32_t v2): t(_t), a(v0), b(v1), c(v2) {}
};

template <typename S> struct Attribute {
	float q;       //quantization
	S o;           //origin
	uint32_t size; //for stats
	std::vector<S> values;
	std::vector<S> diffs;
	void *buffer;
	Attribute(): q(0.0f), buffer(nullptr) {}
	Attribute(float _q, S _o): q(_q), o(_o), size(0), buffer(nullptr) {}
};



/*
//in -> internal representation Point3i, OUT external representation (Point3f)
template <typename IN, typename OUT> struct Attribute1 {
	float q;       //quantization
	IN o;          //origin
	uint32_t size; //for stats
	OUT *buffer;
	bool predicted;

	std::vector<IN> values;
	std::vector<IN> diffs;

	Attribute1(): q(0.0f), buffer(nullptr) {}
	Attribute1(float _q, S _o, bool _predicted = false): q(_q), o(_o), size(0), buffer(nullptr), predicted(_predicted) {}

	virtual void add(OUT *buffer) = 0;
	virtual void encode(Stream &stream) {
		std::vector<uchar> diffs;
		BitStream bitstream(diffs.size());

		for(auto &u: diffs)
			encodeDiff(diffs, bitstream, u);

		Tunstall::compress(stream, &*diffs.begin(), diffs.size());
		stream.write(bitstream);
		size = stream.elapsed();
	}

	virtual void diff(int diffindex, int valueindex, int lastindex) {
		diffs[diffindex] = values[valueindex];
		if(last != -1) {
			if(predicted)
				diffs[diffindex] = values[valueindex] - (values[v0] + values[v1] - values[v2]);
			else
				diffs[diffindex] -= diffs[lastindex];
		}
	}


};

class ColordAttribute: public Attribute<Point3i> {

};
*/

} //namespace

#endif // NXZ_H
