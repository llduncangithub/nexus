#ifndef NXZ_H
#define NXZ_H

#include <vector>

namespace nx {

enum Attributes { INDEX = 1, COORD = 2, NORMAL = 4, COLOR = 8, UV = 16, DATA = 32 };
enum Clers { VERTEX = 0, LEFT = 1, RIGHT = 2, END = 3, BOUNDARY = 4, DELAY = 5, SPLIT = 6};
enum Entropy { NONE = 0, TUNSTALL = 1, HUFFMAN = 2, ZLW = 3, LZMA = 4 };
enum Normals { DIFF = 0,      //do not estimate normals, use diffs to previous
			   ESTIMATED = 1, //estimate normals then encode differences
			   BORDER = 2 };  //encode differences only on the boundary


template <typename S> struct Attribute {
	float q; //quantization
	S o;     //origin
	uint32_t size; //for stats
	std::vector<S> values;
	std::vector<S> diffs;
	void *buffer;
	Attribute(): q(0.0f), buffer(nullptr) {}
	Attribute(float _q, S _o): q(_q), o(_o), size(0), buffer(nullptr) {}
};

} //namespace

#endif // NXZ_H
