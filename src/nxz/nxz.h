#ifndef NXZ_H
#define NXZ_H

#include "cstream.h"

namespace nx {

enum Attributes { INDEX = 1, COORD = 2, NORMAL = 4, COLOR = 8, UV = 16, DATA = 32 };
enum Clers { VERTEX = 0, LEFT = 1, RIGHT = 2, END = 3, BOUNDARY = 4, DELAY = 5, SPLIT = 6};

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
	uint32_t size; //for stats
	std::vector<S> values;
	std::vector<S> diffs;
	void *buffer;
	Attribute(): q(0.0f), buffer(nullptr) {}
	Attribute(float _q): q(_q), size(0), buffer(nullptr) {}
};

/*
void encodeDiff(Stream &stream, uint32_t size, int *values);
void encodeDiff(Stream &stream, uint32_t size, Point2i *values);
void encodeDiff(Stream &stream, uint32_t size, Point3i *values);

int decodeDiff(Stream &stream, int *values);
int decodeDiff(Stream &stream, Point2i *values);
int decodeDiff(Stream &stream, Point3i *values); */

} //namespace

#endif // NXZ_H
