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
	uint32_t size; //for stats
	std::vector<S> values;
	std::vector<S> diffs;
	void *buffer;
	Attribute(): q(0.0f), buffer(nullptr) {}
	Attribute(float _q): q(_q), size(0), buffer(nullptr) {}
};

void encodeDiffs(Stream &stream, uint32_t size, int *values);
void encodeDiffs(Stream &stream, uint32_t size, short *values);
void encodeDiffs(Stream &stream, uint32_t size, Point2i *values);
void encodeDiffs(Stream &stream, uint32_t size, Point2s *values);
void encodeDiffs(Stream &stream, uint32_t size, Point3i *values);
void encodeDiffs(Stream &stream, uint32_t size, Point3s *values);

int decodeDiff(Stream &stream, int *values);
int decodeDiff(Stream &stream, short *values);
int decodeDiff(Stream &stream, Point2i *values);
int decodeDiff(Stream &stream, Point2s *values);
int decodeDiff(Stream &stream, Point3i *values);
int decodeDiff(Stream &stream, Point3s *values);


class Attribute23 {
public:
	float q;         //quantization step
	char *buffer;    //input data buffer
	uint32_t size;   //compressed size

	Attribute23(float _q = 0.0f): q(_q), buffer(nullptr), size(0) {}
	//quantize and store as values
	virtual void quantize(uint32_t nvert, char *buffer) = 0;
	//use parallelogram prediction or just diff from v0
	virtual void diff(std::vector<Quad> &context, std::vector<Point3i> &coords, std::vector<int> &faces) = 0;
	//compress diffs and write to stream
	virtual void encode(Stream &stream) = 0;

	//read quantized data from stream
	virtual void decode(uint32_t nvert, Stream &stream) = 0;
	//use parallelogram prediction to recover values
	virtual void dediff(std::vector<Point3i> &coords, std::vector<int> &index, std::vector<Face> &faces) = 0;
	//reverse quantization operations
	virtual void dequantize(uint32_t nvert) = 0;
};

//T internal format, S external format
template <class T, class S> class Data: public Attribute23 {
public:
	std::vector<T> values, diffs;
	virtual void quantize(uint32_t nvert, S *buffer) {
		values.resize(nvert);
		diffs.resize(nvert);
		for(uint32_t i = 0; i < nvert; i++)
			values[i] = buffer[i]/q;
	}

	virtual void diff(std::vector<Quad> &context, std::vector<Point3i> &/*coords*/, std::vector<int> &/*faces*/) {
		diffs[0] = values[context[0].t];
		for(uint32_t i = 1; i < context.size(); i++) {
			Quad &q = context[i];
			diffs[i] = values[q.t] - (values[q.a] + values[q.b] - values[q.c]);
		}
	}

	virtual void encode(Stream &stream) {
		stream.restart();
		encodeDiff(stream, diffs);
		size = stream.elapsed();
	}

	virtual void decode(uint32_t nvert, Stream &stream) {
		T *coords = (T *)buffer;
		int readed = decodeDiff(stream, coords);
	}

	virtual void dediff(std::vector<Point3i> &/*coords*/, std::vector<int> &/*index*/, std::vector<Face> &context) {
		T *values = (T *)buffer;
		for(uint32_t i = 1; i < context.size(); i++) {
			Face &f = context[i];
			values[i] += values[f.a] + values[f.b] - values[f.c];
		}
	}

	virtual void dequantize(uint32_t nvert) {
		T *coords = (T *)buffer;
		S *points = (S *)buffer;
		for(uint32_t i = 0; i < nvert; i++) {
			T &v = coords[i];
			S &p = points[i];
			p = ((S)v)*q;
		}
	}
};

} //namespace

#endif // NXZ_H
