#ifndef NXZ_H
#define NXZ_H

#include <vector>
#include <map>
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

/*
void encodeDiff(Stream &stream, uint32_t size, int *values);
void encodeDiff(Stream &stream, uint32_t size, Point2i *values);
void encodeDiff(Stream &stream, uint32_t size, Point3i *values);

int decodeDiff(Stream &stream, int *values);
int decodeDiff(Stream &stream, Point2i *values);
int decodeDiff(Stream &stream, Point3i *values); */


class Attribute23 {
public:
	enum Format { INT32, INT16, INT8, FLOAT };
	enum Strategy { PARALLEL = 0x1, CORRELATED = 0x2, CUSTOM = 0x4 };

	char *buffer;     //output data buffer, input is not needed

	float q;          //quantization step
	int components;   //number of components
	uchar strategy;
	Format format;    //input or output format
	uint32_t size;    //compressed size

	Attribute23(float _q = 0.0f): buffer(nullptr), q(_q), components(1), strategy(0), format(FLOAT), size(0) {}

	//quantize and store as values
	virtual void quantize(uint32_t nvert, char *buffer) = 0;
	virtual void preDelta(std::map<std::string, Attribute23 *> &attrs, std::vector<uint32_t> &index) {}
	//use parallelogram prediction or just diff from v0
	virtual void deltaEncode(std::vector<Quad> &context) = 0;
	//compress diffs and write to stream
	virtual void encode(uint32_t nvert, Stream &stream) = 0;

	//read quantized data from stream
	virtual void decode(uint32_t nvert, Stream &stream) = 0;
	//use parallelogram prediction to recover values
	virtual void deltaDecode(std::vector<Face> &faces) = 0;
	//use other attributes to estimate (normals for example)
	virtual void postDelta(std::map<std::string, Attribute23 *> &attrs, std::vector<uint32_t> &index) {}
	//reverse quantization operations
	virtual void dequantize(uint32_t nvert) = 0;
};

//T internal format, S external format
template <class T, int N> class Data: public Attribute23 {
public:
	std::vector<T> values, diffs;

	Data(float q = 0.0f): Attribute23(q) {}

	virtual void quantize(uint32_t nvert, char *buffer) {
		uint32_t n = components*nvert;

		values.resize(n);
		diffs.resize(n);
		int *vals = (int *)&*values.begin();
		switch(format) {
		case INT32:
			for(uint32_t i = 0; i < n; i++)
				vals[i] = ((int32_t *)buffer)[i]/q;
			break;
		case INT16:
			for(uint32_t i = 0; i < n; i++)
				vals[i] = ((int16_t *)buffer)[i]/q;
			break;
		case INT8:
			for(uint32_t i = 0; i < n; i++)
				vals[i] = ((int16_t *)buffer)[i]/q;
			break;
		case FLOAT:
			for(uint32_t i = 0; i < n; i++)
				vals[i] = ((float *)buffer)[i]/q;
			break;
		}

	}

	virtual void deltaEncode(std::vector<Quad> &context) {
		for(int c = 0; c < N; c++)
			diffs[c] = values[context[0].t*N + c];
		for(uint32_t i = 1; i < context.size(); i++) {
			Quad &q = context[i];
			if(q.a != q.b && (strategy & PARALLEL)) {
				for(int c = 0; c < N; c++)
					diffs[i*N + c] = values[q.t*N + c] - (values[q.a*N + c] + values[q.b*N + c] - values[q.c*N + c]);
			} else {
				for(int c = 0; c < N; c++)
					diffs[i*N + c] = values[q.t*N + c] - values[q.a*N + c];
			}
		}
	}

	virtual void encode(uint32_t nvert, Stream &stream) {
		stream.restart();
		stream.write<float>(q);
		//encodeDiff(stream, diffs.size(), &*diffs.begin());
		if(strategy & CORRELATED)
			stream.encodeArray<T, N>(nvert, &*diffs.begin());
		else
			stream.encodeValues<T, N>(nvert, &*diffs.begin());

		size = stream.elapsed();
	}

	virtual void decode(uint32_t nvert, Stream &stream) {
		q = stream.read<float>();
		T *coords = (T *)buffer;
		int readed;
		if(strategy & CORRELATED)
			readed = stream.decodeArray<T, N>(coords);
		else
			readed = stream.decodeValues<T, N>(coords);
	}

	virtual void deltaDecode(std::vector<Face> &context) {
		T *values = (T *)buffer;

		if(strategy & PARALLEL) {
			for(uint32_t i = 1; i < context.size(); i++) {
				Face &f = context[i];
				for(int c = 0; c < N; c++)
					values[i*N + c] += values[f.a*N + c] + values[f.b*N + c] - values[f.c*N + c];
			}
		} else {
			for(uint32_t i = 1; i < context.size(); i++) {
				Face &f = context[i];
				for(int c = 0; c < N; c++)
					values[i*N + c] += values[f.a*N + c];
			}
		}
	}

	virtual void dequantize(uint32_t nvert) {
		T *coords = (T *)buffer;
		uint32_t n = nvert*components;
		switch(format) {
		case FLOAT:
			for(uint32_t i = 0; i < n; i++)
				((float *)buffer)[i] = coords[i]*q;
			break;

		case INT16: //do nothing;
			for(uint32_t i = 0; i < n; i++)
				((uint16_t *)buffer)[i] *= q;
			break;

		case INT32:
			for(uint32_t i = 0; i < n; i++)
				((uint32_t *)buffer)[i] *= q;
			break;

		case INT8:
			for(uint32_t i = 0; i < n; i++)
				((char *)buffer)[i] *= q;
			break;
		}
	}
};

class Position: public Data<int, 3> {
	Position(float _q): Data<int, 3>(_q) {}
};

class Normal1: public Data<int, 2> {
public:
	Normal1(float _q): Data<int, 2>(_q) {}
};

class Color: public Data<uchar, 4> {
public:
	Color(float _q): Data<uchar, 4>(_q) {}
};

class Uv: public Data<int, 2> {
public:
	Uv(float _q): Data<int, 2>(_q) {}
};

} //namespace

#endif // NXZ_H
