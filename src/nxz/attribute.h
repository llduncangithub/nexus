#ifndef NX_ATTRIBUTE_H
#define NX_ATTRIBUTE_H

#include <map>
#include "nxz.h"
#include "cstream.h"

namespace nx {

class Attribute23 {
public:
	enum Format { UINT32, INT32, UINT16, INT16, UINT8, INT8, FLOAT, DOUBLE };
	enum Strategy { PARALLEL = 0x1, CORRELATED = 0x2 };

	char *buffer;     //output data buffer, input is not needed
	int N;            //number of components
	float q;          //quantization step
	int strategy;

	Format format;    //input or output format
	uint32_t size;    //compressed size (for stats and other nefarious purpouses)


	Attribute23(): buffer(nullptr), N(0), q(0.0f), strategy(0), format(INT32), size(0) {}
	virtual ~Attribute23(){}

	//quantize and store as values
	virtual void quantize(uint32_t nvert, char *buffer) = 0;
	//used by attributes which leverage other attributes
	virtual void preDelta(uint32_t /*nvert*/, std::map<std::string, Attribute23 *> &/*attrs*/, std::vector<uint32_t> &/*index*/) {}
	//use parallelogram prediction or just diff from v0
	virtual void deltaEncode(std::vector<Quad> &context) = 0;
	//compress diffs and write to stream
	virtual void encode(uint32_t nvert, Stream &stream) = 0;

	//read quantized data from stream
	virtual void decode(uint32_t nvert, Stream &stream) = 0;
	//use parallelogram prediction to recover values
	virtual void deltaDecode(uint32_t nvert, std::vector<Face> &faces) = 0;
	//use other attributes to estimate (normals for example)
	virtual void postDelta(uint32_t /*nvert*/, std::map<std::string, Attribute23 *> &/*attrs*/, std::vector<uint32_t> &/*index*/) {}
	//reverse quantization operations
	virtual void dequantize(uint32_t nvert) = 0;
};

//T is int short or char (for colors for example)
template <class T> class GenericAttr: public Attribute23 {
public:
	std::vector<T> values, diffs;

	GenericAttr(int _N) { N = _N; }
	virtual ~GenericAttr(){}


	virtual void quantize(uint32_t nvert, char *buffer) {
		uint32_t n = N*nvert;

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
		case DOUBLE:
			for(uint32_t i = 0; i < n; i++)
				vals[i] = ((double *)buffer)[i]/q;
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
		if(strategy & CORRELATED)
			stream.encodeArray<T>(nvert, &*diffs.begin(), N);
		else
			stream.encodeValues<T>(nvert, &*diffs.begin(), N);

		size = stream.elapsed();
	}

	virtual void decode(uint32_t /*nvert */, Stream &stream) {
		//TODO ensure FORMAT has enough space to store an INT
		int readed;
		if(strategy & CORRELATED)
			readed = stream.decodeArray<T>((T *)buffer, N);
		else
			readed = stream.decodeValues<T>((T *)buffer, N);
	}

	virtual void deltaDecode(uint32_t nvert, std::vector<Face> &context) {
		T *values = (T *)buffer;

		if(strategy & PARALLEL) {
			for(uint32_t i = 1; i < context.size(); i++) {
				Face &f = context[i];
				for(int c = 0; c < N; c++)
					values[i*N + c] += values[f.a*N + c] + values[f.b*N + c] - values[f.c*N + c];
			}
		} else if(context.size()) {
			for(uint32_t i = 1; i < context.size(); i++) {
				Face &f = context[i];
				for(int c = 0; c < N; c++)
					values[i*N + c] += values[f.a*N + c];
			}
		} else { //point clouds assuming values are already sorted by proximity.
			for(uint32_t i = N; i < nvert; i += N)
				values[i] += values[i - N];
		}
	}

	virtual void dequantize(uint32_t nvert) {
		T *coords = (T *)buffer;
		uint32_t n = N*nvert;
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

		case DOUBLE:
			for(uint32_t i = 0; i < n; i++)
				((double *)buffer)[i] = coords[i]*q;
			break;

		case UINT16: //do nothing;
			for(uint32_t i = 0; i < n; i++)
				((uint16_t *)buffer)[i] *= q;
			break;

		case UINT32:
			for(uint32_t i = 0; i < n; i++)
				((uint32_t *)buffer)[i] *= q;
			break;

		case UINT8:
			for(uint32_t i = 0; i < n; i++)
				((char *)buffer)[i] *= q;
			break;
		}
	}
};

}

#endif // NX_ATTRIBUTE_H
