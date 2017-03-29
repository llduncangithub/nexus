#ifndef NX_COLOR_ATTRIBUTE_H
#define NX_COLOR_ATTRIBUTE_H

#include "attribute.h"
#include "point.h"

namespace nx {

class ColorAttr: public GenericAttr<uchar> {
public:
	int qc[4];
	ColorAttr(): GenericAttr<uchar>(4) {
		qc[0] = qc[1] = qc[2] = 4;
		qc[3] = 8;
	}

	void setQ(int lumabits, int chromabits, int alphabits) {
		qc[0] = 1<<(8 - lumabits);
		qc[1] = qc[2] = 1<<(8 - chromabits);
		qc[3] = 1<<(8-alphabits);
	}

	virtual void quantize(uint32_t nvert, char *buffer);
	virtual void dequantize(uint32_t nvert);

	virtual void encode(uint32_t nvert, Stream &stream) {
		stream.restart();
		for(int c = 0; c < 4; c++)
			stream.write<uchar>(qc[c]);

		stream.encodeValues<char>(nvert, (char *)&*diffs.begin(), N);
		size = stream.elapsed();
	}
	virtual void decode(uint32_t /*nvert*/, Stream &stream) {
		for(int c = 0; c < 4; c++)
			qc[c] = stream.read<uchar>();
		stream.decodeValues<uchar>((uchar *)buffer, N);
	}
};

} //namespace

#endif // NX_COLOR_ATTRIBUTE_H
