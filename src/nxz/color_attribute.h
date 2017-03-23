#ifndef NX_COLOR_ATTRIBUTE_H
#define NX_COLOR_ATTRIBUTE_H

#include "attribute.h"
#include "point.h"

namespace nx {

class ColorAttr: public GenericAttr<uchar> {
public:
	int qc[4];
	ColorAttr(): GenericAttr<uchar>(4) {
		qc[0] = qc[1] = qc[2] = 6;
		qc[3] = 5;
	}

	void setQ(int lumabits, int chromabits, int alphabits) {
		qc[0] = 1<<(8 - lumabits);
		qc[1] = qc[2] = 1<<(8 - chromabits);
		qc[3] = 1<<(8-alphabits);
	}

	virtual void quantize(uint32_t nvert, char *buffer) {
		uint32_t n = N*nvert;

		values.resize(n);
		diffs.resize(n);
		switch(format) {
		case INT8:
		{
			Color4b *colors = (Color4b *)buffer;
			for(uint32_t i = 0; i < nvert; i++) {
				Color4b c = colors[i];
				Color4b &y = *(Color4b *)&(values[i]);
				for(int k = 0; k < 4; k++)
					y[k] = c[k]/qc[k];
				y = y.toYCC();
			}
		}
		break;

		case FLOAT:
		{
			float *colors = (float *)buffer;
			for(uint32_t i = 0; i < nvert; i++) {
				Color4b &y = *(Color4b *)&(values[i]);
				for(int k = 0; k < 4; k++)
					y[k] = ((int)(colors[i*N + k]*255.0f))/qc[k];
				y = y.toYCC();
			}
		}
		break;

		default: throw "Unsupported color input format.";
		}
	}

	virtual void dequantize(uint32_t nvert) {
		switch(format) {
		case INT8:
		{
			Color4b *colors = (Color4b *)buffer;
			for(uint32_t i = 0; i < nvert; i++) {
				Color4b &c = colors[i];
				c = c.toRGB();
				for(int k = 0; k < 4; k++)
					c[k] *= qc[k];
			}
			break;
		}
		case FLOAT:
		{
			std::vector<Color4b> colors(nvert);
			memcpy(&*colors.begin(), buffer, nvert*sizeof(Color4b));
			for(uint32_t i = 0; i < nvert; i++) {
				Color4b &c = colors[i];
				c = c.toRGB();
				for(int k = 0; k < 4; k++)
					((float *)buffer)[i*4 +k] = (c[k]*qc[k])/255.0f;
			}
			break;
		}
		default: throw "Unsupported color output format.";
		}
	}
};

} //namespace

#endif // NX_COLOR_ATTRIBUTE_H
