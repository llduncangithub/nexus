#ifndef NX_NORMAL_ATTRIBUTE_H
#define NX_NORMAL_ATTRIBUTE_H

#include "attribute.h"
#include "point.h"

namespace nx {

class NormalAttr: public GenericAttr<int> {
public:
	enum Prediction { DIFF = 0x0,      //do not estimate normals, use diffs to previous
					  ESTIMATED = 0x1, //estimate normals then encode differences
					  BORDER = 0x2 };  //encode differences only on the boundary
	uint32_t prediction;
	std::vector<bool> boundary;

	NormalAttr(int bits = 10): GenericAttr<int>(2) {
		q = pow(2, bits-1);
		prediction = DIFF;
		strategy |= Attribute23::CORRELATED;
	}

	virtual void quantize(uint32_t nvert, char *buffer) {
		uint32_t n = N*nvert;

		values.resize(n);
		diffs.resize(n);

		Point2i *normals = (Point2i *)&*values.begin();
		switch(format) {
		case FLOAT:
			for(uint32_t i = 0; i < nvert; i++)
				normals[i] = toOcta(((Point3f *)buffer)[i], (int)q);
			break;
		case INT32:
			for(uint32_t i = 0; i < nvert; i++)
				normals[i] =  toOcta(((Point3i *)buffer)[i], (int)q);
			break;
		case INT16:
		{
			Point3<int16_t> *s = (Point3<int16_t> *)buffer;
			for(uint32_t i = 0; i < nvert; i++)
				normals[i] = toOcta(Point3i(s[i][0], s[i][1], s[i][2]), (int)q);
			break;
		}
		case INT8:
		{
			Point3<int8_t> *s = (Point3<int8_t> *)buffer;
			for(uint32_t i = 0; i < nvert; i++)
				normals[i] = toOcta(Point3i(s[i][0], s[i][1], s[i][2]), (int)q);
			break;
		}
		default: throw "Unsigned types not supported for normals";
		}
	}


	virtual void preDelta(uint32_t nvert, std::map<std::string, Attribute23 *> &attrs, std::vector<uint32_t> &index) {
		if(prediction == DIFF)
			return;

		if(attrs.find("position") == attrs.end())
			throw "No position attribute found. Use DIFF normal strategy instead.";

		GenericAttr<int> *coord = dynamic_cast<GenericAttr<int> *>(attrs["position"]);
		if(!coord)
			throw "Position attr has been overloaded, Use DIFF normal strategy instead.";

		//estimate normals using vertices and faces existing.
		std::vector<Point3f> estimated;
		estimateNormals(coord->values, index, estimated);

		for(uint32_t i = 0; i < nvert; i++) {
			Point2i o = toOcta(estimated[i], q);
			diffs[i*2] = o[0];
			diffs[i*2+1] = o[1];
		}

		if(prediction == BORDER)
			markBoundary(values.size()/2, index); //mark boundary points on original vertices.
		for(uint32_t i = 0; i < values.size(); i++)
			values[i] -= diffs[i];
	}

	virtual void deltaEncode(std::vector<Quad> &context) {
		for(int c = 0; c < N; c++)
			diffs[c] = values[context[0].t*N + c];

		for(uint32_t i = 1; i < context.size(); i++) {
			Quad &quad = context[i];
			if(prediction == DIFF) {
				for(int c = 0; c < N; c++) {
					int &d = diffs[i*N + c];
					d = values[quad.t*N + c] - values[quad.a*N + c];
					if(d < -q)      d += 2*q;
					else if(d > q) d -= 2*q;
				}
			} else //just reorder diffs.
				for(int c = 0; c < N; c++)
					diffs[i*N + c] = values[quad.t*N + c];
		}
	}

	virtual void decode(uint32_t nvert, Stream &stream) {
		diffs.resize(nvert*2);
		int readed = stream.decodeArray<int32_t>(&*diffs.begin(), N);

		if(prediction == BORDER)
			diffs.resize(readed*2);
	}

	template <class T> void deltaDecodeT(uint32_t nvert, std::vector<Face> &context) {
		T *values = (T *)buffer;

		if(context.size()) {
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

	virtual void deltaDecode(uint32_t nvert, std::vector<Face> &context) {
		if(prediction != DIFF)
			return;

		switch(format) {
		case FLOAT:
		case INT32:
			deltaDecodeT<int32_t>(nvert, context);
			break;
		case INT16:
			deltaDecodeT<int16_t>(nvert, context);
		default: throw "Format not supported for normal attribute (float, int32 or int16 only)";
		}
	}

	virtual void postDelta(uint32_t nvert, std::map<std::string, Attribute23 *> &attrs, std::vector<uint32_t> &index) {
		//for border and estimate we need the position already deltadecoded but before dequantized
		if(prediction == DIFF)
			return;

		if(attrs.find("position") == attrs.end())
			throw "No position attribute found. Use DIFF normal strategy instead.";

		GenericAttr<int> *coord = dynamic_cast<GenericAttr<int> *>(attrs["position"]);
		if(!coord)
			throw "Position attr has been overloaded, Use DIFF normal strategy instead.";


		vector<Point3f> estimated(nvert, Point3f(0, 0, 0));
		estimateNormals(coord->values, index, estimated);

		if(prediction == BORDER)
			markBoundary(nvert, index);

		switch(format) {
		case FLOAT:
		case INT32:
			computeNormals((Point3f *)buffer, estimated);
			break;
		case INT16:
			computeNormals((Point3s *)buffer, estimated);
			break;
		default: throw "Format not supported for normal attribute (float, int32 or int16 only)";
		}
	}

	virtual void dequantize(uint32_t nvert) {
		if(prediction != DIFF)
			return;

		switch(format) {
		case FLOAT:
		case INT32:
			for(uint32_t i = 0; i < nvert; i++)
				*(Point3f *)buffer = toSphere(Point2i(diffs[i*2], diffs[i*2 + 1]), (int)q);
			break;
		case INT16:
			for(uint32_t i = 0; i < nvert; i++)
				*(Point3s *)buffer = toSphere(Point2s(diffs[i*2], diffs[i*2 + 1]), (int)q);
		default: throw "Format not supported for normal attribute (float, int32 or int16 only)";
		}
	}

	//Normal estimation

	void markBoundary(uint32_t nvert, std::vector<uint32_t> &index) {
		uint32_t nface = index.size()/3;
		boundary.resize(nvert, false);

		vector<int> count(nvert, 0);
		for(uint32_t i = 0; i < nface; i++) {
			uint32_t *f = &index[i*3];
			count[f[0]] += (int)f[1] - (int)f[2];
			count[f[1]] += (int)f[2] - (int)f[0];
			count[f[2]] += (int)f[0] - (int)f[1];
		}
		for(uint32_t i = 0; i < nvert; i++)
			if(count[i] != 0)
				boundary[i] = true;
	}


	void estimateNormals(std::vector<int32_t> &_coords, std::vector<uint32_t> &index, std::vector<Point3f> &estimated) {
		uint32_t nvert = _coords.size()/3;
		uint32_t nface = index.size()/3;
		Point3i *coords = (Point3i *)&*_coords.begin();

		estimated.resize(nvert);
		for(uint32_t i = 0; i < nface; i++) {
			uint32_t *f = &index[i*3];
			Point3i &p0 = coords[f[0]];
			Point3i &p1 = coords[f[1]];
			Point3i &p2 = coords[f[2]];
			Point3i qn = (( p1 - p0) ^ (p2 - p0));
			Point3f n(qn[0], qn[1], qn[2]);
			estimated[f[0]] += n;
			estimated[f[1]] += n;
			estimated[f[2]] += n;
		}
		for(uint32_t i = 0; i < nvert; i++) {
			Point3f &n = estimated[i];
			float len = n.norm();
			if(len < 0.00001) {
				continue;
			}
			n /= len;
		}
	}

	void computeNormals(Point3s *normals, std::vector<Point3f> &estimated) {
		uint32_t nvert = estimated.size();

		int count = 0; //here for the border.
		for(unsigned int i = 0; i < nvert; i++) {
			Point3f &e = estimated[i];
			int32_t *d = &diffs[count*2];
			Point3s &n = normals[i];

			if(prediction == ESTIMATED || boundary[i]) {
				Point2i qn = toOcta(e, (int)q);
				n = toSphere(Point2s(qn[0] + d[0], qn[0] + d[1]), (int)q);
				count++;
			} else {//no correction
				for(int k = 0; k < 3; k++)
					n[k] = (int16_t)(e[k]*32767.0f);
			}
		}
	}

	void computeNormals(Point3f *normals, std::vector<Point3f> &estimated) {
		uint32_t nvert = estimated.size();

		int count = 0; //here for the border.
		for(unsigned int i = 0; i < nvert; i++) {
			Point3f &e = estimated[i];
			int32_t *d = &diffs[count*2];
			Point3f &n = normals[i];
			if(prediction == ESTIMATED || boundary[i]) {
				Point2i qn = toOcta(e, (int)q);
				n = toSphere(Point2i(qn[0] + d[0], qn[0] + d[1]), (int)q);
				count++;
			} else //no correction
				n = e;
		}
	}



	//Conversion to Octahedron encoding.

	static Point2i toOcta(Point3f v, int unit) {
		Point2f p(v[0], v[1]);
		p /= (fabs(v[0]) + fabs(v[1]) + fabs(v[2]));

		if(v[2] < 0) {
			p = Point2f(1.0f - fabs(p[1]), 1.0f - fabs(p[0]));
			if(v[0] < 0) p[0] = -p[0];
			if(v[1] < 0) p[1] = -p[1];
		}
		return Point2i(p[0]*unit, p[1]*unit);
	}

	static Point2i toOcta(Point3i v, int unit) {
		float length = (float)v.norm();
		Point3f n(v[0]/length, v[1]/length, v[2]/length);
		return toOcta(n, unit);
	}

	static Point3f toSphere(Point2i v, int unit) {
		Point3f n(v[0], v[1], unit - abs(v[0]) -abs(v[1]));
		if (n[2] < 0) {
			n[0] = ((v[0] > 0)? 1 : -1)*(unit - abs(v[1]));
			n[1] = ((v[1] > 0)? 1 : -1)*(unit - abs(v[0]));
		}
		n /= n.norm();
		return n;
	}

	static Point3s toSphere(Point2s v, int unit) {
		Point3f n(v[0], v[1], unit - abs(v[0]) -abs(v[1]));
		if (n[2] < 0) {
			n[0] = ((v[0] > 0)? 1 : -1)*(unit - abs(v[1]));
			n[1] = ((v[1] > 0)? 1 : -1)*(unit - abs(v[0]));
		}
		n /= n.norm();
		return Point3s(n[0]*32767, n[1]*32767, n[2]*32767);
	}

};

} //namespace
#endif // NX_NORMAL_ATTRIBUTE_H
