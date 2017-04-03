#ifndef NX_NORMAL_ATTRIBUTE_H
#define NX_NORMAL_ATTRIBUTE_H

#include "attribute.h"
#include "point.h"

namespace nx {

class NormalAttr: public Attribute23 {
public:
	enum Prediction { DIFF = 0x0,      //do not estimate normals, use diffs to previous
					  ESTIMATED = 0x1, //estimate normals then encode differences
					  BORDER = 0x2 };  //encode differences only on the boundary
	uint32_t prediction;
	std::vector<int32_t> boundary;
	std::vector<int32_t> values, diffs;

	NormalAttr(int bits = 10) {
		N = 3;
		q = pow(2, bits-1);
		prediction = DIFF;
		strategy |= Attribute23::CORRELATED;
	}

	virtual void quantize(uint32_t nvert, char *buffer);
	virtual void preDelta(uint32_t nvert,  uint32_t nface, std::map<std::string, Attribute23 *> &attrs, IndexAttr &index);
	virtual void deltaEncode(std::vector<Quad> &context);
	virtual void encode(uint32_t nvert, Stream &stream);

	virtual void decode(uint32_t nvert, Stream &stream);
	virtual void deltaDecode(uint32_t nvert, std::vector<Face> &context);
	virtual void postDelta(uint32_t nvert,  uint32_t nface, std::map<std::string, Attribute23 *> &attrs, IndexAttr &index);
	virtual void dequantize(uint32_t nvert);

	//Normal estimation
	void computeNormals(Point3s *normals, std::vector<Point3i> &estimated);
	void computeNormals(Point3f *normals, std::vector<Point3i> &estimated);


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

		Point2i p(v[0]*unit, v[1]*unit);
		p /= (fabs(v[0]) + fabs(v[1]) + fabs(v[2]));

		if(v[2] < 0) {
			p = Point2i(unit - fabs(p[1]), unit - fabs(p[0]));
			if(v[0] < 0) p[0] = -p[0];
			if(v[1] < 0) p[1] = -p[1];
		}
		return p;
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
