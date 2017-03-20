#ifndef NX_POINT_H
#define NX_POINT_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef unsigned char uchar;

namespace nx {

template <typename S> class Point2 {
private:
	S v[2];

public:
	Point2() {}
	Point2(S x, S y) { v[0] = x; v[1] = y; }
	Point2(S *x) { v[0] = x[0]; v[1] = x[1]; }
	explicit Point2(S x) { v[0] = v[1] = x; }

	S &operator[](int k) { return v[k]; }
	const S &operator[](int k) const { return v[k]; }

	void operator=(const Point2<float> &p) { v[0] = (S)p[0]; v[1] = (S)p[1]; }


	bool operator==(const Point2 &c) const { return v[0] == c(0) && v[1] == c(1); }
	bool operator!=(const Point2 &c) const { return v[0] != c(0) || v[1] != c(1); }
	bool operator<(const Point2 &c) const {
		if(v[0] == c[0]) return v[1] < c[1];
		return v[0] < c[0];
	}
	void setMin(const Point2 &c) {
		if(c[0] < v[0]) v[0] = c[0];
		if(c[1] < v[1]) v[1] = c[1];
	}
	void setMax(const Point2 &c) {
		if(c[0] > v[0]) v[0] = c[0];
		if(c[1] > v[1]) v[1] = c[1];
	}

	S operator*(const Point2 &c) { return v[0]*c[0] + v[1]*c[1]; }
	S norm() { return (S)sqrt((double)(v[0]*v[0] + v[1]*v[1])); }

	Point2 operator-() const     { return Point2(-v[0], -v[1]); }
	Point2 operator+(const Point2 &c) const { return Point2(v[0] + c[0], v[1] + c[1]); }
	Point2 operator-(const Point2 &c) const { return Point2(v[0] - c[0], v[1] - c[1]); }
	Point2 operator*(S c) const { return Point2(v[0]*c, v[1]*c); }
	Point2 operator/(S c) const { return Point2(v[0]/c, v[1]/c); }
	Point2 operator>>(int i) const { return Point2(v[0]>>i, v[1]>>i); }
	Point2 operator<<(int i) const { return Point2(v[0]<<i, v[1]<<i); }

	Point2 &operator+=(const Point2 &c) {
		v[0] += c[0];
		v[1] += c[1];
		return *this;
	}

	Point2 &operator-=(const Point2 &c) {
		v[0] -= c[0];
		v[1] -= c[1];
		return *this;
	}

	Point2 &operator*=(S c)  {
		v[0] *= c;
		v[1] *= c;
		return *this;
	}
	Point2 &operator/=(S c)  {
		v[0] /= c;
		v[1] /= c;
		return *this;
	}
};

typedef Point2<int16_t> Point2s;
typedef Point2<int32_t> Point2i;
typedef Point2<float>   Point2f;
typedef Point2<double>  Point2d;


template <typename S> class Point3 {
private:
	S v[3];
public:
	Point3() {}
	Point3(S x, S y, S z) { v[0] = x; v[1] = y; v[2] = z; }
	Point3(S *x) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; }
	explicit Point3(S x) { v[0] = v[1] = v[2] = x; }
	//explicit Point3(Point3<float> &p) { v[0] = (S)p[0]; v[1] = (S)p[1]; v[2] = (S)p[2]; }
	void operator=(const Point3<float> &p) { v[0] = (S)p[0]; v[1] = (S)p[1]; v[2] = (S)p[2]; }

	S &operator[](int k) { return v[k]; }
	const S &operator[](int k) const { return v[k]; }

	S norm() { return (S)sqrt((double)(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])); }
	S operator*(const Point3 &c) { return v[0]*c[0] + v[1]*c[1] + v[2]*c[2]; }
	Point3 operator^(const Point3 &p) const {
		return Point3(v[1]*p.v[2] - v[2]*p.v[1], v[2]*p.v[0] - v[0]*p.v[2], v[0]*p.v[1] - v[1]*p.v[0]);
	}

	bool operator==(const Point3 &c) const { return v[0] == c[0] && v[1] == c[1] && v[2] == c[2]; }
	bool operator!=(const Point3 &c) const { return !(c == *this); }
	bool operator<(const Point3 &c) const {
		if(v[0] == c[0]) {
			if(v[1] == c[1])
				return v[2] < c[2];
			return v[1] < c[1];
		}
		return v[0] < c[0];
	}
	void setMin(const Point3 &c) {
		if(c[0] < v[0]) v[0] = c[0];
		if(c[1] < v[1]) v[1] = c[1];
		if(c[2] < v[2]) v[2] = c[2];
	}
	void setMax(const Point3 &c) {
		if(c[0] > v[0]) v[0] = c[0];
		if(c[1] > v[1]) v[1] = c[1];
		if(c[2] > v[2]) v[2] = c[2];
	}

	Point3 operator-() const { return Point3(-v[0], -v[1], -v[2]); }

	Point3 operator+(const Point3 &c) const { return Point3(v[0] + c[0], v[1] + c[1], v[2] + c[2]); }
	Point3 operator-(const Point3 &c) const { return Point3(v[0] - c[0], v[1] - c[1], v[2] - c[2]); }
	Point3 operator*(S c) const { return Point3(v[0]*c, v[1]*c, v[2]*c); }
	Point3 operator/(S c) const { return Point3(v[0]/c, v[1]/c, v[2]/c); }

	Point3 &operator-=(const Point3 &c) {
		v[0] -= c[0];
		v[1] -= c[1];
		v[2] -= c[2];
		return *this;
	}
	Point3 &operator+=(const Point3 &c) {
		v[0] += c[0];
		v[1] += c[1];
		v[2] += c[2];
		return *this;
	}
	Point3 &operator*=(S c)  {
		v[0] *= c;
		v[1] *= c;
		v[2] *= c;
		return *this;
	}
	Point3 &operator/=(S c)  {
		v[0] /= c;
		v[1] /= c;
		v[2] /= c;
		return *this;
	}
};

typedef Point3<int16_t> Point3s;
typedef Point3<int32_t> Point3i;
typedef Point3<float>   Point3f;
typedef Point3<double>  Point3d;

template <typename S> class Point4 {
protected:
	S v[4];
public:
	Point4() {}
	Point4(S x, S y, S z, S w) { v[0] = x; v[1] = y; v[2] = z; v[3] = w; }
	Point4(S *x) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; v[3] = x[3]; }
	explicit Point4(S x) { v[0] = v[1] = v[2] = v[3] = x; }

	S &operator[](int k) { return v[k]; }
	const S &operator[](int k) const { return v[k]; }
	Point4 &operator-=(const Point4 &c) {
		v[0] -= c[0];
		v[1] -= c[1];
		v[2] -= c[2];
		v[3] -= c[3];
		return *this;
	}
	Point4 &operator+=(const Point4 &c) {
		v[0] += c[0];
		v[1] += c[1];
		v[2] += c[2];
		v[3] += c[3];
		return *this;
	}
};

class Color4b: public Point4<unsigned char> {
public:
	Color4b() {}
	Color4b(uchar r, uchar g, uchar b, uchar a) { v[0] = r; v[1] = g; v[2] = b; v[3] = a; }
	Color4b toYCC() { return Color4b(v[1], (uchar)(v[2] - v[1]), (uchar)(v[0] - v[1]), v[3]); }
	Color4b toRGB() { return Color4b(v[2] + v[0], v[0], v[1] + v[0], v[3]); }
};

/* convenience methods */

/*Point2i lroundf(Point2f &p) { return Point2i(::lroundf(p[0]), ::lroundf(p[1])); }
Point3i lroundf(Point3f &p) { return Point3i(::lroundf(p[0]), ::lroundf(p[1]), ::lroundf(p[2])); } */

/* mapping sphere to octahedron */
class Normal {
public:
	static Point2i encode(Point3f v, int unit) {
		Point2f p(v[0], v[1]);
		p /= (fabs(v[0]) + fabs(v[1]) + fabs(v[2]));

		if(v[2] < 0) {
			p = Point2f(1.0f - fabs(p[1]), 1.0f - fabs(p[0]));
			if(v[0] < 0) p[0] = -p[0];
			if(v[1] < 0) p[1] = -p[1];
		}
		return Point2i(p[0]*unit, p[1]*unit);
	}

	static Point3f decode(Point2i v, int unit) {
		Point3f n(v[0], v[1], unit - abs(v[0]) -abs(v[1]));
		if (n[2] < 0) {
			n[0] = ((v[0] > 0)? 1 : -1)*(unit - abs(v[1]));
			n[1] = ((v[1] > 0)? 1 : -1)*(unit - abs(v[0]));
		}
		n /= n.norm();
		return n;
	}

	static Point3s decode(Point2s v, int unit) {
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
#endif // NX_POINT_H
