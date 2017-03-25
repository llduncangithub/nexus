#include "normal_attribute.h"

using namespace nx;



void NormalAttr::quantize(uint32_t nvert, char *buffer) {
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


void NormalAttr::preDelta(uint32_t nvert, std::map<std::string, Attribute23 *> &attrs, std::vector<uint32_t> &index) {
	if(prediction == DIFF)
		return;

	if(attrs.find("position") == attrs.end())
		throw "No position attribute found. Use DIFF normal strategy instead.";

	GenericAttr<int> *coord = dynamic_cast<GenericAttr<int> *>(attrs["position"]);
	if(!coord)
		throw "Position attr has been overloaded, Use DIFF normal strategy instead.";

	//estimate normals using vertices and faces existing.
	std::vector<Point3f> estimated;
	estimateNormals(nvert, (Point3i *)&*coord->values.begin(), index, estimated);

	Point2i *d = (Point2i *)&*diffs.begin();
	for(uint32_t i = 0; i < nvert; i++)
		d[i] = toOcta(estimated[i], q);

	if(prediction == BORDER)
		markBoundary(values.size()/2, index); //mark boundary points on original vertices.

	for(uint32_t i = 0; i < values.size(); i++)
		values[i] -= diffs[i];
}

void NormalAttr::deltaEncode(std::vector<Quad> &context) {
	diffs[0] = values[context[0].t*2];
	diffs[1] = values[context[0].t*2+1];

	if(prediction == DIFF) {

		for(uint32_t i = 1; i < context.size(); i++) {
			Quad &quad = context[i];
			for(int c = 0; c < N; c++) {
				int &d = diffs[i*N + c];
				d = values[quad.t*N + c] - values[quad.a*N + c];
				if(d < -q)     d += 2*q;
				else if(d > q) d -= 2*q;
			}
		}
		diffs.resize(context.size()*N); //unreferenced vertices

	} else  {//just reorder diffs, for border story only boundary diffs
		uint32_t count = 1;

		for(uint32_t i = 1; i < context.size(); i++) {
			Quad &quad = context[i];
			if(prediction != BORDER || boundary[i]) {
				for(int c = 0; c < N; c++)
					diffs[count*N + c] = values[quad.t*N + c];
				count++;
			}
		}
		diffs.resize(count*N); //unreferenced vertices and borders
	}
}

void NormalAttr::decode(uint32_t nvert, Stream &stream) {
	prediction = stream.read<uchar>();
	diffs.resize(nvert*2);
	int readed = stream.decodeArray<int32_t>(&*diffs.begin(), N);

	if(prediction == BORDER)
		diffs.resize(readed*2);
}

void NormalAttr::deltaDecode(uint32_t nvert, std::vector<Face> &context) {
	if(prediction != DIFF)
		return;

	if(context.size()) {
		for(uint32_t i = 1; i < context.size(); i++) {
			Face &f = context[i];
			for(int c = 0; c < 2; c++) {
				int &d = diffs[i*2 + c];
				d += diffs[f.a*2 + c];
				if(d < -q)     d += 2*q;
				else if(d > q) d -= 2*q;
			}

		}
	} else { //point clouds assuming values are already sorted by proximity.
		for(uint32_t i = 2; i < nvert*2; i++) {
			int &d = diffs[i];
			d += diffs[i-2];
			if(d < -q)     d += 2*q;
			else if(d > q) d -= 2*q;
		}
	}
}

void NormalAttr::postDelta(uint32_t nvert,
						   std::map<std::string, Attribute23 *> &attrs,
						   std::vector<uint32_t> &index) {
	//for border and estimate we need the position already deltadecoded but before dequantized
	if(prediction == DIFF)
		return;

	if(attrs.find("position") == attrs.end())
		throw "No position attribute found. Use DIFF normal strategy instead.";

	GenericAttr<int> *coord = dynamic_cast<GenericAttr<int> *>(attrs["position"]);
	if(!coord)
		throw "Position attr has been overloaded, Use DIFF normal strategy instead.";


	vector<Point3f> estimated(nvert, Point3f(0, 0, 0));
	estimateNormals(nvert, (Point3i *)coord->buffer, index, estimated);

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

void NormalAttr::dequantize(uint32_t nvert) {
	if(prediction != DIFF)
		return;

	switch(format) {
	case FLOAT:
	case INT32:
		for(uint32_t i = 0; i < nvert; i++)
			((Point3f *)buffer)[i] = toSphere(Point2i(diffs[i*2], diffs[i*2 + 1]), (int)q);
		break;
	case INT16:
		for(uint32_t i = 0; i < nvert; i++)
			((Point3s *)buffer)[i] = toSphere(Point2s(diffs[i*2], diffs[i*2 + 1]), (int)q);
	default: throw "Format not supported for normal attribute (float, int32 or int16 only)";
	}
}

//Normal estimation

void NormalAttr::markBoundary(uint32_t nvert, std::vector<uint32_t> &index) {
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


void NormalAttr::estimateNormals(uint32_t nvert, Point3i *coords, std::vector<uint32_t> &index, std::vector<Point3f> &estimated) {
	uint32_t nface = index.size()/3;

	estimated.resize(nvert, Point3f(0.0f, 0.0f, 0.0f));
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
		if(len < 0.00001f) {
			n = Point3f(0.0f, 0.0f, 1.0f);
			continue;
		}
		n /= len;
	}
}

void NormalAttr::computeNormals(Point3s *normals, std::vector<Point3f> &estimated) {
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

void NormalAttr::computeNormals(Point3f *normals, std::vector<Point3f> &estimated) {
	uint32_t nvert = estimated.size();

	/*	Point3f tn(0.7, 0.7, 0);
	Point2i tq = toOcta(tn, (int)q);
	tn = toSphere(tq, (int)q);
	cout << "TQ: " <<	tq[0] << " " << tq[1] << " N: " << tn[0] << " " << tn[1] << " " << tn[2] << endl;*/

	Point2i *diffp = (Point2i *)&*diffs.begin();
	int count = 0; //here for the border.
	for(unsigned int i = 0; i < nvert; i++) {
		Point3f &e = estimated[i];
		Point2i &d = diffp[count];
		Point3f &n = normals[i];
		if(prediction == ESTIMATED || boundary[i]) {
			Point2i qn = toOcta(e, (int)q);
			n = toSphere(qn + d, (int)q);
			count++;
		} else //no correction
			n = e;
	}
}
