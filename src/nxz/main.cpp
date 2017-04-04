/*
Nexus

Copyright(C) 2012 - Federico Ponchio
ISTI - Italian National Research Council - Visual Computing Lab

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
for more details.
*/
#include <assert.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>

#include "nxzencoder.h"
#include "nxzdecoder.h"

#include "timer.h"
#include "tinyply.h"

using namespace tinyply;

using namespace nx;
using namespace std;

int main(int argc, char *argv[]) {
	if(argc != 2) {
		cerr << "Usage: " << argv[0] << " [file.ply]\n";
		return 0;
	}


	std::ifstream ss(argv[1], std::ios::binary);
	PlyFile ply(ss);

	std::vector<float> coords;
	std::vector<float> norms;
	std::vector<uint8_t> colors;
	std::vector<float> uvs;
	std::vector<float> radiuses;
	std::vector<uint32_t> index;
	std::vector<float> wedge_uvs;


	ply.request_properties_from_element("vertex", { "x", "y", "z" }, coords);
	ply.request_properties_from_element("vertex", { "nx", "ny", "nz" }, norms);
	ply.request_properties_from_element("vertex", { "red", "green", "blue", "alpha" }, colors);
	ply.request_properties_from_element("vertex", { "texcoord" }, uvs);
	ply.request_properties_from_element("vertex", { "radius" }, radiuses);
	ply.request_properties_from_element("face", { "vertex_indices" }, index, 3);
	ply.request_properties_from_element("face", { "texcoord" }, wedge_uvs, 6);
	ply.read(ss);

	uint32_t nface = index.size()/3;
	uint32_t nvert = coords.size()/3;

	if(wedge_uvs.size()) {
		std::vector<bool> visited(nvert, false);
		uvs.resize(nvert*2);
		std::multimap<uint32_t, uint32_t> duplicated;
		for(uint32_t i = 0; i < index.size(); i++) {
			uint32_t k = index[i];
			Point3f p = *(Point3f *)&coords[k*3];
			Point2f wuv = *(Point2f *)&wedge_uvs[i*2]; //wedge uv
			Point2f &uv = *(Point2f *)&uvs[k*2];
			if(!visited[k]) {
				visited[k] = true;
				uv = wuv;
				continue;
			}
			if(uv == wuv)
				continue;

			uint32_t found = 0xffffffff;
			auto d = duplicated.find(k);
			if(d != duplicated.end()) {
				auto range = duplicated.equal_range(k);
				for (auto i = range.first; i != range.second; ++i) {
					uint32_t j = i->second;
					if(*(Point2f *)&uvs[j*2] == wuv) {
						found = i->first;
						break;
					}
				}
			}
			if(found == 0xffffffff) {
				found = coords.size()/3;
				coords.push_back(p[0]);
				coords.push_back(p[1]);
				coords.push_back(p[2]);

				uvs.push_back(wuv[0]);
				uvs.push_back(wuv[1]);

				if(colors.size())
					for(int j = 0; j < 4; j++)
						colors.push_back(colors[k*4 + j]);

				//TODO do the same thing for all other attributes

				duplicated.insert(std::make_pair(k, found));
			}
			assert(found < coords.size()/3);
			index[i] = found;
		}
	}
	for(int i = 0; i < index.size(); i++)
		assert(index[i] < coords.size()/3);

	nvert = coords.size()/3;
	nface = index.size()/3;


	bool pointcloud = (nface == 0);
	//if(force_pointcloud)
	//	nface = 0;

	nx::NxzEncoder encoder(nvert, nface, nx::Stream::TUNSTALL);
	if(pointcloud)
		encoder.addPositions(&*coords.begin());
	else
		encoder.addPositions(&*coords.begin(), &*index.begin());
	//TODO add suppor for wedge and face attributes adding simplex attribute
	if(norms.size())
		encoder.addNormals(&*norms.begin(), 10, nx::NormalAttr::BORDER);

	if(colors.size())
		encoder.addColors(&*colors.begin());

	if(uvs.size())
		encoder.addUvs((float *)&*uvs.begin(), 1.0f/1024.0f);

	if(radiuses.size())
		encoder.addAttribute("radius", (char *)&*radiuses.begin(), nx::VertexAttribute::FLOAT, 1, 1.0f);
	encoder.encode();

	nvert = encoder.nvert;
	nface = encoder.nface;

	cout << "Nvert: " << nvert << " Nface: " << nface << endl;
	cout << "Compressed to: " << encoder.stream.size() << endl;
	cout << "Ratio: " << 100.0f*encoder.stream.size()/(nvert*12 + nface*12) << "%" << endl;
	cout << "Bpv: " << 8.0f*encoder.stream.size()/nvert << endl << endl;

	cout << "Header: " << encoder.header_size << " bpv: " << (float)encoder.header_size/nvert << endl;

	nx::VertexAttribute *coord = encoder.data["position"];
	cout << "Coord bpv; " << 8.0f*coord->size/nvert << endl;
	cout << "Coord q: " << coord->q << endl << endl;

	nx::VertexAttribute *norm = encoder.data["normal"];
	if(norm) {
		cout << "Normal bpv; " << 8.0f*norm->size/nvert << endl;
		cout << "Normal q: " << norm->q << endl << endl;
	}

	nx::ColorAttr *color = dynamic_cast<nx::ColorAttr *>(encoder.data["color"]);
	if(color) {
		cout << "Color bpv; " << 8.0f*color->size/nvert << endl;
		cout << "Color q: " << color->qc[0] << " " << color->qc[1] << " " << color->qc[2] << " " << color->qc[3] << endl << endl;
	}


	nx::GenericAttr<float> *radius = dynamic_cast<nx::GenericAttr<float> *>(encoder.data["radius"]);
	if(radius) {
		cout << "Radius  bpv; " << 8.0f*radius->size/nvert << endl;
		cout << "Radius q: " << radius->q << endl << endl;
	}


	cout << "Face bpv; " << 8.0f*encoder.index.size/nvert << endl;


	std::vector<float>    recoords(nvert*3);
	std::vector<float>    renorms(nvert*3);
	std::vector<uint8_t>  recolors(nvert*4);
	std::vector<float>    reuvs(nvert*2);
	std::vector<float>    reradius(nvert);
	std::vector<uint32_t> reindex(nface*3);
	std::vector<float> tex(nface*6);


	nx::Timer timer;

	int iter = 1;
	for(int i = 0; i < iter; i++) {
		nx::NxzDecoder decoder(encoder.stream.size(), encoder.stream.data());
		assert(decoder.nface == nface);
		assert(decoder.nvert = nvert);

		decoder.setPositions(&*recoords.begin());
		if(decoder.data.count("normal"))
			decoder.setNormals(&*renorms.begin());
		if(decoder.data.count("color"))
			decoder.setColors(&*recolors.begin());
		if(decoder.data.count("uv"))
			decoder.setUvs(&*reuvs.begin());
		if(decoder.data.count("radius"))
			decoder.setAttribute("radius", (char *)&*reradius.begin(), nx::VertexAttribute::FLOAT);

		decoder.setIndex(&*reindex.begin());
		decoder.decode();
	}

	int delta = timer.elapsed();
	if(nface) {
		float mfaces = nface*iter/1000000.0f;
		cout << "TOT M faces: " << mfaces << " in: " << delta << "ms or " << 1000*mfaces/delta << " MT/s" << endl;
	} else {
		float mverts = nvert * iter /1000000.0f;
		cout << "TOT M verts: " << mverts << " in: " << delta << "ms or " << 1000*mverts/delta << " MT/s" << endl;
	}

	FILE *file = fopen("test.nxz", "w");
	if(!file) {
		cerr << "Couldl not open file: " << "test.nxz" << endl;
		return 1;
	}
	size_t count = encoder.stream.size();
	size_t written = fwrite ( encoder.stream.data(), 1, count, file);
	if(written != count) {
		cerr << "Failed saving file: " << "test.nxz" << endl;
		return 1;
	}


	nx::NxzDecoder test_decoder(encoder.stream.size(), encoder.stream.data());

	std::filebuf fb;
	fb.open("test.ply", std::ios::out | std::ios::binary);
	std::ostream outputStream(&fb);
	PlyFile out;
	out.comments = ply.comments;

	out.add_properties_to_element("vertex", { "x", "y", "z" }, recoords);
	if(test_decoder.data.count("normal"))
		out.add_properties_to_element("vertex", { "nx", "ny", "nz" }, renorms);
	if(test_decoder.data.count("color"))
		out.add_properties_to_element("vertex", { "red", "green", "blue", "alpha" }, recolors);
	if(test_decoder.data.count("uv")) {
		for(int i = 0; i < reindex.size(); i++) {
			tex[i*2] = reuvs[reindex[i]*2];
			tex[i*2+1] = reuvs[reindex[i]*2+1];
		}
//		out.add_properties_to_element("texcoord", { "u", "v" }, reuvs);
		out.add_properties_to_element("face", { "texcoord" }, tex, 6, PlyProperty::Type::UINT8);
	}

	if(test_decoder.nface > 0)
		out.add_properties_to_element("face", { "vertex_indices" }, reindex, 3, PlyProperty::Type::UINT8);



	out.comments.push_back("generated by tinyply");
	out.write(outputStream, true);
	fb.close();

	return 0;
}

