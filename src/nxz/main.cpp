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
#include "nxzencoder.h"
#include "nxzdecoder.h"

#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/import_ply.h>
#include<vcg/complex/algorithms/update/normal.h>

#include<wrap/io_trimesh/export_ply.h>

#include <QTime>
#include <QFile>

#include "cstream.h"

class NxVertex; class NxEdge; class NxFace;
struct NxUsedTypes : public vcg::UsedTypes<vcg::Use<NxVertex>   ::AsVertexType,
											  vcg::Use<NxEdge>     ::AsEdgeType,
											  vcg::Use<NxFace>     ::AsFaceType>{};
class NxVertex  : public vcg::Vertex< NxUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b,
		vcg::vertex::TexCoord2f, vcg::vertex::Radiusf, vcg::vertex::BitFlags  >{};
class NxFace    : public vcg::Face<   NxUsedTypes, vcg::face::VertexRef, vcg::face::BitFlags > {};
class NxEdge    : public vcg::Edge<NxUsedTypes> {};
class NxMesh    : public vcg::tri::TriMesh< std::vector<NxVertex>, std::vector<NxFace> , std::vector<NxEdge>  > {};

using namespace std;
using namespace vcg;

int main(int argc, char *argv[]) {
	if(argc != 2) {
		cerr << "Usage: " << argv[0] << " [file.ply]\n";
		return 0;
	}
/*
	for(int skip = 2; skip < 200; skip += 1) {
		nx::Stream stream;
		vector<uchar> a;
		a.resize(1000, 0);
		for(int i = 0; i < a.size(); i+= skip)
			a[i] = 1;

		stream.compress(a.size(), &*a.begin());

		stream.rewind();
		vector<uchar> b;
		stream.decompress(b);
	} */

	NxMesh mesh;
	int loadmask = 0;
//		vcg::tri::io::Mask::IOM_VERTCOORD, IOM_VERTCOLOR, IOM_VERTQUALITY, IOM_VERTNORMAL, IOM_VERTTEXCOORD, IOM_VERTRADIUS,
//				IOM_FACECOLOR, IOM_FACEQUALITY, IOM_WEDGTEXCOORD

	if(tri::io::ImporterPLY<NxMesh>::Open(mesh, argv[1], loadmask) != 0) {
		printf("Error reading file  %s\n",argv[1]);
		exit(0);
	}
	tri::UpdateNormal<NxMesh>::PerVertexNormalized(mesh);

	//create the buffers
	vector<vcg::Point3f> coords;
	vector<vcg::Point3f> normals;
	vector<vcg::Color4b> colors;
	vector<vcg::TexCoord2f> uvs;
	vector<float> radiuses;
	vector<unsigned int> index;

	for(auto &v: mesh.vert) {
		coords.push_back(v.P());
		normals.push_back(v.N());
		colors.push_back(v.C());
		uvs.push_back(v.T());
		radiuses.push_back(v.R());
	}

	NxMesh::VertexType *start = &*mesh.vert.begin();
	for(auto &f: mesh.face) {
		index.push_back(f.V(0) - start);
		index.push_back(f.V(1) - start);
		index.push_back(f.V(2) - start);
	}


	uint32_t nvert = mesh.vert.size();
	uint32_t nface = mesh.face.size();

	//TODO checks for: strategy cannot be parallel in point clouds (ALL attributes)
	//                 normals strategy must be DIFF

	bool pointcloud = false;
	if(pointcloud)
		nface = 0;

	nx::NxzEncoder encoder(nvert, nface, nx::Stream::TUNSTALL);
	if(pointcloud)
		encoder.addPositions((float *)&*coords.begin());
	else
		encoder.addPositions((float *)&*coords.begin(), &*index.begin());
	//		vcg::tri::io::Mask::IOM_VERTCOORD, IOM_VERTCOLOR, IOM_VERTQUALITY, IOM_VERTNORMAL, IOM_VERTTEXCOORD, IOM_VERTRADIUS,
	//				IOM_FACECOLOR, IOM_FACEQUALITY, IOM_WEDGTEXCOORD
	if(loadmask & vcg::tri::io::Mask::IOM_VERTNORMAL)
		encoder.addNormals((float *)&*normals.begin(), 10, nx::NormalAttr::BORDER);
	if(loadmask & vcg::tri::io::Mask::IOM_VERTCOLOR)
		encoder.addColors((unsigned char *)&*colors.begin());
	if(loadmask & vcg::tri::io::Mask::IOM_VERTTEXCOORD)
		encoder.addUvs((float *)&*uvs.begin(), 2.0f/1024.0f);
	if(loadmask & vcg::tri::io::Mask::IOM_VERTRADIUS)
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

	std::vector<vcg::Point3f> recoords(nvert);
	std::vector<vcg::Point3f> renorms(nvert, vcg::Point3f(0.0f, 0.0f, 0.0f));
	std::vector<vcg::Color4b> recolors(nvert, vcg::Color4b(255, 255, 255, 255));
	std::vector<vcg::TexCoord2f> reuvs(nvert, vcg::TexCoord2f(0.0f, 0.0f));
	std::vector<uint32_t> reindex(nface*3);


	QTime time;
	time.start();

	int iter = 10;
	for(int i = 0; i < iter; i++) {
		nx::NxzDecoder decoder(encoder.stream.size(), encoder.stream.data());
		decoder.setPositions((float *)&*recoords.begin());
		if(decoder.data.count("normal"))
			decoder.setNormals((float *)&*renorms.begin());
		if(decoder.data.count("color"))
			decoder.setColors((uchar *)&*recolors.begin());
		if(decoder.data.count("uv"))
			decoder.setUvs((float *)&*reuvs.begin());


#ifndef POINTCLOUD
		decoder.setIndex(&*reindex.begin());
#endif
		decoder.decode();
	}

	int delta = time.elapsed();
	if(nface) {
		float mfaces = nface*iter/1000000.0f;
		cout << "TOT M faces: " << mfaces << " in: " << delta << "ms or " << 1000*mfaces/delta << " MT/s" << endl;
	} else {
		float mverts = nvert * iter /1000000.0f;
		cout << "TOT M verts: " << mverts << " in: " << delta << "ms or " << 1000*mverts/delta << " MT/s" << endl;
	}

	QFile file("test.nxz");
	file.open(QFile::WriteOnly);
	file.write((char *)encoder.stream.data(), encoder.stream.size());


	NxMesh out;
	NxMesh::VertexIterator vi = vcg::tri::Allocator<NxMesh>::AddVertices(out, nvert);
	vcg::tri::Allocator<NxMesh>::AddFaces(out, nface);
	for(uint32_t i = 0; i < nvert; i++) {
		out.vert[i].P() = recoords[i];
		out.vert[i].N() = renorms[i];
		out.vert[i].C() = recolors[i];
		out.vert[i].T() = reuvs[i];
	}

	cout << "Nvert: " << nvert << " nface: " << nface << endl;
	for(uint32_t i = 0; i < nface; i++) {
		assert(reindex[i*3] != reindex[i*3+1]);
		assert(reindex[i*3+1] != reindex[i*3+2]);
		assert(reindex[i*3+2] != reindex[i*3]);

		assert(reindex[i*3] < nvert);
		assert(reindex[i*3+1] < nvert);
		assert(reindex[i*3+2] < nvert);

		out.face[i].V(0) = &*vi + reindex[i*3+0];
		out.face[i].V(1) = &*vi + reindex[i*3+1];
		out.face[i].V(2) = &*vi + reindex[i*3+2];
	}

	vcg::tri::io::ExporterPLY<NxMesh>::Save(out, "test.ply",
		tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTNORMAL | tri::io::Mask::IOM_VERTTEXCOORD );


	return 0;
}

