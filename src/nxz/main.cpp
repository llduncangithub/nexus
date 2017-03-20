#include "nxzencoder.h"
#include "nxzdecoder.h"

#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/import_ply.h>
#include<vcg/complex/algorithms/update/normal.h>

#include<wrap/io_trimesh/export_ply.h>

#include <QTime>
#include <QFile>

class NxVertex; class NxEdge; class NxFace;
struct NxUsedTypes : public vcg::UsedTypes<vcg::Use<NxVertex>   ::AsVertexType,
											  vcg::Use<NxEdge>     ::AsEdgeType,
											  vcg::Use<NxFace>     ::AsFaceType>{};
class NxVertex  : public vcg::Vertex< NxUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::TexCoord2f, vcg::vertex::BitFlags  >{};
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
	NxMesh mesh;
	if(tri::io::ImporterPLY<NxMesh>::Open(mesh, argv[1]) != 0) {
		printf("Error reading file  %s\n",argv[1]);
		exit(0);
	}
	tri::UpdateNormal<NxMesh>::PerVertexNormalized(mesh);

	//create the buffers
	vector<vcg::Point3f> coords;
	vector<vcg::Point3f> normals;
	vector<vcg::Color4b> colors;
	vector<unsigned int> index;

	for(auto &v: mesh.vert) {
		coords.push_back(v.P());
		normals.push_back(v.N());
		colors.push_back(v.C());
	}

	NxMesh::VertexType *start = &*mesh.vert.begin();
	for(auto &f: mesh.face) {
		index.push_back(f.V(0) - start);
		index.push_back(f.V(1) - start);
		index.push_back(f.V(2) - start);
	}


	uint32_t nvert = mesh.vert.size();
	uint32_t nface = mesh.face.size();



//	nx::NxzEncoder encoder(nvert, 0, nx::Stream::TUNSTALL);
//	encoder.addCoords((float *)&*coords.begin());

	nx::NxzEncoder encoder(nvert, nface, nx::Stream::TUNSTALL);
	encoder.addCoords((float *)&*coords.begin(), &*index.begin());
//	encoder.addNormals((float *)&*normals.begin(), 10, nx::DIFF);
//	encoder.addColors((unsigned char *)&*colors.begin());
	encoder.encode();

	nvert = encoder.nvert;
	nface = encoder.nface;

	cout << "Nvert: " << nvert << " Nface: " << nface << endl;
	cout << "Compressed to: " << encoder.stream.size() << endl;
	cout << "Ratio: " << 100.0f*encoder.stream.size()/(nvert*12 + nface*12) << "%" << endl;
	cout << "Bpv: " << 8.0f*encoder.stream.size()/nvert << endl << endl;

	cout << "Header: " << encoder.header_size << " bpv: " << (float)encoder.header_size/nvert << endl;

	cout << "Coord bpv; " << 8.0f*encoder.coord.size/nvert << endl;
	cout << "Coord q: " << encoder.coord.q << endl << endl;

	cout << "Normal bpv; " << 8.0f*encoder.norm.size/nvert << endl;
	cout << "Normal q: " << encoder.norm.q << endl << endl;

	cout << "Color LCCA bpv; " << 8.0f*encoder.color[0].size/nvert << " "
			<< 8.0f*encoder.color[1].size/nvert << " "
			<< 8.0f*encoder.color[2].size/nvert << " "
			<< 8.0f*encoder.color[3].size/nvert << " " << endl;
	cout << "Color LCCA q: " << encoder.color[0].q << " "
		 << encoder.color[1].q << " "
		 << encoder.color[2].q << " "
		 << encoder.color[3].q << " " << endl << endl;


	cout << "Face bpv; " << 8.0f*encoder.face.size/nvert << endl;

	std::vector<vcg::Point3f> recoords(nvert);
	std::vector<vcg::Point3f> renorms(nvert, Point3f(0, 0, 0));
	std::vector<vcg::Color4b> recolors(nvert, Color4b(255, 255, 255, 255));
	std::vector<uint32_t> reindex(nface*3);

	QTime time;
	time.start();

	int iter = 20;
	for(int i = 0; i < iter; i++) {
		nx::NxzDecoder decoder(encoder.stream.size(), encoder.stream.data());
		decoder.setCoords((float *)&*recoords.begin());
//		decoder.setNormals((float *)&*renorms.begin());
//		decoder.setColors((uchar *)&*recolors.begin());
		decoder.setIndex(&*reindex.begin());
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
	}

	cout << "Nvert: " << nvert << " nface: " << nface << endl;
	for(uint32_t i = 0; i < nface; i++) {
		assert(reindex[i*3] < nvert);
		assert(reindex[i*3+1] < nvert);
		assert(reindex[i*3+2] < nvert);

		out.face[i].V(0) = &*vi + reindex[i*3+0];
		out.face[i].V(1) = &*vi + reindex[i*3+1];
		out.face[i].V(2) = &*vi + reindex[i*3+2];
	}

	vcg::tri::io::ExporterPLY<NxMesh>::Save(out, "test.ply",
		tri::io::Mask::IOM_VERTCOLOR | tri::io::Mask::IOM_VERTNORMAL );


	return 0;
}

