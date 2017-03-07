#include "nxzencoder.h"
#include "nxzdecoder.h"


#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/import_ply.h>
#include<vcg/complex/algorithms/update/normal.h>

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

	nx::NxzEncoder encoder(coords.size(), index.size()/3);
	encoder.addCoords((float *)&*coords.begin(), &*index.begin());
	encoder.addNormals((float *)&*normals.begin(), 10, nx::ESTIMATED);
	encoder.addColors((unsigned char *)&*colors.begin());
	encoder.encode();

	cout << "Nvert: " << coords.size() << " Nface: " << index.size()/3 << endl;
	cout << "Compressed to: " << encoder.stream.size() << endl;
	cout << "Ratio: " << 100.0f*encoder.stream.size()/(coords.size()*12 + index.size()*12) << endl;
	cout << "Bpv: " << 8.0f*encoder.stream.size()/coords.size() << endl << endl;

	cout << "Header: " << encoder.header_size << " bpv: " << (float)encoder.header_size/coords.size() << endl;

	cout << "Coord bpv; " << 8.0f*encoder.coord.size/coords.size() << endl;
	cout << "Coord q: " << encoder.coord.q << endl << endl;

	cout << "Normal bpv; " << 8.0f*encoder.norm.size/coords.size() << endl;
	cout << "Normal q: " << encoder.norm.q << endl << endl;

	cout << "Color LCCA bpv; " << 8.0f*encoder.color[0].size/coords.size() << " "
			<< 8.0f*encoder.color[1].size/coords.size() << " "
			<< 8.0f*encoder.color[2].size/coords.size() << " "
			<< 8.0f*encoder.color[3].size/coords.size() << " " << endl;
	cout << "Color LCCA q: " << encoder.color[0].q << " "
		 << encoder.color[1].q << " "
		 << encoder.color[2].q << " "
		 << encoder.color[3].q << " " << endl << endl;


	cout << "Face bpv; " << 8.0f*encoder.face.size/coords.size() << endl;

	nx::NxzDecoder decoder(encoder.stream.size(), encoder.stream.buffer);
	return 0;
}

