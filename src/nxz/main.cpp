#include "nxzencoder.h"

#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/import_ply.h>
#include<vcg/complex/algorithms/update/normal.h>

class CoordVertex; class CoordEdge; class CoordFace;
struct CoordUsedTypes : public vcg::UsedTypes<vcg::Use<CoordVertex>   ::AsVertexType,
											  vcg::Use<CoordEdge>     ::AsEdgeType,
											  vcg::Use<CoordFace>     ::AsFaceType>{};
class CoordVertex  : public vcg::Vertex< CoordUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >{};
class CoordFace    : public vcg::Face<   CoordUsedTypes, vcg::face::VertexRef, vcg::face::BitFlags > {};
class CoordEdge    : public vcg::Edge<CoordUsedTypes> {};
class CoordMesh    : public vcg::tri::TriMesh< std::vector<CoordVertex>, std::vector<CoordFace> , std::vector<CoordEdge>  > {};

/*
class ColorVertex; class ColorEdge; class ColorFace;
struct ColorUsedTypes : public vcg::UsedTypes<vcg::Use<ColorVertex>   ::AsVertexType,
										   vcg::Use<ColorEdge>     ::AsEdgeType,
										   vcg::Use<ColorFace>     ::AsFaceType>{};
class ColorEdge    : public vcg::Edge<ColorUsedTypes> {};

class ColorVertex  : public vcg::Vertex< ColorUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags  >{};
class ColorFace    : public vcg::Face<   ColorUsedTypes, vcg::face::VertexRef, vcg::face::BitFlags > {};
class ColorMesh    : public vcg::tri::TriMesh< std::vector<ColorVertex>, std::vector<CoordFace> , std::vector<ColorEdge>  > {};


class UvVertex; class UvEdge; class UvFace;
struct UvUsedTypes : public vcg::UsedTypes<vcg::Use<UvVertex>   ::AsVertexType,
										   vcg::Use<UvEdge>     ::AsEdgeType,
										   vcg::Use<UvFace>     ::AsFaceType>{};
class UvEdge    : public vcg::Edge<UvUsedTypes> {};

class UvVertex  : public vcg::Vertex< UvUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::TexCoord2f, vcg::vertex::BitFlags  >{};
class UvFace    : public vcg::Face<   UvUsedTypes, vcg::face::VertexRef, vcg::face::BitFlags > {};
class UvMesh    : public vcg::tri::TriMesh< std::vector<UvVertex>, std::vector<UvFace> , std::vector<UvEdge>  > {}; */

using namespace std;
using namespace vcg;

int main(int argc, char *argv[]) {

	CoordMesh mesh;
	if(tri::io::ImporterPLY<CoordMesh>::Open(mesh, argv[1]) != 0) {
		printf("Error reading file  %s\n",argv[1]);
		exit(0);
	}
	tri::UpdateNormal<CoordMesh>::PerVertexNormalized(mesh);

	//create the buffers
	vector<Point3f> coords;
	vector<Point3f> normals;
	vector<unsigned int> index;

	for(auto &v: mesh.vert)
		coords.push_back(v.P());

	return 0;
}

