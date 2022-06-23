#include "ply.h"
#include "base.h"
#include <glog/logging.h>

// generate polygon mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <boost/function_output_iterator.hpp>


typedef struct Vertex
{
	float x, y, z;
	float nx, ny, nz;
	uint8_t r, g, b;
} Vertex;

typedef struct Face
{
	unsigned char nverts; 		// number of vertex indices in list 
	int *verts; 		// vertex index list 
} Face;

char *elem_names[] = {(char*)"vertex", (char*)"face"};

PlyProperty vert_props[] = 
{
	// list of property information for a vertex
	{(char*)"x", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, x), 0, 0, 0, 0},
	{(char*)"y", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, y), 0, 0, 0, 0},
	{(char*)"z", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, z), 0, 0, 0, 0},
	{(char*)"nx", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, nx), 0, 0, 0, 0},
	{(char*)"ny", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, ny), 0, 0, 0, 0},
	{(char*)"nz", PLY_FLOAT, PLY_FLOAT, offsetof(Vertex, nz), 0, 0, 0, 0},
	{(char*)"red", PLY_UCHAR, PLY_UCHAR, offsetof(Vertex, r), 0, 0, 0, 0},
	{(char*)"green", PLY_UCHAR, PLY_UCHAR, offsetof(Vertex, g), 0, 0, 0, 0},
	{(char*)"blue", PLY_UCHAR, PLY_UCHAR, offsetof(Vertex, b), 0, 0, 0, 0},
};

PlyProperty face_props[] = 
{
	// list of property information for a vertex 
	{(char*)"vertex_indices", PLY_INT, PLY_INT, offsetof(Face, verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(Face, nverts)},
};

bool PLYPointLoad3(const string fileName, Pwn_vector& points)
{
	if (fileName.empty())
		return false;
	PlyFile *ply;
	PlyProperty **plist;
	int nelems;
	int nprops;
	int num_elems;
  	char **elist;
	int file_type;
	float version;
	char *elem_name;
	ply = ply_open_for_reading(fileName.c_str(), &nelems, &elist, &file_type, &version);

	if (ply == NULL || nelems <= 0){
		return 0;}

	for (int i = 0; i < nelems; i++)
	{
		elem_name = elist[i];
		plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);
		
		if (equal_strings((char*)"vertex", elem_name))
		{
			ply_get_property(ply, (char*)"vertex", &vert_props[0]);
			ply_get_property(ply, (char*)"vertex", &vert_props[1]);
			ply_get_property(ply, (char*)"vertex", &vert_props[2]);
			ply_get_property(ply, (char*)"vertex", &vert_props[3]);
			ply_get_property(ply, (char*)"vertex", &vert_props[4]);
			ply_get_property(ply, (char*)"vertex", &vert_props[5]);

			Vertex vertex;			
			for (unsigned long long j = 0; j < num_elems; j++)
			{
				ply_get_element(ply, (void *)&vertex);
                points.push_back(make_pair(Point_3(vertex.x, vertex.y, vertex.z),Vector_3(vertex.nx, vertex.ny, vertex.nz)));
			}
		}
	}
	ply_close(ply);
	return true;
}

bool PLYMeshLoad(string fileName, Mesh<Point_3>& mesh)
{
	if (fileName.empty())
		return false;
	mesh.vertices.clear();
	mesh.faces.clear();

	PlyFile *ply;
	PlyProperty **plist;
	int nelems;
	int nprops;
	int num_elems;
  	char **elist;
	int file_type;
	float version;
	char *elem_name;
	
	ply = ply_open_for_reading(fileName.c_str(), &nelems, &elist, &file_type, &version);

	if (ply == NULL || nelems <= 0)
		return 0;

	for (unsigned long long i = 0; i < nelems; i++)
	{
		elem_name = elist[i];
		plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);
		
		if (equal_strings("vertex", elem_name))
		{
			//mesh.vertices.resize(num_elems);
			ply_get_property(ply, "vertex", &vert_props[0]);
			ply_get_property(ply, "vertex", &vert_props[1]);
			ply_get_property(ply, "vertex", &vert_props[2]);

			Vertex vertex;			
			for (int j = 0; j < num_elems; j++)
			{
				ply_get_element(ply, (void *)&vertex);
				mesh.vertices.push_back(Point_3(vertex.x,vertex.y,vertex.z));
			}
		}

		if (equal_strings("face", elem_name))
		{
			//mesh.faces.resize(num_elems);
			ply_get_property(ply, "face", &face_props[0]);

			Face face;
			for (unsigned long long j = 0; j < num_elems; j++)
			{
				ply_get_element(ply, (void *)&face);
				mesh.faces.push_back(Point_3(face.verts[0],face.verts[1],face.verts[2]));
			}
		}
	}
	ply_close(ply);
	return true;
}

bool PLYPointSave(const string fileName, Pwn_vector& points, int type)
{
	if (points.empty() || fileName.empty())
		return false;

	unsigned long long pointsNum = points.size();
	PlyFile *ply;
	float version;
		
	ply = ply_open_for_writing(fileName.c_str(), 1, elem_names, type, &version);

	ply_element_count(ply, (char*)"vertex", pointsNum);
	ply_describe_property(ply, (char*)"vertex", &vert_props[0]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[1]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[2]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[3]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[4]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[5]);

	ply_header_complete(ply);
	
	Vertex vertex;
	ply_put_element_setup(ply, (char*)"vertex");
	for (unsigned long long i = 0; i < pointsNum; i++)
	{
		vertex.x = points[i].first.x();
		vertex.y = points[i].first.y();
		vertex.z = points[i].first.z();

		vertex.nx = points[i].second.x();
		vertex.ny = points[i].second.y();
		vertex.nz = points[i].second.z();

		ply_put_element(ply, (void *)&vertex);
	}	
	
	ply_close(ply);
	return true;
}

bool PLYMeshSave(const string fileName, Mesh<Point_3>& mesh, int type)
{
	if (mesh.vertices.empty() || fileName.empty())
		return false;

	unsigned long long vertexNum = mesh.vertices.size();
	unsigned long long faceNum = mesh.faces.size();

	PlyFile *ply;
	float version;
		
	ply = ply_open_for_writing(fileName.c_str(), 2, elem_names, type, &version);

	ply_element_count(ply, (char*)"vertex", vertexNum);
	ply_describe_property(ply, (char*)"vertex", &vert_props[0]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[1]);
	ply_describe_property(ply, (char*)"vertex", &vert_props[2]);

	ply_element_count(ply, (char*)"face", faceNum);
	ply_describe_property(ply, (char*)"face", &face_props[0]);

	ply_header_complete(ply);

	ply_put_element_setup(ply, (char*)"vertex");
	Vertex vertex;
	for (unsigned long long i = 0; i < vertexNum; i++)
	{
		vertex.x = mesh.vertices[i].x();
		vertex.y = mesh.vertices[i].y();
		vertex.z = mesh.vertices[i].z();
		ply_put_element(ply, (void *)&vertex);
	}
	
	ply_put_element_setup (ply, (char*)"face");
	Face face;
	for (unsigned long long i = 0; i < faceNum; i++)
	{
		face.nverts = 3;
		int v[3];
		v[0] = mesh.faces[i].x();
		v[1] = mesh.faces[i].y();
		v[2] = mesh.faces[i].z();
		face.verts = v;
		ply_put_element(ply, (void *) &face);
	}
	
	ply_close(ply);
	return true;
}

bool OFFMeshSave(const string fileName, Mesh<K_epec::Point_3>& mesh){
    // merge
    ofstream ofs(fileName);
    ofs << "COFF\n" << mesh.vertices.size() << ' ' << mesh.faces.size() << ' ' << "0\n";
     for (const auto &p : mesh.vertices)
       ofs << setiosflags(ios::fixed) << setprecision(8) << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
     for(auto poly: mesh.faces)
        ofs << 3 << " " << static_cast<int>(CGAL::to_double(poly[0])) << " " << static_cast<int>(CGAL::to_double(poly[1])) << " " << static_cast<int>(CGAL::to_double(poly[2])) << "\n";
     ofs.close();
}

// polygon mesh
typedef CGAL::Surface_mesh<Point_3> 		Mesh3;
typedef boost::graph_traits<Mesh3>::halfedge_descriptor 		halfedge_descriptor;
typedef boost::graph_traits<Mesh3>::edge_descriptor 		edge_descriptor;
typedef Mesh3::Vertex_index 		vertex_descriptor;
typedef Mesh3::Face_index 		face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

struct halfedge2edge
{
  halfedge2edge(const Mesh3& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Mesh3& m_mesh;
  std::vector<edge_descriptor>& m_edges;
};

bool PLYTriMeshSave(Mesh<Point_3>& ori_mesh, const cm::Config& config, string& name){
	Mesh3 tri_mesh;
	vector<vertex_descriptor> vertexmap(ori_mesh.vertices.size());
	for(int i = 0; i < ori_mesh.vertices.size(); i++){
		vertex_descriptor vi = tri_mesh.add_vertex(ori_mesh.vertices[i]);
		vertexmap[i] = vi; 
	}
    for(auto &f: ori_mesh.faces){
		tri_mesh.add_face(vertexmap[f.x()], vertexmap[f.y()], vertexmap[f.z()]);
	}
	LOG(INFO) << "Before remshing: #face: " << tri_mesh.number_of_faces() << ", #vertices: " << tri_mesh.number_of_vertices();	
	if (!CGAL::is_triangle_mesh(tri_mesh)) {
       LOG(INFO) << "Error: not a valid triangle mesh.";
       return false;
    }
	auto target_edge_length = config.get_optional<double>("mesh.target_edge_length");
	auto nb_iter = config.get_optional<unsigned int>("mesh.nb_iter");
	LOG(INFO) << "target_edge_length: " << target_edge_length;
	std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(tri_mesh),
		tri_mesh,
		boost::make_function_output_iterator(halfedge2edge(tri_mesh, border)));
	// recondition: split long edges
    PMP::split_long_edges(border, *target_edge_length, tri_mesh);
	// isotropic remeshing
	PMP::isotropic_remeshing(
		faces(tri_mesh),
		*target_edge_length,
		tri_mesh,
		PMP::parameters::number_of_iterations(*nb_iter)
		.protect_constraints(true)//i.e. protect border, here
	);
	LOG(INFO) << "After remeshing: #face: " << tri_mesh.number_of_faces() << " , #vertices: " << tri_mesh.number_of_vertices();

	Mesh<Point_3> mesh;
	for (auto &f : tri_mesh.faces()) {
		CGAL::Vertex_around_face_iterator<Mesh3> vbegin, vend;
		vector<int> idx;
		for (boost::tie(vbegin, vend) = CGAL::vertices_around_face(tri_mesh.halfedge(f), tri_mesh);
		  vbegin != vend;
		  ++vbegin) {
			auto it = find(mesh.vertices.begin(), mesh.vertices.end(), tri_mesh.point(*vbegin));
            if(it == mesh.vertices.end()){
				idx.push_back(mesh.vertices.size());
				mesh.vertices.push_back(tri_mesh.point(*vbegin));
			}
			else{
				idx.push_back(it - mesh.vertices.begin());
			}
		}
		if(idx.size() != 3){
			LOG(INFO) << "Facde num error.";
			return false;
		}
		mesh.faces.push_back({idx[0], idx[1], idx[2]});
	}
	PLYMeshSave(name, mesh, 1);
	return true;
}




