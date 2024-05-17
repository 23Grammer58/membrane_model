//
// Created by alex on 13.06.2020.
//

#ifndef AORTIC_VALVE_TRIANGULARMESHHELPERS_H
#define AORTIC_VALVE_TRIANGULARMESHHELPERS_H
#include "AVMesh.h"
#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <tuple>
#include <cfloat>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

using namespace std;
using namespace World3d;

struct AniMesh{
    std::vector<double> vertices;
    std::vector<int> faces;
    std::vector<int> edges;
    std::vector<int> vertex_label;
    std::vector<int> face_label;
    std::vector<int> edge_label;
public:
    AniMesh(){}
    AniMesh(std::vector<double> vertices, std::vector<int> faces, std::vector<int> edges, std::vector<int> vertex_label,
            std::vector<int> face_label, std::vector<int> edge_label):
            vertices{vertices},
            faces(faces),
            edges(edges),
            vertex_label(vertex_label),
            face_label(face_label),
            edge_label(edge_label)
            {}
    void clear(){
        vertices.clear();
        faces.clear();
        edges.clear();
        vertex_label.clear();
        face_label.clear();
        edge_label.clear();
    }
};

Mesh makeMesh(std::vector<array<double, 3>>& points, std::vector<std::vector<std::size_t>>& polygons);
void set_edges_animesh(AniMesh& am);
Mesh split_twice(const Mesh& mesh);
AniMesh convert_to_animesh(const Mesh& mesh, std::string v_lbl = "", std::string f_lbl = "", std::string e_lbl = "");
void invertFaceOrientation(AniMesh& am);
Mesh convert_to_Mesh(const AniMesh& am, std::string v_lbl = "", std::string f_lbl = "", std::string e_lbl = "");
[[maybe_unused]] Mesh::Property_map<V_ind, int> remark_boundary(Mesh& m, int bnd_mrk = 1, std::string v_lbl = "v:boundary_lbl");
void translate_Mesh(Mesh& mesh, Vector translate);
bool write_STL(std::vector<CGAL::Simple_cartesian<double>::Triangle_3> triags, std::string filename);

//attention: this function doesn't copy property map, excluding "v:boundary_lbl"
Mesh get_invert_Mesh(const Mesh &mesh);
void divideBorderElems(AniMesh& am);

#endif //AORTIC_VALVE_TRIANGULARMESHHELPERS_H
