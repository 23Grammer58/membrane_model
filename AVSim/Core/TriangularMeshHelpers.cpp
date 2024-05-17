//
// Created by alex on 14.06.2020.
//

#include "TriangularMeshHelpers.h"
#include "AVMeshConnectivityHelpers.h"
#include <valarray>
#include <iterator>

using namespace std;

AniMesh convert_to_animesh(const Mesh &mesh, std::string v_lbl, std::string f_lbl, std::string e_lbl) {
    using namespace std;
    int nV = mesh.num_vertices(), nF = mesh.num_faces(), nE = mesh.num_edges();
    vector<double> v; v.reserve(nV*3);
    vector<int> f; f.reserve(3*nF);
    vector<int> e; e.reserve(2*nE);
    for (auto vert: mesh.vertices()){
        for (int i = 0; i < 3; ++i)
            v.push_back(mesh.point(vert)[i]);
    }
    for (auto face: mesh.faces()){
        auto v = vert_around(mesh, face);
        for (auto i: v)
            f.push_back(i+1);
    }
    for (auto edge: mesh.edges()){
        auto v = vert_around(mesh, edge);
        for (auto i: v)
            e.push_back(i+1);
    }
    std::vector<int> vertex_label(nV);
    std::vector<int> face_label(nF);
    std::vector<int> edge_label(nE);
    auto vlp = mesh.property_map<V_ind, int>(v_lbl);
    auto flp = mesh.property_map<F_ind, int>(f_lbl);
    auto elp = mesh.property_map<E_ind, int>(e_lbl);
    if (vlp.second) {
        copy(vlp.first.begin(), vlp.first.end(), vertex_label.begin());
    }
    if (flp.second) {
        copy(flp.first.begin(), flp.first.end(), face_label.begin());
    }
    if (elp.second) {
        copy(elp.first.begin(), elp.first.end(), edge_label.begin());
    }
    AniMesh am{v, f, e, vertex_label, face_label, edge_label};
    return am;
}

Mesh makeMesh(std::vector<array<double, 3>>& points, std::vector<std::vector<std::size_t>>& polygons){
    typedef std::vector<std::size_t> CGAL_Polygon;
    Mesh m;
    namespace PMP = CGAL::Polygon_mesh_processing;
    struct Array_traits
    {
        struct Equal_3
        {
            bool operator()(const array<double, 3>& p, const array<double, 3>& q) const {
                return (p == q);
            }
        };
        struct Less_xyz_3
        {
            bool operator()(const array<double, 3>& p, const array<double, 3>& q) const {
                return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
            }
        };
        Equal_3 equal_3_object() const { return Equal_3(); }
        Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
    };
    PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(Array_traits()));
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, m);
    return m;
}

Mesh split_twice(const Mesh& mesh){
    typedef std::vector<std::size_t> CGAL_Polygon;
    std::vector<CGAL_Polygon> polygons;
    std::vector<array<double, 3>> points;
    std::map<E_ind, int> edge_remap;
    std::map<V_ind, int> vert_remap;
    polygons.resize(mesh.num_faces()*4);
    points.reserve(mesh.num_vertices() + mesh.num_edges());
    for (auto v: mesh.vertices()){
        auto& p = mesh.point(v);
        points.push_back({p[0], p[1], p[2]});
        vert_remap[v] = points.size()-1;
    }
    for (auto& e: mesh.edges()) {
        auto v = vert_around(mesh, e);
        auto p0 = mesh.point(v[0]), p1 = mesh.point(v[1]);
        points.push_back({(p0[0] + p1[0])/2, (p0[1] + p1[1])/2, (p0[2] + p1[2])/2});
        edge_remap[e] = points.size()-1;
    }
    int nf = mesh.num_faces();
    int fi = 0;
    for (auto& f: mesh.faces()){
        auto e = edge_around(mesh, f);
        auto vv = vert_around(mesh, f);
//        std::array<std::array<V_ind, 2>, 3> ev;
//        for (int i = 0; i < 3; ++i) ev[i] = vert_around(mesh, e[i]);
//        assert(std::min(vert_around(mesh, e[0])[0], vert_around(mesh, e[0])[1]) == std::min(vv[0], vv[2])
//            && std::max(vert_around(mesh, e[0])[0], vert_around(mesh, e[0])[1])  == std::max(vv[0], vv[2])
//            && std::min(vert_around(mesh, e[1])[0], vert_around(mesh, e[1])[1]) == std::min(vv[0], vv[1])
//            && std::max(vert_around(mesh, e[1])[0], vert_around(mesh, e[1])[1])  == std::max(vv[0], vv[1]));
        std::array<std::size_t, 6> v;
        for (int i = 0; i < 3; ++i) v[(i+2)%3] = vert_remap[vv[i]];
        for (int i = 0; i < 3; ++i) v[i+3] = edge_remap[e[i]];
        polygons[fi] = std::vector<std::size_t>({v[3], v[4], v[5]});
        polygons[nf + fi] = std::vector<std::size_t>({v[0], v[5], v[4]});
        polygons[2*nf + fi] = std::vector<std::size_t>({v[1], v[3], v[5]});
        polygons[3*nf + fi] = std::vector<std::size_t>({v[2], v[4], v[3]});
        fi++;
    }
    return makeMesh(points, polygons);
}

Mesh convert_to_Mesh(const AniMesh &am, std::string v_lbl, std::string f_lbl, std::string e_lbl) {
    typedef std::vector<std::size_t>   CGAL_Polygon;
    std::vector<CGAL_Polygon> polygons;
    polygons.resize(am.faces.size() / 3);
    for (int i = 0; i < am.faces.size() / 3; ++i) {
        polygons[i].resize(3);
        for (int k = 0; k < 3; ++k) polygons[i][k] = am.faces[3 * i + k] - 1;
    }
    std::vector<array<double, 3>> points(am.vertices.size()/3);
    if (am.vertices.size()/3 > 0)
        std::memcpy(points[0].data(), am.vertices.data(), am.vertices.size() * sizeof(double));
    Mesh m = makeMesh(points, polygons);


//    for (unsigned i = 0; i < am.vertices.size()/3; ++i)
//        m.add_vertex(Point(am.vertices[3*i], am.vertices[3*i+1], am.vertices[3*i+2]));
//    for (unsigned i = 0; i < am.faces.size()/3; ++i)
//        if (m.add_face(static_cast<V_ind>(am.faces[3*i]-1), static_cast<V_ind>(am.faces[3*i+1]-1), static_cast<V_ind>(am.faces[3*i+2]-1)) == Mesh::null_face())
//            throw std::runtime_error("Can't add face");
    if (v_lbl != "" && am.vertex_label.size() == am.vertices.size()/3){
        auto it = m.add_property_map<V_ind, int>(v_lbl);
        for (unsigned i = 0; i < am.vertex_label.size(); ++i)
            it.first[static_cast<V_ind>(i)] = am.vertex_label[i];
    }
    if (f_lbl != "" && am.face_label.size() == am.faces.size()/3){
        auto it = m.add_property_map<F_ind, int>(f_lbl);
        for (unsigned i = 0; i < am.face_label.size(); ++i)
            it.first[static_cast<F_ind>(i)] = am.face_label[i];
    }
    if (e_lbl != "" && am.edge_label.size() == am.edges.size()/3){
        auto it = m.add_property_map<E_ind, int>(e_lbl);
        for (unsigned i = 0; i < am.edge_label.size(); ++i) {
            E_ind eid = m.edge(m.halfedge(static_cast<V_ind>(am.edges[2*i]-1), static_cast<V_ind>(am.edges[2*i+1]-1)));
            it.first[eid] = am.edge_label[i];
        }
    }

    return m;
}

void divideBorderElems(AniMesh& am){
    set<std::array<int, 2>> bnd;
    auto& edge = am.edges;
    for (int i = 0; i < am.edges.size()/2; ++i)
        bnd.insert({edge[2*i], edge[2*i+1]}), bnd.insert({edge[2*i+1], edge[2*i]});
    struct BndPatch{
        std::array<int, 2> fi;
        int com;
    };
    std::vector<BndPatch> wf;

    for (int f = 0; f < am.faces.size()/3; ++f){
        int ebnd = 0;
        int eind[3] = {0};
        for (int k = 0; k < 3; ++k)
            if (bnd.count({am.faces[3*f+k%3], am.faces[3*f+(k + 1)%3]})) {
                eind[ebnd++] = k;
            }
        if (ebnd >= 2) {
            wf.push_back({{f, -1}, -1});
            if (eind[0] == 0 && eind[1] == 1) wf.back().com = 1;
            else if (eind[0] == 1 && eind[1] == 2) wf.back().com  = 2;
            else if (eind[0] == 0 && eind[1] == 2) wf.back().com  = 0;
        }
    }

    for (int i = 0; i < am.faces.size()/3; ++i){
        std::set<int> lf{am.faces[3*i], am.faces[3*i+1], am.faces[3*i+2]};
        for (int j = 0; j < wf.size(); ++j){
            auto& bb = wf[j];
            if (    lf.count(am.faces[bb.fi[0]*3 + (bb.com + 1)%3])
                    &&  lf.count(am.faces[bb.fi[0]*3 + (bb.com + 2)%3])
                    && !lf.count(am.faces[bb.fi[0]*3 + (bb.com + 0)%3])){
                bb.fi[1] = i;
            }
        }
    }
    for (int l = 0; l < wf.size(); ++l){
        auto& bb = wf[l];
        int lind[4] = {0};
        for (int k = 0; k < 3; ++k){
            lind[k] = am.faces[bb.fi[0]*3 + (bb.com + k)%3];
        }
        for (int k = 0; k < 3; ++k) {
            if (am.faces[bb.fi[1]*3 + k] != lind[1] && am.faces[bb.fi[1]*3 + k] != lind[2])
                lind[3] = am.faces[3 * bb.fi[1] + k];
        }
        am.faces[3 * bb.fi[0] + 0] = lind[0], am.faces[3 * bb.fi[0] + 1] = lind[1], am.faces[3 * bb.fi[0] + 2] = lind[3];
        am.faces[3 * bb.fi[1] + 0] = lind[0], am.faces[3 * bb.fi[1] + 1] = lind[3], am.faces[3 * bb.fi[1] + 2] = lind[2];
    }
}

void translate_Mesh(Mesh &mesh, Vector translate) {
    for (auto v: mesh.vertices())
        mesh.points()[v] += translate;
}

Mesh get_invert_Mesh(const Mesh &mesh) {
    Mesh new_mesh;
    map<V_ind, V_ind> vv;
    for (auto v: mesh.vertices())
        vv.insert({v, new_mesh.add_vertex(mesh.point(v))});
    for (auto f: mesh.faces()){
        auto v = vert_around(mesh, f);
        new_mesh.add_face(v[0], v[2], v[1]);
    }
    auto vlb = mesh.property_map<V_ind, int>("v:boundary_lbl");
    if (vlb.second){
        auto it = new_mesh.add_property_map<V_ind, int>("v:boundary_lbl");
        for (auto& v: new_mesh.vertices())
            it.first[v] = vlb.first[v];
    }

    return new_mesh;
}

Mesh::Property_map<V_ind, int> remark_boundary(Mesh &m, int bnd_mrk, std::string v_lbl) {
    auto it = m.add_property_map<V_ind, int>(v_lbl);
    Mesh::Property_map<V_ind, int> lbl = it.first;
    std::fill(lbl.begin(), lbl.end(), 0);
    for (auto e: m.edges()){
        auto f = face_around(m, e);
        if (!f[1].second){
            auto v = vert_around(m, e);
            lbl[v[0]] = lbl[v[1]] = bnd_mrk;
        }
    }
    return lbl;
}

void set_edges_animesh(AniMesh &am) {
    std::set<std::pair<int, int>> ee;
    for (int i = 0; i < am.faces.size()/3; ++i)
        for (int k = 0; k < 3; ++k) {
            std::pair<int, int> edge = {am.faces[3*i + k], am.faces[3*i + (k + 1)%3]};
            if (edge.first > edge.second) std::swap(edge.first, edge.second);
            ee.insert(edge);
        }
    am.edges.resize(0);
    am.edges.reserve(ee.size()*2);
    for (auto& i: ee) {
        am.edges.push_back(i.first);
        am.edges.push_back(i.second);
    }
    return ;
}

bool write_STL(std::vector<CGAL::Simple_cartesian<double>::Triangle_3> triags, std::string filename) {
    ofstream ob(filename, ios::binary | ios::trunc);
    if (!ob) return false;
    array<uint8_t, 80> zero = {0};
    ob.write(reinterpret_cast<const char *> ((zero.data())), 80);
    uint32_t nt = triags.size();
    ob.write(reinterpret_cast<const char *> ((&nt)), sizeof(uint32_t));
    array<float, 3> v[4];
    for (auto &t: triags) {
        auto n = CGAL::cross_product(t[1] - t[0], t[2] - t[0]);
        for (int i = 0; i < 3; ++i)
            v[0][i] = n[i];
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                v[j + 1][i] = t[j][i];
        for (int j = 0; j < 4; ++j) {
            ob.write(reinterpret_cast<const char *> ((v[j].data())), 3 * sizeof(float));
        }
        uint16_t zz = 0;
        ob.write(reinterpret_cast<const char *> ((&zz)), sizeof(uint16_t));
    }
    ob.close();
    return true;
}

void invertFaceOrientation(AniMesh &am) { for (int i = 0; i < am.faces.size() / 3; ++i) std::swap(am.faces[3*i + 1], am.faces[3*i + 2]); }


