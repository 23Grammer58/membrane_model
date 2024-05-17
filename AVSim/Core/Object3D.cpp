//
// Created by alex on 13.06.2020.
//

#include "Object3D.h"
#include "TriangularMeshHelpers.h"

#include <cstdlib>
#include <fcntl.h>

#include <Eigen/Dense>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/IO/OFF.h>

#include "IO/Obj3dVTKSaver.h"

#ifdef USE_ANI3D
extern "C" {
#include "libaft.h"
}
#endif

using namespace std;

namespace World3d {

    static auto set_force_label(Mesh &m) {
        auto res = m.add_property_map<V_ind, Vector>("v:force");
        if (res.second) {
            for (auto v: m.vertices())
                res.first[v] = CGAL::NULL_VECTOR;
        }
        return res.first;
    }

    static auto set_next_x_label(Mesh &m) {
        auto res = m.add_property_map<V_ind, Point>("v:next_x");
        if (res.second) {
            for (auto i: m.vertices())
                res.first[i] = m.point(i);
        }
        return res.first;
    }

    static auto set_normal_label(Mesh &m) {
        auto res = m.add_property_map<F_ind, Vector>("f:normal");
        if (res.second) {
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                Vector v1(m.point(v[0]), m.point(v[1]));
                Vector v2(m.point(v[0]), m.point(v[2]));
                DReal v1_m = 0, v2_m = 0;
                for (int i = 0; i < 3; ++i){
                    if (fabs(v1[i]) > v1_m) v1_m = fabs(v1[i]);
                    if (fabs(v2[i]) > v2_m) v2_m = fabs(v2[i]);
                }
                v1 /= v1_m;
                v2 /= v2_m;

                auto n = CGAL::cross_product(v1, v2);
                n /= sqrt(n.squared_length());
                res.first[f] = n;
            }
        }
        return res.first;
    }

    static auto set_x_label(Mesh &m) {
        auto res = m.add_property_map<V_ind, Point>("v:x");
        if (res.second) {
            for (auto i: m.vertices())
                res.first[i] = m.point(i);
        }
        return res.first;
    }

    static auto set_x0_label(Mesh &m) {
        auto res = m.property_map<V_ind, Point>("v:point");
        return res.first;
    }

    static auto set_boundary_label(Mesh &m) {
        auto res = m.add_property_map<V_ind, int>("v:boundary_lbl");
        if (res.second) {
            for (auto i: m.vertices())
                res.first[i] = 1;
        }
        return res.first;
    }

    void Object3D::repair_connectivity(){
        auto it1 = m_mesh.add_property_map<V_ind, Vector>("v:force");
        auto it2 = m_mesh.add_property_map<V_ind, Point>("v:next_x");
        auto it3 = m_mesh.add_property_map<V_ind, Point>("v:x");
        auto it6 = m_mesh.property_map<V_ind, Point>("v:point");
        auto it4 = m_mesh.add_property_map<V_ind, int>("v:boundary_lbl");
        auto it5 = m_mesh.add_property_map<F_ind, Vector>("f:normal");
        for (auto v: m_mesh.vertices()){
            it1.first[v] = m_F[v];
            it2.first[v] = m_next_x[v];
            it3.first[v] = m_x[v];
            it4.first[v] = m_boundary[v];
            it6.first[v] = m_x0[v];
        }
        for (auto f: m_mesh.faces()){
            it5.first[f] = m_normal[f];
        }

        m_F = std::move(it1.first);
        m_next_x = std::move(it2.first);
        m_x = std::move(it3.first);
        m_x0 = std::move(it6.first);
        m_boundary = std::move(it4.first);
        m_normal.first = std::move(it5.first);
        auto normal_updater = [](Object3D &obj) -> int {
            auto &m = obj.m_mesh;
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                Vector v1(obj.m_x[v[0]], obj.m_x[v[1]]);
                Vector v2(obj.m_x[v[0]], obj.m_x[v[2]]);
                DReal v1_m = 0, v2_m = 0;
                for (int i = 0; i < 3; ++i){
                    if (fabs(v1[i]) > v1_m) v1_m = fabs(v1[i]);
                    if (fabs(v2[i]) > v2_m) v2_m = fabs(v2[i]);
                }
                v1 /= v1_m;
                v2 /= v2_m;

                auto n = CGAL::cross_product(v1, v2);
                n /= sqrt(n.squared_length());

                obj.m_normal.first[f] = n;
            }
            return 0;
        };
        m_normal.second = add_updater(normal_updater);
    }

    static void copy_labels_from_mesh(const Mesh &src, Mesh &dest){
        assert( src.num_vertices() == dest.num_vertices()
                && src.num_edges() == dest.num_edges()
                && src.num_faces() == dest.num_faces());
        auto l1d = dest.add_property_map<V_ind, Vector>("v:force");
        auto l1s = src.property_map<V_ind, Vector>("v:force");
        auto l2d = dest.add_property_map<V_ind, Point>("v:next_x");
        auto l2s = src.property_map<V_ind, Point>("v:next_x");
        auto l4d = dest.add_property_map<V_ind, Point>("v:x");
        auto l4s = src.property_map<V_ind, Point>("v:x");
        auto l5d = dest.add_property_map<V_ind, int>("v:boundary_lbl");
        auto l5s = src.property_map<V_ind, int>("v:boundary_lbl");
        for (auto v: src.vertices()){
            l1d.first[v] = l1s.first[v];
            l2d.first[v] = l2s.first[v];
            l4d.first[v] = l4s.first[v];
            l5d.first[v] = l5s.first[v];
        }
        auto l3d = dest.add_property_map<F_ind, Vector>("f:normal");
        auto l3s = src.property_map<F_ind, Vector>("f:normal");
        for (auto f: src.faces()){
            l3d.first[f] = l3s.first[f];
        }
    }


    static auto Mesh_init(Mesh &m) {
        return tuple<   Mesh::Property_map<V_ind, Vector>,
                        Mesh::Property_map<V_ind, Point>,
                        Mesh::Property_map<F_ind, Vector>,
                        Mesh::Property_map<V_ind, Point>,
                        Mesh::Property_map<V_ind, Point>,
                        Mesh::Property_map<V_ind, int>>
                    {set_force_label(m),
                     set_next_x_label(m),
                     set_normal_label(m),
                     set_x_label(m),
                     set_x0_label(m),
                     set_boundary_label(m)};
    }

    void Object3D::reset_mesh() {
        auto t = Mesh_init(m_mesh);
        m_F = std::move(get<0>(t));
        m_next_x = std::move(get<1>(t));
        m_x = std::move(get<3>(t));
        m_x0 = std::move(get<4>(t));
        m_boundary = std::move(get<5>(t));
        m_normal.first = std::move(get<2>(t));
        auto normal_updater = [](Object3D &obj) -> int {
            auto &m = obj.m_mesh;
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                Vector v1(obj.m_x[v[0]], obj.m_x[v[1]]);
                Vector v2(obj.m_x[v[0]], obj.m_x[v[2]]);
                DReal v1_m = 0, v2_m = 0;
                for (int i = 0; i < 3; ++i){
                    if (fabs(v1[i]) > v1_m) v1_m = fabs(v1[i]);
                    if (fabs(v2[i]) > v2_m) v2_m = fabs(v2[i]);
                }
                v1 /= v1_m;
                v2 /= v2_m;

                auto n = CGAL::cross_product(v1, v2);
                n /= sqrt(n.squared_length());

                obj.m_normal.first[f] = n;
            }
            return 0;
        };
        m_normal.second = add_updater(normal_updater);
    }

    void Object3D::init() {
        reset_mesh();
    }

    Object3D::Object3D(Object3D&& obj):
        name{std::move(obj.name)},
        __id{std::move(obj.__id)},
        m_mesh{std::move(obj.m_mesh)},
        m_forces(std::move(obj.m_forces)),
        _next_x_updater(std::move(obj._next_x_updater)),
        m_apply_next_x(std::move(obj.m_apply_next_x)),
        m_normal{std::move(obj.m_normal)},
        m_F{std::move(obj.m_F)},
        m_x0{std::move(obj.m_x0)},
        m_x(std::move(obj.m_x)),
        m_next_x(std::move(obj.m_next_x)),
        m_boundary(std::move(obj.m_boundary)),
        _renderer{std::move(obj._renderer)},
        _bc{std::move(obj._bc)}
    {
        repair_connectivity();
    }
    Object3D::Object3D(const Object3D& obj):
        name{obj.name},
        __id{obj.__id},
        m_mesh{obj.m_mesh},
        m_forces(obj.m_forces),
        m_apply_next_x(obj.m_apply_next_x),
        _next_x_updater(obj._next_x_updater),
        m_normal{obj.m_normal},
        m_F{obj.m_F},
        m_x0{obj.m_x0},
        m_x(obj.m_x),
        m_next_x(obj.m_next_x),
        m_boundary(obj.m_boundary),
        _renderer{obj._renderer},
        _bc{obj._bc}
    {
        repair_connectivity();
    }
    Object3D& Object3D::operator=(Object3D&& obj){
        if (&obj == this) return *this;
        name = std::move(obj.name);
        __id = std::move(obj.__id);
        m_mesh = std::move(obj.m_mesh);
        m_forces = std::move(obj.m_forces);
        m_apply_next_x = std::move(obj.m_apply_next_x);
        _next_x_updater = std::move(obj._next_x_updater);
        m_normal = std::move(obj.m_normal);
        m_F = std::move(obj.m_F);
        m_x0 = std::move(obj.m_x0);
        m_x = std::move(obj.m_x);
        m_next_x = std::move(obj.m_next_x);
        m_boundary = std::move(obj.m_boundary);
        _renderer = std::move(obj._renderer);
        _bc = std::move(obj._bc);
        repair_connectivity();

        return *this;
    }
    Object3D& Object3D::operator=(const Object3D& obj){
        if (&obj == this) return *this;
        name = obj.name;
        __id = obj.__id;
        m_mesh = obj.m_mesh;
        m_forces = obj.m_forces;
        m_apply_next_x = obj.m_apply_next_x;
        _next_x_updater = obj._next_x_updater;
        m_normal = obj.m_normal;
        m_F = obj.m_F;
        m_x0 = obj.m_x0;
        m_x = obj.m_x;
        m_next_x = obj.m_next_x;
        m_boundary = obj.m_boundary;
        _renderer = obj._renderer;
        _bc = obj._bc;
        repair_connectivity();

        return *this;
    }

    static bool write_STL(std::string filename, Mesh &mesh) {
        ofstream ob(filename, ios::binary | ios::trunc);
        if (!ob)
            return false;
        auto n = mesh.property_map<F_ind, Vector>("f:normal");
        if (!n.second) {
            cout << "The mesh doesn't have normal tag!";
            return false;
        }
        array<uint8_t, 80> zero = {0};
        ob.write(reinterpret_cast<const char *> ((zero.data())), 80);
        uint32_t nt = mesh.num_faces();
        ob.write(reinterpret_cast<const char *> ((&nt)), sizeof(uint32_t));
        array<float, 3> v[4];
        for (auto &f: mesh.faces()) {
            for (int i = 0; i < 3; ++i)
                v[0][i] = n.first[f][i];
            auto itt = vertices_around_face(mesh.halfedge(f), mesh);
            if (itt.size() != 3)
                throw runtime_error("write_STL() is implemented only for triangulations now\n");
            auto it = itt.first;
            for (int j = 0; j < 3; ++j) {
                auto p = mesh.point(*(it++));
                for (int i = 0; i < 3; ++i)
                    v[j + 1][i] = p[i];
            }
            for (int j = 0; j < 4; ++j) {
                ob.write(reinterpret_cast<const char *> ((v[j].data())), 3 * sizeof(float));
            }
            uint16_t zz = 0;
            ob.write(reinterpret_cast<const char *> ((&zz)), sizeof(uint16_t));
        }
        ob.close();
        return true;
    }

    static bool write_custom(std::string filename, vector<Mesh> &vmesh) {
        ofstream ob(filename);
        if (!ob)
            return false;
        ob << vmesh.size() << "\n";
        for (unsigned i = 0; i < vmesh.size(); ++i) {
            auto &m = vmesh[i];
            ob << m.num_vertices() << " " << m.num_faces() << "\n";
            auto blb_p = vmesh[i].property_map<V_ind, int>("v:boundary_lbl");
            auto norm_p = vmesh[i].property_map<F_ind, Vector>("f:normal");
            for (auto v: m.vertices()) {
                std::array<double, 3> p;
                for (int i = 0; i < 3; ++i)
                    p[i] = m.point(v)[i];
                ob << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                    << p[0] << " " << p[1] << " " << p[2] << " ";
                if (blb_p.second)
                    ob << blb_p.first[v] << "\n";
                else
                    ob << 0 << "\n";
            }
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                Vector n;
                if (norm_p.second) {
                    n = norm_p.first[f];
                } else {
                    n = CGAL::cross_product(Vector(m.point(v[0]), m.point(v[1])),
                                            Vector(m.point(v[0]), m.point(v[2])));
                    n /= sqrt(n.squared_length());
                }
                ob << v[0] + 1 << " " << v[1] + 1 << " " << v[2] + 1 << " ";
                std::array<double, 3> np;
                for (int i = 0; i < 3; ++i)
                    np[i] = n[i];
                ob << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                   << np[0] << " " << np[1] << " " << np[2] << "\n";
            }
        }
        ob.close();
        return true;
    }

    static bool read_custom(std::string filename, vector<Mesh> &vmesh) {
        ifstream in(filename);
        if (!in)
            return false;
        int nm = 1;
        {
            std::string line;
            getline(in, line);
            auto pos = line.find_first_of(" ");
            while (pos == line.size() - 1) {
                line.erase(pos, 1);
            }
            if (pos == line.npos)
                nm = stoi(line);
            else in.seekg(0);
        }
        vmesh.resize(nm);
//        vector<V_ind> vertices;
//        for (int i = 0; i < nm; ++i) {
//            vertices.clear();
//            auto lbl_p = vmesh[i].add_property_map<V_ind, int>("v:boundary_lbl");
//            auto norm_p = vmesh[i].add_property_map<F_ind, Vector>("f:normal");
//            if (!lbl_p.second || !norm_p.second) return false;
//            auto lbl = lbl_p.first;
//            auto norm = norm_p.first;
//            int nv = 0, nt = 0;
//            in >> nv >> nt;
//            for (int v = 0; v < nv; ++v) {
//                double x, y, z;
//                int l;
//                in >> x >> y >> z >> l;
//                auto id = vmesh[i].add_vertex(Point(x, y, z));
//                lbl[id] = l;
//                vertices.push_back(id);
//            }
//            for (int f = 0; f < nt; ++f) {
//                int v1, v2, v3;
//                double n1, n2, n3;
//                in >> v1 >> v2 >> v3 >> n1 >> n2 >> n3;
//                --v1;
//                --v2;
//                --v3;
//                Vector normal = Vector(n1, n2, n3);
//                if (CGAL::cross_product(Vector(vmesh[i].point(vertices[v1]), vmesh[i].point(vertices[v2])),
//                                        Vector(vmesh[i].point(vertices[v1]), vmesh[i].point(vertices[v3]))) *
//                    normal < 0)
//                    swap(v3, v2);
//                auto id = vmesh[i].add_face(vertices[v1], vertices[v2], vertices[v3]);
//                norm[id] = normal;
//            }
//        }
        vector<array<double, 3>> points;
        vector<int> lbl;
        std::vector<std::vector<std::size_t>> polygons;
//        std::vector<Vector> normals;
        for (int i = 0; i < nm; ++i) {
            int nv = 0, nt = 0;
            in >> nv >> nt;
            points.resize(nv); lbl.resize(nv);
            for (int v = 0; v < nv; ++v)
                in >> points[v][0] >> points[v][1] >> points[v][2] >> lbl[v];

            polygons.resize(nt);
//            normals.resize(nt);
            auto get_p = [&points](int v) { return Point(points[v][0], points[v][1], points[v][2]); };
            for (int f = 0; f < nt; ++f) {
                int v1=0, v2=0, v3=0;
                double n1=0, n2=0, n3=0;
                in >> v1 >> v2 >> v3 >> n1 >> n2 >> n3;
                --v1; --v2; --v3;
                Vector normal = Vector(n1, n2, n3);
                if (CGAL::cross_product(Vector(get_p(v1), get_p(v2)),
                                        Vector(get_p(v1), get_p(v3))) *
                    normal < 0)
                    swap(v3, v2);
                polygons[f].resize(3);
                polygons[f][0] = v1, polygons[f][1] = v2, polygons[f][2] = v3;
//                normals[f] = normal;
            }
            vmesh[i] = std::move(makeMesh(points, polygons));
            if (vmesh[i].num_vertices() != nv) throw std::runtime_error("Wrong number of vertices");
            if (vmesh[i].num_faces() != nt) throw std::runtime_error("Wrong number of faces");
            auto lbl_p = vmesh[i].add_property_map<V_ind, int>("v:boundary_lbl");
            auto norm_p = vmesh[i].add_property_map<F_ind, Vector>("f:normal");
            if (!lbl_p.second || !norm_p.second) return false;
            for (auto v: vmesh[i].vertices()) lbl_p.first[v] = lbl[v.idx()];
            for (auto f: vmesh[i].faces()) {
                auto v = vert_around(vmesh[i], f);
                Vector lnorm = CGAL::cross_product(vmesh[i].point(v[1]) - vmesh[i].point(v[0]), vmesh[i].point(v[2]) - vmesh[i].point(v[0]));
                lnorm /= sqrt(lnorm.squared_length());
//                if ((lnorm - normals[f.idx()]).squared_length() > 1e-6) throw std::runtime_error("Wrong faces numeration");
                norm_p.first[f] = lnorm;//normals[f.idx()];
            }
        }
        return true;
    }

    static bool read_STL(std::string filename, Mesh &mesh) {
        ifstream in(filename);
        if (!in)
            return false;

        std::string line;
        getline(in, line);
        const char *c = line.c_str();
        bool bin = false;
        if (strncmp(c, "solid ", 6))
            bin = true;

        struct Comp {
            bool operator()(const array<float, 3> &a1,
                            const array<float, 3> &a2) const {
                for (unsigned i = 0; i < a1.size(); ++i)
                    if (a1[i] != a2[i])
                        return a1[i] < a2[i];
                return false;
            }
        };
        auto normals = mesh.add_property_map<F_ind, Vector>("f:normal");

        std::vector<array<float, 3>> points;
        typedef std::vector<std::size_t>   CGAL_Polygon;
        std::vector<CGAL_Polygon> polygons;
        namespace PMP = CGAL::Polygon_mesh_processing;
        struct Array_traits
        {
            struct Equal_3
            {
                bool operator()(const array<float, 3>& p, const array<float, 3>& q) const {
                    return (p == q);
                }
            };
            struct Less_xyz_3
            {
                bool operator()(const array<float, 3>& p, const array<float, 3>& q) const {
                    return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
                }
            };
            Equal_3 equal_3_object() const { return Equal_3(); }
            Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
        };

        auto process_triangle = [&points, &polygons](array<float, 3>& v0, array<float, 3>& v1, array<float, 3>& v2, array<float, 3>& n){
            Vector nv[3] = {Vector(v0[0], v0[1], v0[2]),
                            Vector(v1[0], v1[1], v1[2]),
                            Vector(v2[0], v2[1], v2[2])};
            Vector norm(n[0], n[1], n[2]);
            if (CGAL::cross_product((nv[1] - nv[0]), (nv[2] - nv[0])) * norm < 0) {
                std::swap(v2, v1);
            }
            std::size_t id = points.size();
            points.push_back(v0);
            points.push_back(v1);
            points.push_back(v2);
            polygons.push_back({id, id+1, id+2});
        };

        if (bin) {
            in.close();
            in.open(filename, std::ios::binary);
            in.seekg(80, ios_base::beg);
            uint32_t nt = 0;
            in.read(reinterpret_cast<char *> ((&nt)), sizeof(uint32_t));

            array<float, 3> v[4];
            for (uint32_t i = 0; i < nt; ++i) {
                uint16_t color;
                for (int j = 0; j < 4; ++j)
                    in.read(reinterpret_cast<char *> ((v[j].begin())), 3 * sizeof(float));
                in.read(reinterpret_cast<char *> ((&color)), sizeof(uint16_t));
                process_triangle(v[1], v[2], v[3], v[0]);
            }
            in.close();
        } else {
            map<array<float, 3>, V_ind, Comp> connect;
            array<float, 3> n, v[3];
            std::vector<string> words;
            auto divide_line_to_words = [](string s, std::vector<string>& out){
                out.resize(0);
                std::stringstream ss(s);
                while (!ss.eof()){
                    out.resize(out.size()+1);
                    ss >> out[out.size()-1];
                }
            };

            while (!in.eof()) {
                getline(in, line);
                if (!strncmp(line.c_str(), "endsolid", 8)) break;
                divide_line_to_words(line, words);
                for (int i = 0; i < 3; ++i)
                    n[2 - i] = std::stof(words[words.size()-1-i]);
                in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (int j = 0; j < 3; ++j) {
                    getline(in, line);
                    divide_line_to_words(line, words);
                    for (int i = 0; i < 3; ++i)
                        v[j][2 - i] = std::stod(words[words.size()-1-i]);
                }
                in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                process_triangle(v[0], v[1], v[2], n);
            }

//            Vector nv[3] = {Vector(v[0][0], v[0][1], v[0][2]),
//                            Vector(v[1][0], v[1][1], v[1][2]),
//                            Vector(v[2][0], v[2][1], v[2][2])};
//            Vector norm(n[0], n[1], n[2]);
//            if (CGAL::cross_product((nv[1] - nv[0]), (nv[2] - nv[0])) * norm < 0)
//                std::swap(v[2], v[1]);
//
//            V_ind id[3];
//            for (int i = 0; i < 3; ++i) {
//                auto it = connect.find(v[i]);
//                if (it != connect.end())
//                    id[i] = it->second;
//                else {
//                    id[i] = mesh.add_vertex(Point(v[i][0], v[i][1], v[i][2]));
//                    connect.insert({v[i], id[i]});
//                }
//            }
//            auto f = mesh.add_face(id[0], id[1], id[2]);
//            normals.first[f] = Vector(n[0], n[1], n[2]);
            in.close();
        }
        PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(Array_traits()));
        PMP::orient_polygon_soup(points, polygons);
        PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
        for (auto f: mesh.faces()){
            auto v = vert_around(mesh, f);
            normals.first[f] = CGAL::cross_product(mesh.point(v[1]) - mesh.point(v[0]), mesh.point(v[2]) - mesh.point(v[0]));
            normals.first[f] /= sqrt(normals.first[f].squared_length());
        }

        return true;
    }
#ifdef USE_ANI3D
    static bool read_gmv(std::string filename, Mesh &mesh) {
        AniMesh am;
        int maxnV = 10000000, maxnF = 10000000;
        int nV = 0, nF = 0;
        am.vertices.resize(3 * maxnV);
        am.faces.resize(3 * maxnF);
        am.face_label.resize(maxnF);
        int status = read_front(const_cast<char *> ((filename).c_str()), &nV, am.vertices.data(), &nF, am.faces.data(),
                                am.face_label.data(),
                                maxnV, maxnF, 1, 0, 0);
        if (status != 0)
            return false;
        am.vertices.resize(3 * nV);
        am.faces.resize(3 * nF);
        am.face_label.resize(nF);
        mesh = convert_to_Mesh(am, "", "f:ani_label", "");
        return true;
    }

    static bool write_gmv(std::string filename, Mesh &mesh) {
        AniMesh am = convert_to_animesh(mesh);
        int st = write_front_gmv(const_cast<char *> ((filename).c_str()), am.vertices.size(), am.vertices.data(),
                                 am.faces.size(), am.faces.data(), am.face_label.data());
        return st == 0;
    }

    static bool write_smv(std::string filename, Mesh &mesh) {
        AniMesh am = convert_to_animesh(mesh);
        int st = write_front(const_cast<char *> ((filename).c_str()), am.vertices.size() / 3, am.vertices.data(),
                             am.faces.size() / 3, am.faces.data(), am.face_label.data());
        return st == 0;
    }
#endif
    bool Object3D::read(std::string filename) {
        std::string::size_type dot(filename.rfind("."));
        if (dot == std::string::npos) throw runtime_error("Can't define file extension for\"" + filename + "\"");
        std::string ext = filename.substr(dot + 1, filename.length() - dot - 1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        Mesh mesh;
//        Mesh_init(mesh);
        bool status = false;
        if (ext == "ply") {
            ifstream in(filename);
            std::string comment;
            status = CGAL::IO::read_PLY(in, mesh, comment);
        }
#ifdef USE_ANI3D
        else if (ext == "gmv")
            status = read_gmv(filename, mesh);
#endif
        else if (ext == "off"){
            ifstream in(filename);
            status = CGAL::IO::read_OFF(in, mesh);
        } else if (ext == "stl")
            status = read_STL(filename, mesh);
        else if (ext == "txt") {
            vector<Mesh> vm;
            status = read_custom(filename, vm);
            mesh = vm[0];
        }
//    else if (ext == "3mf") {
//        vector<Mesh> vm;
//         status = CGAL::read_3mf(filename, vm);
//         mesh = vm[0];
//    }
        else
            throw runtime_error("Unknown extension " + ext);

        if (status) {
            m_mesh = std::move(mesh);
            reset_mesh();
        }

        return status;
    }

    void save_mesh(std::string filename, Mesh& m){
        std::string::size_type dot(filename.rfind("."));
        if (dot == std::string::npos) throw runtime_error("Can't define file extension for\"" + filename + "\"");
        std::string ext = filename.substr(dot + 1, filename.length() - dot - 1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == "ply") {
            ofstream in(filename);
            CGAL::IO::write_PLY(in, m, "");
        } else if (ext == "off"){
            ofstream in(filename);
            CGAL::IO::write_OFF(in, m);
        } else if (ext == "stl")
            write_STL(filename, m);
#ifdef USE_ANI3D
        else if (ext == "gmv")
            write_gmv(filename, m);
        else if (ext == "smv")
            write_smv(filename, m);
#endif
        else if (ext == "txt") {
            vector<Mesh> vm(1);
            vm[0] = m;
            write_custom(filename, vm);
        } else
            throw runtime_error("Unknown extension " + ext);
    }

    void Object3D::save(std::string file, const Mesh::Property_map<V_ind, Point>& view_x, const Mesh::Property_map<F_ind, Vector>& view_n) const{
        std::string::size_type dot(file.rfind("."));
        if (dot == std::string::npos) throw runtime_error("Can't define file extension for\"" + file + "\"");
        std::string ext = file.substr(dot + 1, file.length() - dot - 1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == "vtk"){
            VTKSaver vs;
            vs.setObj(*this);
            vs.setX(view_x);
            auto nodetg = m_mesh.properties<V_ind>();
            for (auto& n: nodetg) vs.setTagByName(n, n);
            auto facetg = m_mesh.properties<F_ind>();
            for (auto& n: facetg) vs.setTagByName(n, n);
            vs.save(file);
        } else {
            Mesh m = m_mesh;
            auto &x0 = m.points();
            auto it2x = view_x.begin();
            for (auto it1 = x0.begin(); it1 != x0.end() && it2x != view_x.end(); ++it1, ++it2x) {
                *it1 = *it2x;
            }
            auto t = Mesh_init(m);
            auto &n = get<2>(t);
            auto &lbl = get<5>(t);
            auto it2n = view_n.begin();
            for (auto it1 = n.begin(); it1 != n.end() && it2n != view_n.end(); ++it1, ++it2n) {
                *it1 = *it2n;
            }
            auto &view_lbl = m_boundary;
            auto it1 = lbl.begin();
            auto it2 = view_lbl.begin();
            for (; it1 != lbl.end() && it2 != view_lbl.end(); ++it1, ++it2) {
                *it1 = *it2;
            }
            save_mesh(file, m);
        }
    }

    void Object3D::save(std::string file, const Mesh::Property_map<V_ind, Point>& view_x) const{
        std::string::size_type dot(file.rfind("."));
        if (dot == std::string::npos) throw runtime_error("Can't define file extension for\"" + file + "\"");
        std::string ext = file.substr(dot + 1, file.length() - dot - 1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == "vtk"){
            VTKSaver vs;
            vs.setObj(*this);
            vs.setX(view_x);
            auto nodetg = m_mesh.properties<V_ind>();
            for (auto& n: nodetg) vs.setTagByName(n, n);
            auto facetg = m_mesh.properties<F_ind>();
            for (auto& n: facetg) vs.setTagByName(n, n);
            vs.save(file);
        } else {
            Mesh m = m_mesh;
            auto &x0 = m.points();
            auto it2x = view_x.begin();
            for (auto it1 = x0.begin(); it1 != x0.end() && it2x != view_x.end(); ++it1, ++it2x) {
                *it1 = *it2x;
            }
            auto t = Mesh_init(m);
            auto &n = get<2>(t);
            auto &lbl = get<5>(t);

            auto &view_lbl = m_boundary;
            auto it1 = lbl.begin();
            auto it2 = view_lbl.begin();
            for (; it1 != lbl.end() && it2 != view_lbl.end(); ++it1, ++it2) {
                *it1 = *it2;
            }
            save_mesh(file, m);
        }
    }

    void Object3D::save(std::string filename) const {
        save(filename, m_x, m_normal.first);
    }

    void Object3D::set_zero_residual(){ 
        for (auto v: m_mesh.vertices()) 
            m_F[v] = CGAL::NULL_VECTOR; 
    }

    Force_ID Object3D::add_force(Force f) {
        static int force_id = 0;
        m_forces.insert({force_id, f});
        return force_id++;
    }

    Updater_ID Object3D::add_updater(Processor f) {
        static int update_id = 0;
        m_updaters.insert({update_id, f});
        return update_id++;
    }

    int Object3D::apply_forces() {
        int status = 0;
        for (auto &i: m_forces)
            if ((status = i.second(*this)) < 0)
                return status;
        return status;
    }

    int Object3D::apply_force(int fid){
        return m_forces[fid](*this);
    }

    void Object3D::apply_updaters() {
        set_zero_residual();

        for (auto &i: m_updaters)
            i.second(*this);
    }

    void Object3D::apply_render() {
        if (_renderer)
            _renderer->operator()(*this);
    }

    void Object3D::remove_force(Force_ID force_id) {
        m_forces.erase(force_id);
    }

    void Object3D::remove_updater(Updater_ID updater_id) {
        m_updaters.erase(updater_id);
    }

    ConstProperty<E_ind, DReal> set_l0(Object3D &obj, bool forced) {
        auto it = obj.m_mesh.add_property_map<E_ind, DReal>("e:l0");
        if (forced || it.second) {
            auto &m = obj.m_mesh;
            for (auto e: m.edges()) {
                auto v = vert_around(m, e);
                DReal l0 = sqrt((obj.m_x0[v[1]] - obj.m_x0[v[0]]).squared_length());
                it.first[e] = l0;
            }
        }
        return it;
    }

    ConstProperty<F_ind, DReal> set_Ap(Object3D &obj, bool forced) {
        auto it = obj.m_mesh.add_property_map<F_ind, DReal>("f:Ap");
        if (forced || it.second) {
            auto &m = obj.m_mesh;
            auto &Ap = it.first;
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                Eigen::Matrix<DReal, 3, 1> P[3];
                for (int i = 0; i < 3; ++i) {
                    P[i] << obj.m_x0[v[i]][0], obj.m_x0[v[i]][1], obj.m_x0[v[i]][2];
                }
                Ap[f] = ((P[1] - P[0]).cross(P[2] - P[0])).stableNorm() / 2;
//                Ap[f] = 0.5 * sqrt((CGAL::cross_product(Vector(obj.m_x0[v[0]], obj.m_x0[v[1]]),
//                                                        Vector(obj.m_x0[v[0]], obj.m_x0[v[2]]))).squared_length());
            }
        }
        return it;
    }

    UpdatableProperty<E_ind, DReal> set_l(Object3D &obj) {
        auto it = obj.m_mesh.add_property_map<E_ind, DReal>("e:l");
        Updater_ID updaterId = -1;
        if (it.second) {
            auto l_updater = [l = it.first](Object3D &obj) {
                auto &m = obj.m_mesh;
                for (auto e: m.edges()) {
                    auto v = vert_around(m, e);
                    DReal _l = sqrt((obj.m_x[v[1]] - obj.m_x[v[0]]).squared_length());
                    l[e] = _l;
                }
                return 0;
            };
            l_updater(obj);
            updaterId = obj.add_updater(l_updater);
        }
        return UpdatableProperty<E_ind, DReal>{pair<Mesh::Property_map<E_ind, DReal>, Updater_ID>{it.first, updaterId}};
    }

    UpdatableProperty<F_ind, Vector> set_S(Object3D &obj) {
        auto it = obj.m_mesh.add_property_map<F_ind, Vector>("f:VectorArea");
        Updater_ID updaterId = -1;
        if (it.second) {
            auto S_updater = [S = it.first](Object3D &obj) -> int {
                auto &m = obj.m_mesh;
                for (auto f: m.faces()) {
                    auto v = vert_around(m, f);
                    Eigen::Matrix<DReal, 3, 1> Q[3], OrArea;
                    for (int i = 0; i < 3; ++i)
                        Q[i] << obj.m_x[v[i]][0], obj.m_x[v[i]][1], obj.m_x[v[i]][2];
                    OrArea = ((Q[1] - Q[0]).cross(Q[2] - Q[0]))/2;
                    S[f] = Vector(OrArea[0], OrArea[1], OrArea[2]);
//                    S[f] = 0.5 * CGAL::cross_product(Vector(obj.m_x[v[0]], obj.m_x[v[1]]),
//                                                     Vector(obj.m_x[v[0]], obj.m_x[v[2]]));
                }
                return 0;
            };
            S_updater(obj);
            updaterId = obj.add_updater(S_updater);
        }
        return UpdatableProperty<F_ind, Vector>{pair<Mesh::Property_map<F_ind, Vector>, Updater_ID>{it.first, updaterId}};
    }

    ConstProperty<F_ind, Vector> set_n0(Object3D &obj, bool forced) {
        auto it = obj.m_mesh.add_property_map<F_ind, Vector>("f:normal0");
        if (forced || it.second) {
            auto &m = obj.m_mesh;
            auto &n0 = it.first;
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                n0[f] = CGAL::cross_product(Vector(obj.m_x0[v[0]], obj.m_x0[v[1]]),
                                            Vector(obj.m_x0[v[0]], obj.m_x0[v[2]]));
                n0[f] /= sqrt(n0[f].squared_length());
            }
        }
        return it;
    }

    ConstProperty<F_ind, std::array<Vector, 3>> set_D_vecs(Object3D& obj, bool forced){
        auto it = obj.m_mesh.add_property_map<F_ind, std::array<Vector, 3>>("f:D_vecs");
        if (forced || it.second) {
            auto Ap = set_Ap(obj, forced);
            auto n0 = set_n0(obj, forced);
            auto &m = obj.m_mesh;
            auto &DD = it.first;
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                Eigen::Matrix<DReal, 3, 1> n;
                n << n0[f][0], n0[f][1], n0[f][2];
                Eigen::Matrix<DReal, 3, 3> P, Dp;
                P <<    obj.m_x0[v[0]][0], obj.m_x0[v[1]][0], obj.m_x0[v[2]][0],
                        obj.m_x0[v[0]][1], obj.m_x0[v[1]][1], obj.m_x0[v[2]][1],
                        obj.m_x0[v[0]][2], obj.m_x0[v[1]][2], obj.m_x0[v[2]][2];
                for (int i = 0; i < 3; ++i) {
                    Dp.col(i) = n.cross(P.col((i + 2) % 3) - P.col((i + 1) % 3));
                    assert(Dp.col(i).dot(P.col((i) % 3) - P.col((i + 1) % 3)) >= 0);
                }
                Dp /= 2 * Ap[f];
                for (int i = 0; i < 3; ++i)
                        DD[f][i] = Vector(Dp(0, i), Dp(1, i), Dp(2, i));
            }
            if (n0.second) obj.m_mesh.remove_property_map(n0.first);
            if (Ap.second) obj.m_mesh.remove_property_map(Ap.first);
        }
        return it;
    }

    ConstProperty<F_ind, std::array<DReal, 6>> set_DD_matrix(Object3D &obj, bool forced) {
        auto it = obj.m_mesh.add_property_map<F_ind, std::array<DReal, 6>>("f:DD_matrix");
        if (forced || it.second) {
            auto _D = set_D_vecs(obj, forced);
            auto &m = obj.m_mesh;
            auto &DD_m = it.first;
            for (auto f: m.faces()) {
                auto v = vert_around(m, f);
                Eigen::Matrix<DReal, 3, 1> D[3];
                for (int i = 0; i < 3; ++i)
                    D[i] << _D[f][i][0], _D[f][i][1], _D[f][i][2];
                {
                    DD_m[f][0] = D[0].squaredNorm();//D[0] * D[0];
                    DD_m[f][1] = D[0].dot(D[1]);//D[0] * D[1];
                    DD_m[f][2] = D[0].dot(D[2]);//-(DD_m[f][0] + DD_m[f][1]);
                    DD_m[f][3] = D[1].squaredNorm();//D[1] * D[1];
                    DD_m[f][4] = D[1].dot(D[2]);//-(DD_m[f][3] + DD_m[f][1]);
                    DD_m[f][5] = D[2].squaredNorm();//DD_m[f][0] + 2 * DD_m[f][1] + DD_m[f][3];
                }
            }
            if (_D.second) obj.m_mesh.remove_property_map(_D.first);
        }
        return ConstProperty < F_ind, std::array<DReal, 6>>
        (std::move(it));
    }

    void Object3D::apply_next_x() {
        m_apply_next_x(*this);
    }

    double Object3D::residual(unsigned norm1, unsigned norm2, std::string on_tag){
        auto dat = m_mesh.property_map<V_ind, Vector>(on_tag);
        if (!dat.second) throw std::runtime_error("Object \"" + name + "\" doesn't have <V_ind, Vector> property \"" + on_tag + "\"");
        return residual<Vector>(dat.first, norm1, norm2);
    }

//    double Object3D::residual(unsigned norm1, unsigned norm2, std::string on_tag){
//        auto dat = m_mesh.property_map<V_ind, Vector>(on_tag);
//        if (!dat.second) throw std::runtime_error("Object \"" + name + "\" doesn't have <V_ind, Vector> property \"" + on_tag + "\"");
//        std::function<double(Mesh::Property_map<V_ind, Vector>&)> nrm;
//        if (norm1 != norm2) {
//            std::function<double(Mesh::Property_map < V_ind, Vector > &)> &nrm1 = nrm;
//            std::function<double(Vector &)> nrm2;
//            switch (norm2) {
//                case 0:
//                    nrm2 = [](Vector &v) {
//                        double res = 0;
//                        for (int i = 0; i < v.dimension(); ++i)
//                            res = std::max(res, static_cast<double>(fabs(v[i])));
//                        return res;
//                    };
//                    break;
//                case 1:
//                    nrm2 = [](Vector &v) {
//                        double res = 0;
//                        for (int i = 0; i < v.dimension(); ++i)
//                            res += fabs(v[i]);
//                        return res;
//                    };
//                    break;
//                case 2:
//                    nrm2 = [](Vector &v) {
//                        double res = 0;
//                        for (int i = 0; i < v.dimension(); ++i) {
//                            res += v[i] * v[i];
//                        }
//                        return std::sqrt(res);
//                    };
//                    break;
//                default:
//                    nrm2 = [norm2](Vector &v) {
//                        double res = 0;
//                        for (int i = 0; i < v.dimension(); ++i)
//                            res += pow(fabs(v[i]), norm2);
//                        return pow(res, 1.0 / norm2);
//                    };
//            }
//            switch (norm1) {
//                case 0:
//                    nrm1 = [&nrm2](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v)
//                            res = std::max(res, fabs(nrm2(i)));
//                        return res;
//                    };
//                    break;
//                case 1:
//                    nrm1 = [&nrm2](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v)
//                            res += fabs(nrm2(i));
//                        return res;
//                    };
//                    break;
//                case 2:
//                    nrm1 = [&nrm2](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v) {
//                            double q = fabs(nrm2(i));
//                            res += q * q;
//                        }
//                        return sqrt(res);
//                    };
//                    break;
//                default:
//                    nrm1 = [norm1, &nrm2](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v) {
//                            double q = nrm2(i);
//                            res += pow(q, norm1);
//                        }
//                        return pow(res, 1.0 / norm1);
//                    };
//                    break;
//            }
//        } else {
//            switch (norm1) {
//                case 0:
//                    nrm = [](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v)
//                            for (int j = 0; j < i.dimension(); ++j)
//                                res = std::max(res, static_cast<double>(fabs(i[j])));
//                        return res;
//                    };
//                    break;
//                case 1:
//                    nrm = [](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v)
//                            for (int j = 0; j < i.dimension(); ++j)
//                                res += fabs(i[j]);
//                        return res;
//                    };
//                    break;
//                case 2:
//                    nrm = [](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v) {
//                            for (int j = 0; j < i.dimension(); ++j)
//                                res += i[j] * i[j];
//                        }
//                        return sqrt(res);
//                    };
//                    break;
//                default:
//                    nrm = [norm1](Mesh::Property_map <V_ind, Vector> &v) {
//                        double res = 0;
//                        for (auto &i: v) {
//                            for (int j = 0; j < i.dimension(); ++j)
//                                res += pow(fabs(i[j]), norm1);
//                        }
//                        return pow(res, 1.0 / norm1);
//                    };
//                    break;
//            }
//        }
//        return nrm(dat.first);
//    }

    bool check_flat_initial_template(Object3D& obj, double tol){
        if (obj.m_mesh.num_faces() <= 0) return true;
        auto n = [&obj](F_ind f){
            auto v = vert_around(obj.m_mesh, f);
            auto n = CGAL::cross_product(Vector(obj.m_x0[v[0]], obj.m_x0[v[1]]),
                                Vector(obj.m_x0[v[0]], obj.m_x0[v[2]]));
            return n/sqrt(n.squared_length());
        };
        auto n0 = n(*obj.m_mesh.faces().begin());
        for (auto f: obj.m_mesh.faces())
            if (fabs(fabs(n0 * n(f)) - 1) > tol) return false;
        return true;
    }

};
