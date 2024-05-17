//
// Created by alex on 17.11.2021.
//

#ifndef AORTIC_VALVE_CGALCOLLISSIONMANAGER_H
#define AORTIC_VALVE_CGALCOLLISSIONMANAGER_H

#include "../CollisionManagerInterface.h"

#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/squared_distance_3.h>

template <typename NT>
static inline NT Clamp(const NT& x, const NT& l, const NT& h)
{
    return (x < l ? l : x > h ? h : x);
}
template <class Vector3d, typename NT>
static inline void ProjectOrigin(const Vector3d& a,
                                 const Vector3d& b,
                                 Vector3d& prj,
                                 NT& sqd)
{
    const Vector3d d = b - a;
    auto squared_length = [](const Vector3d& d) { return d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; };
    const NT m2 = squared_length(d);
    if (m2 > DBL_EPSILON)
    {
        const NT t = Clamp<NT>(-a * d / m2, 0, 1);
        const Vector3d p = a + d * t;
        const NT l2 = squared_length(p);
        if (l2 < sqd)
        {
            prj = p;
            sqd = l2;
        }
    }
    else {
        const NT l2 = squared_length(a);
        if (l2 < sqd)
        {
            prj = a;
            sqd = l2;
        }
    }
}
template <class Vector3d, typename NT>
static inline void ProjectOrigin(const Vector3d& a,
                                 const Vector3d& b,
                                 const Vector3d& c,
                                 Vector3d& prj,
                                 NT& sqd)
{
    auto squared_length = [](const Vector3d& d) { return d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; };
    auto cross_product = [](const Vector3d& a, const Vector3d& b){
        return Vector3d(a[1]*b[2] - a[2]*b[1],a[2]*b[0] - a[0]*b[2],a[0]*b[1] - a[1]*b[0]);
    };
    const Vector3d q = cross_product(b - a, c - a);
    const double m2 = squared_length(q);
    if (m2 > DBL_EPSILON)
    {
        const Vector3d n = q / sqrt(m2);
        const double k = a * n;
        const double k2 = k * k;
        if (k2 < sqd)
        {
            const Vector3d p = n * k;
            if ((cross_product(a - p, b - p) * q > 0) &&
                (cross_product(b - p, c - p) * q > 0) &&
                (cross_product(c - p, a - p) * q > 0))
            {
                prj = p;
                sqd = k2;
            }
            else
            {
                ProjectOrigin(a, b, prj, sqd);
                ProjectOrigin(b, c, prj, sqd);
                ProjectOrigin(c, a, prj, sqd);
            }
        }
    }
    else
        ProjectOrigin(a, b, prj, sqd);
}
template <class Vector3d, typename NT>
static inline void ProjectOnTriangle(const Vector3d& a,
                                     const Vector3d& b,
                                     const Vector3d& c,
                                     const Vector3d& o,
                                     Vector3d& prj,
                                     NT& sqd)
{
    sqd = DBL_MAX;
    auto vprj = prj - o;
    ProjectOrigin(a - o, b - o, c - o, vprj, sqd);
    prj = o + vprj;
}

template<class Vector3d>
static inline Vector3d BaryCoord(const Vector3d& a,
                                  const Vector3d& b,
                                  const Vector3d& c,
                                  const Vector3d& p)
{
    auto cross_product = [](const auto& a, const auto& b){
        return Vector3d(a[1]*b[2] - a[2]*b[1],a[2]*b[0] - a[0]*b[2],a[0]*b[1] - a[1]*b[0]);
    };
    auto norm = [](const auto& a) { return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]); };
    const auto w0 = norm(cross_product(a - p, b - p)), w1 = norm(cross_product(b - p, c - p)), w2 = norm(cross_product(c - p, a - p));
    const auto isum = 1 / (w0 + w1 + w2);
    return (Vector3d(w1 * isum, w2 * isum, w0 * isum));
}
template <typename T, typename R>
static inline T BaryEval(const T& a,
                         const T& b,
                         const T& c,
                         const R& coord)
{
    return (a * coord[0] + b * coord[1] + c * coord[2]);
}

template<class NT_, class Handle_>
class Box_with_handle
        : public CGAL::Box_intersection_d::Box_d< NT_, 3> {
protected:
    Handle_ m_handle;
public:
    typedef CGAL::Box_intersection_d::Box_d< NT_, 3>    Base;
    typedef CGAL::Bbox_3                                Bbox_3;
    typedef NT_                                         NT;
    typedef Handle_                                     Handle;
    typedef Handle                                      ID;

    Box_with_handle() {}
    Box_with_handle( Handle h) : m_handle(h) {}
    Box_with_handle( bool complete, Handle h): Base(complete), m_handle(h) {}
    Box_with_handle(NT l[3], NT h[3], Handle n) : Base( l, h), m_handle(n) {}
    Box_with_handle( const Bbox_3& b, Handle h) : Base( b), m_handle(h) {}
    Handle handle() const { return m_handle; }
    ID  id() const { return m_handle; }
};

template <class Point, typename NT>
CGAL::Bbox_3 make_bbox(const Point &p, NT size) {
    NT eps = size / 2;
    return CGAL::Bbox_3(p[0] - eps, p[1] - eps, p[2] - eps, p[0] + eps, p[1] + eps, p[2] + eps);
}
template <class Point, typename NT, int N>
CGAL::Bbox_3 make_bbox(const std::array<Point, N>& pc, NT size){
    std::array<NT, 3> p[2] = {{pc[0][0], pc[0][1], pc[0][2]}, {pc[0][0], pc[0][1], pc[0][2]}};
    for (int i = 1; i < N; ++i)
        for (int j = 0; j < 3; ++j)
            if (p[0][j] > pc[i][j]) p[0][j] = pc[i][j];
            else if (p[1][j] < pc[i][j]) p[1][j] = pc[i][j];
    NT eps = size / 2;
    for (int j = 0; j < 3; ++j){
        p[0][j] -= eps;
        p[1][j] += eps;
    }
    return CGAL::Bbox_3(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]);
}
template <class Point, typename NT>
CGAL::Bbox_3 make_bbox(const Point &p1, const Point &p2, const Point &p3, NT size) {
    std::array<NT, 3> p[2] = {{p1[0], p1[1], p1[2]}, {p1[0], p1[1], p1[2]}};
    for (int j = 0; j < 3; ++j)
        if (p[0][j] > p2[j]) p[0][j] = p2[j];
        else if (p[1][j] < p2[j]) p[1][j] = p2[j];
    for (int j = 0; j < 3; ++j)
        if (p[0][j] > p3[j]) p[0][j] = p3[j];
        else if (p[1][j] < p3[j]) p[1][j] = p3[j];
    NT eps = size / 2;
    for (int j = 0; j < 3; ++j){
        p[0][j] -= eps;
        p[1][j] += eps;
    }
    return CGAL::Bbox_3(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]);
}
template <class Point, typename NT>
CGAL::Bbox_3 make_bbox(const Point &p1, const Point &p2, NT size) {
    std::array<NT, 3> p[2] = {{p1[0], p1[1], p1[2]}, {p2[0], p2[1], p2[2]}};
    for (int j = 0; j < 3; ++j) if (p[0][j] > p[1][j]) std::swap(p[0][j], p[1][j]);
    NT eps = size / 2;
    for (int j = 0; j < 3; ++j){
        p[0][j] -= eps;
        p[1][j] += eps;
    }
    return CGAL::Bbox_3(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]);
}

static CGAL::Bbox_3 make_bbox(World3d::Object3D& obj, World3d::V_ind v, double size){ return make_bbox(obj.m_x[v], size); }
static CGAL::Bbox_3 make_bbox(World3d::Object3D& obj, World3d::E_ind e, double size){
    auto v = World3d::vert_around(obj.m_mesh, e);
    return make_bbox(obj.m_x[v[0]], obj.m_x[v[1]], size);
}
static CGAL::Bbox_3 make_bbox(World3d::Object3D& obj, World3d::F_ind f, double size){
    auto v = World3d::vert_around(obj.m_mesh, f);
    return make_bbox(obj.m_x[v[0]], obj.m_x[v[1]], obj.m_x[v[2]], size);
}
static std::vector<CGAL::Bbox_3> make_bbox(World3d::Object3D& obj, unsigned odf_type, double size){
    unsigned sz = ((odf_type & 1) != 0) * obj.m_mesh.num_vertices() + ((odf_type & 2) != 0) * obj.m_mesh.num_edges() + ((odf_type & 4) != 0) * obj.m_mesh.num_faces();
    std::vector<CGAL::Bbox_3> res; res.reserve(sz);
    if (odf_type & 1)
        for (auto v: obj.m_mesh.vertices())
            res.push_back(make_bbox(obj, v, size));
    if (odf_type & 2)
        for (auto e: obj.m_mesh.edges())
            res.push_back(make_bbox(obj, e, size));
    if (odf_type & 4)
        for (auto f: obj.m_mesh.faces())
            res.push_back(make_bbox(obj, f, size));
    return res;
}
static std::vector<CGAL::Bbox_3> make_bbox_on_vertices(World3d::Object3D& obj, double size){
    std::vector<CGAL::Bbox_3> res; res.reserve(obj.m_mesh.num_vertices());
    for (auto v: obj.m_mesh.vertices())
        res.push_back(make_bbox(obj, v, size));
    return res;
}
static std::vector<CGAL::Bbox_3> make_bbox_on_edges(World3d::Object3D& obj, double size){
    std::vector<CGAL::Bbox_3> res; res.reserve(obj.m_mesh.num_edges());
    for (auto e: obj.m_mesh.edges())
        res.push_back(make_bbox(obj, e, size));
    return res;
}
static std::vector<CGAL::Bbox_3> make_bbox_on_faces(World3d::Object3D& obj, double size){
    std::vector<CGAL::Bbox_3> res; res.reserve(obj.m_mesh.num_faces());
    for (auto f: obj.m_mesh.faces())
        res.push_back(make_bbox(obj, f, size));
    return res;
}

namespace World3d {
    //cgal collision manager vertex vs face with saving of all collision pairs
    class CGALColManVFfull: public CollisionManagerBase {
    public:
        using VBBox = Box_with_handle<double, std::size_t/*V_ind*/>;
        using FBBox = Box_with_handle<double, std::size_t/*F_ind*/>;
        struct ObjColInfo;
        struct ContactData{
            V_ind v; //on v
            F_ind f; //from f
            Point proj;
            double dist2;
            Vector w, nr;
            double margin;
            ObjColInfo* from;
        };
        struct ObjColInfo{
            Object3D* obj;
            int status;
            double margin;
            std::vector<VBBox> vb;
            std::vector<FBBox> fb;
            std::vector<ContactData> cd;
        };
    private:
        std::map<ObjectID, ObjColInfo> m_shape;

        static int default_solve_collision_algo(std::map<ObjectID, ObjColInfo>& shape);
        std::function<int(std::map<ObjectID, ObjColInfo>&)> m_solve_collission = [](std::map<ObjectID, ObjColInfo>& shape){
            return default_solve_collision_algo(shape);
        };
    public:
        void addObj(Object3D& obj, ObjectID id, bool is_dynamic, double margin) override{
            ObjColInfo info;
            info.obj = &obj;
            info.status = is_dynamic;
            info.margin = 2*margin;
            m_shape.insert({id, std::move(info)});
        }
        bool set_margin(ObjectID id, double margin);
        void removeObj(ObjectID id) override { m_shape.erase(id); }
        void findCollisions() override;

        void solveCollisions() override{
            m_solve_collission(m_shape);
            for (auto & i1 : m_shape) i1.second.cd.resize(0);
        }
    };

    //cgal collision manager vertex vs face_center with saving of all collision pairs
    class CGALColManVFcfull: public CollisionManagerBase {
    public:
        using VBBox = Box_with_handle<double, std::size_t/*V_ind*/>;
        using FBBox = Box_with_handle<double, std::size_t/*F_ind*/>;
        struct ObjColInfo;
        struct ContactData{
            V_ind v; //on v
            F_ind f; //from f
            Vector n_cur;
            double dpi_cur;
            ObjColInfo* from;
        };
        struct ObjColInfo{
            Object3D* obj;
            int status;
            double d_e;
            std::vector<VBBox> vb;
            std::vector<FBBox> fb;
            std::vector<ContactData> cd;
            std::vector<Vector> work_array;
        };
    private:
        double m_de_scale = 1.0;
        double m_d_pi = 0.7;
        std::array<double, 2> m_k = {1e4, 1e4};
        std::map<ObjectID, ObjColInfo> m_shape;

        static int default_solve_collision_algo(std::map<ObjectID, ObjColInfo>& shape, const std::array<double, 2>& m_k){
            for (auto& ob: shape){
                Object3D* to = ob.second.obj;
                if (ob.second.cd.empty()) continue;
                ob.second.work_array.resize(to->m_mesh.num_vertices());
                std::fill(ob.second.work_array.begin(), ob.second.work_array.end(), CGAL::NULL_VECTOR);
                auto& temp_F = ob.second.work_array;
                for (auto col: ob.second.cd){
                    Object3D* from = col.from->obj;
                    if (!to->is_movable(col.v)) continue;
                    double sqF = to->m_F[col.v].squared_length();
                    double F = sqrt(sqF);
                    double dpi = col.dpi_cur;
                    Vector F_cont = ((dpi > 0) ? F * exp(-m_k[0]*dpi / F) : (F - m_k[1]*dpi)) * col.n_cur;
                    temp_F[col.v] += F_cont;
                }
                for (auto v: to->m_mesh.vertices())
                    to->m_F[v] += temp_F[v];
            }

            for (auto& ob: shape) {
                Object3D *to = ob.second.obj;
                if (ob.second.cd.empty()) continue;
                to->update_next_x();
            }

            return 0;
        }
        std::function<int(std::map<ObjectID, ObjColInfo>&, const std::array<double, 2>&)> m_solve_collission = [](std::map<ObjectID, ObjColInfo>& shape, const std::array<double, 2> m_k){
            return default_solve_collision_algo(shape, m_k);
        };
    public:
        void addObj(Object3D& obj, ObjectID id, bool is_dynamic, double margin) override{
            ObjColInfo info;
            info.obj = &obj;
            info.status = is_dynamic;
            if (margin > 0) info.d_e = margin;
            else {
                info.d_e = 0;
                for (auto f: obj.m_mesh.faces()){
                    auto v = vert_around(obj.m_mesh, f);
                    Vector c = CGAL::NULL_VECTOR;
                    for (int k = 0; k < 3; ++k) c += obj.m_x0[v[k]] - CGAL::ORIGIN;
                    c /= 3;
                    double de_f2 = (obj.m_x0[v[0]] - (CGAL::ORIGIN + c)).squared_length();
                    for (int k = 1; k < 3; ++k) de_f2 = std::max(de_f2, (obj.m_x0[v[k]] - (CGAL::ORIGIN + c)).squared_length());
                    info.d_e += sqrt(de_f2);
                }
                info.d_e /= obj.m_mesh.num_faces();
            }
            m_shape.insert({id, std::move(info)});
        }
        CGALColManVFcfull& set_de_scale(double c) { m_de_scale = c; return *this; }
        CGALColManVFcfull& set_d_pi_treshhold(double d_pi) { m_d_pi = d_pi; return *this; }
        CGALColManVFcfull& set_params(double k1, double k2) { m_k[0] = k1, m_k[1] = k2; return *this;}
        void removeObj(ObjectID id) override { m_shape.erase(id); }
        void findCollisions() override{
            double gmrg = 0;
            {
                int nf = 0;
                for (auto &i: m_shape)
                    gmrg += i.second.d_e * i.second.obj->m_mesh.num_faces(), nf += i.second.obj->m_mesh.num_faces();
                if (nf) gmrg *= m_de_scale / nf;
            }
            for (auto& i: m_shape){
                auto& obj = *i.second.obj;
                auto getP = [&obj](V_ind v) { return obj.m_x[v] - CGAL::ORIGIN; };
                auto mrg = gmrg;
                auto& fb = i.second.fb;
                if (fb.empty()) {
                    fb.resize(obj.m_mesh.num_faces());
                    int ind = 0;
                    for (auto f: obj.m_mesh.faces()) {
                        auto v = vert_around(obj.m_mesh, f);
                        fb[ind++] = FBBox(make_bbox((getP(v[0])+getP(v[1])+getP(v[2])) / 3, mrg), f);
                    }
                } else {
                    for (auto& fd: fb) {
                        auto f = F_ind(fd.handle());
                        auto v = vert_around(obj.m_mesh, f);
                        fd = FBBox(make_bbox((getP(v[0])+getP(v[1])+getP(v[2])) / 3, mrg), f);
                    }
                }
                if ((i.second.status & 1) != 0){
                    auto& vb = i.second.vb;
                    if (vb.empty()){
                        vb.resize(obj.m_mesh.num_vertices());
                        int ind = 0;
                        for (auto v: obj.m_mesh.vertices()) {
                            vb[ind++] = VBBox(make_bbox(getP(v), mrg), v);
                        }
                    } else {
                        for (auto& vd: vb) {
                            auto v = V_ind(vd.handle());
                            vd = VBBox(make_bbox(getP(v), mrg), v);
                        }
                    }
                }
            }
            for (auto i1 = m_shape.begin(); i1 != m_shape.end(); ++i1) { //from
                auto& fb = i1->second.fb;
                for (auto i2 = m_shape.begin(); i2 != m_shape.end(); ++i2) { //to
                    if (i1 == i2) continue;
                    if ((i2->second.status & 1) == 0) continue;
                    auto& vb = i1->second.vb;
                    //double mrg = (i1->second.margin + i2->second.margin) / 2;
                    auto check_collision = [ob1 = i1->second.obj, ob2 = i2->second.obj, dpi = m_d_pi, &contacts = i2->second.cd, &it = i1](const VBBox &vb, const FBBox &fb) {
                        auto fv = World3d::vert_around(ob1->m_mesh, F_ind(fb.handle()));
                        auto vv = V_ind(vb.handle());
                        std::array<Point, 3> tri_next = {ob1->m_next_x[fv[0]], ob1->m_next_x[fv[1]], ob1->m_next_x[fv[2]]};
                        Vector n_next = CGAL::cross_product(Vector(tri_next[0], tri_next[1]), Vector(tri_next[0], tri_next[2]));
                        n_next /= sqrt(n_next.squared_length());
                        Point p_next = ob2->m_next_x[vv];
                        double ldpi_next = (p_next - tri_next[0]) * n_next;
                        if (ldpi_next < dpi){ //collision occur
                            std::array<Point, 3> tri_cur = {ob1->m_x[fv[0]], ob1->m_x[fv[1]], ob1->m_x[fv[2]]};
                            Vector n_cur = CGAL::cross_product(Vector(tri_cur[0], tri_cur[1]), Vector(tri_cur[0], tri_cur[2]));
                            n_cur /= sqrt(n_cur.squared_length());
                            Point p_cur = ob2->m_x[vv];
                            double ldpi_cur = (p_cur - tri_cur[0]) * n_cur;
                            ContactData cd;
                            cd.v = vv;
                            cd.f = F_ind(fb.handle());
                            cd.dpi_cur = ldpi_cur;
                            cd.n_cur = n_cur;
                            cd.from = &it->second;
                            contacts.push_back(cd);
                        }
                    };
                    CGAL::box_intersection_d(vb.begin(), vb.end(), fb.begin(), fb.end(), check_collision);
                }
            }
        }

        void solveCollisions() override{
            m_solve_collission(m_shape, m_k);
            for (auto & i1 : m_shape) i1.second.cd.resize(0);
        }
    };
}


#endif //AORTIC_VALVE_CGALCOLLISSIONMANAGER_H
