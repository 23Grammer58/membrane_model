//
// Created by alex on 29.11.2020.
//

#include "Messurer.h"
#include <gsl/gsl_multimin.h>

template <typename T>
static inline T Clamp(const T& x, const T& l, const T& h)
{
    return (x < l ? l : x > h ? h : x);
}
static inline void ProjectOrigin(const World3d::Vector& a,
                                 const World3d::Vector& b,
                                 World3d::Vector& prj,
                                 double& sqd)
{
    const World3d::Vector d = b - a;
    const double m2 = d.squared_length();
    if (m2 > DBL_EPSILON)
    {
        const double t = Clamp<double>(-a * d / m2, 0, 1);
        const World3d::Vector p = a + d * t;
        const double l2 = p.squared_length();
        if (l2 < sqd)
        {
            prj = p;
            sqd = l2;
        }
    }
    else {
        const double l2 = a.squared_length();
        if (l2 < sqd)
        {
            prj = a;
            sqd = l2;
        }
    }
}

static inline void ProjectOrigin(const World3d::Vector& a,
                                 const World3d::Vector& b,
                                 const World3d::Vector& c,
                                 World3d::Vector& prj,
                                 double& sqd)
{
    using namespace World3d;
    const Vector& q = CGAL::cross_product(b - a, c - a);
    const double m2 = q.squared_length();
    if (m2 > DBL_EPSILON)
    {
        const Vector n = q / sqrt(m2);
        const double k = a * n;
        const double k2 = k * k;
        if (k2 < sqd)
        {
            const Vector p = n * k;
            if ((CGAL::cross_product(a - p, b - p) * q > 0) &&
                (CGAL::cross_product(b - p, c - p) * q > 0) &&
                (CGAL::cross_product(c - p, a - p) * q > 0))
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

static inline void ProjectOnTriangle(const World3d::Point& a,
                                     const World3d::Point& b,
                                     const World3d::Point& c,
                                     const World3d::Point& o,
                                     World3d::Point& prj,
                                     double& sqd)
{
    sqd = DBL_MAX;
    World3d::Vector vprj = prj - CGAL::ORIGIN;
    ProjectOrigin(a - o, b - o, c - o, vprj, sqd);
    prj = o + vprj;
}

using namespace World3d;
struct CBR{
    using ObjectShape = Messurer::ObjectShape;
    using CollisionMap = Messurer::CollisionMap;
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
    struct CollisionInternalInfo{
        typedef Messurer::CollisionPairData CPD;
        //если когда-то через это будет реализовываваться динамический вариант, то надо будет изменить типы для коробки

        Messurer::CollisionInfo data;
        typedef Box_with_handle<double, std::size_t/*V_ind*/> BoxFrom;
        typedef Box_with_handle<double, std::size_t/*F_ind*/> BoxTo;
        std::array<std::vector<BoxFrom>, 3> bfrom;
        std::array<std::vector<BoxTo>, 3> bto;

        std::vector<CPD>& getCollisionInfo(int from, int to) { return data.getCollisionInfo(from, to); };
        std::vector<CPD>& operator()(int from, int to){ return getCollisionInfo(from, to); }

        void setBoxes(const std::array<ObjectShape, 3>& shapes, const std::array<double, 3>& mrg_from, const std::array<double, 3>& mrg_to)
        {
            for (int i = 0; i < 3; ++i){
                bfrom[i].resize(shapes[i].mesh.points.size());
                bto[i].resize(shapes[i].obj->m_mesh.num_faces());
                for (int v = 0; v < bfrom[i].size(); ++v)
                    bfrom[i][v] = BoxFrom(shapes[i].mesh.bbox(V_ind(v), mrg_from[i]), V_ind(v));
                for (int f = 0; f < bto[i].size(); ++f) {
                    auto v = World3d::vert_around(shapes[i].obj->m_mesh, F_ind(f));
                    auto& p = shapes[i].obj->m_x;
                    bto[i][f] = BoxTo(Messurer::TriangleSoup::bbox(p[v[0]], p[v[1]], p[v[2]], mrg_to[i]), F_ind(f));
                }
            }
        }
        void setBoxes(const std::array<ObjectShape, 3>& shapes, const std::array<double, 3>& mrg){ setBoxes(shapes, mrg, mrg); }
        void setBoxes(const std::array<ObjectShape, 3>& shapes, double mrg){ setBoxes(shapes, {mrg, mrg, mrg}); }
    };

    using CollisionInfo = CollisionInternalInfo;
    using SpecialCase = std::function<bool(const ObjectShape& from, const ObjectShape& to, std::vector<CollisionInfo::CPD>& info)>;

    const std::array<Object3D*, 3>& m_cusps;
    std::array<ObjectShape, 3> m_colission_shapes;
    double m_margin;
    CollisionInfo m_colission_info;
    std::array<CollisionMap, 3> m_colMap;

    std::array<int, 3> prev_npnts = {0, 0, 0};
    std::vector<CollisionInfo::CPD> colPairtmp;

    SpecialCase m_specCase;

    CBR(const std::array<Object3D*, 3>& cusps, double margin): m_cusps{cusps}, m_margin{margin} {
        for (int i = 0; i < 3; ++i)
            m_colission_shapes[i] = ObjectShape(m_cusps[i]);
    }
    void setSpecialCase(SpecialCase cs) { m_specCase = std::move(cs); }
    void prepareBoxes(){
        m_colission_info.setBoxes(m_colission_shapes, m_margin);
    }
    void updateBoxes(){
        updateBoxes(std::array<double, 3>{m_margin, m_margin, m_margin});
    }
    void updateBoxes(const std::array<double, 3>& mrg_from){
        for (int i = 0; i < 3; ++i){
            auto& bfrom = m_colission_info.bfrom[i];
            auto& shape = m_colission_shapes[i];
            bfrom.resize(shape.mesh.points.size() - prev_npnts[i]);
            for (int v = 0; v < bfrom.size(); ++v)
                bfrom[v] = CollisionInfo::BoxFrom(shape.mesh.bbox(V_ind(v+prev_npnts[i]), mrg_from[i]), V_ind(v+prev_npnts[i]));
        }
    }
    void computeCollisionInfo(){
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j){
                if (i == j) continue;
                auto& bfrom = m_colission_info.bfrom[i];
                auto& bto = m_colission_info.bto[j];
                const auto& from = m_colission_shapes[i];
                const auto& to = m_colission_shapes[j];
                colPairtmp.resize(0);

                CGAL::box_intersection_d(bfrom.begin(), bfrom.end(), bto.begin(), bto.end(),
                                         [&from, &to, &info = colPairtmp, mrg2 = m_margin * m_margin](const CollisionInfo::BoxFrom& b1, const CollisionInfo::BoxTo& b2){
                                             auto v = World3d::vert_around(to.obj->m_mesh, F_ind(b2.handle()));
                                             auto& triP = to.obj->m_x;
                                             Point p = CGAL::ORIGIN; double sqrd = DBL_MAX;
                                             ProjectOnTriangle(triP[v[0]], triP[v[1]], triP[v[2]], from.mesh.points[b1.handle()], p, sqrd);
                                             if (sqrd < mrg2) {
                                                 CollisionInfo::CPD c;
                                                 c.v = V_ind(b1.handle());
                                                 c.f = F_ind(b2.handle());
                                                 c.fp = p;
                                                 c.dist2 = sqrd;
                                                 info.push_back(c);
                                             }
                                         });
                if (m_specCase) m_specCase(from, to, colPairtmp);
                std::sort(colPairtmp.begin(), colPairtmp.end(), [](const auto& a, const auto& b){ return (a.v != b.v) ? a.v < b.v : a.dist2 < b.dist2; });
                {
                    int jj = colPairtmp.size() > 0;
                    for (int l = 1; l < colPairtmp.size(); ++l) {
                        if (colPairtmp[l].v != colPairtmp[jj-1].v)
                            colPairtmp[jj++] = colPairtmp[l];
                    }
                    colPairtmp.resize(jj);
                }
                auto& info = m_colission_info(i, j);
                for (auto& dat: colPairtmp) info.push_back(dat);
            }
    }

    void makeCollisionMap() {
        for (int i = 0; i < 3; ++i){
            auto& lbl = m_colMap[i].face_lbl;
            auto& remap = m_colMap[i].remap;
            auto& invmap = m_colMap[i].invmap;
            auto& shape = m_colission_shapes[i];
            lbl.resize(shape.mesh.tria.size());
            invmap.resize(shape.mesh.tria.size());
            std::fill(lbl.begin(), lbl.end(), 0);
            std::vector<bool> vc1(shape.mesh.points.size(), false), vc2(shape.mesh.points.size(), false);
            for (auto vc: m_colission_info(i, (i+1)%3)){ vc1[vc.v] = true; }
            for (auto vc: m_colission_info(i, (i+2)%3)){ vc2[vc.v] = true; }
            {
                int j = 0;
                for (auto& triit: shape.mesh.tria){
                    remap[triit.first] = j;
                    invmap[j] = triit.first;
                    auto& tri = triit.second;
                    for (int k = 0; k < 3; ++k){
                        lbl[j] |= ((1 << (i+1)%3)*vc1[tri[k]] | (1 << (i+2)%3)*vc2[tri[k]]) << 3*k;
                    }
                    ++j;
                }
            }
        }
    }

    void preparePartialyCollideTriags(int i, std::vector<F_ind>& partialy){
        partialy.resize(0);
        auto& lbl = m_colMap[i].face_lbl;
        auto& shape = m_colission_shapes[i];
        for (int j = 0; j < lbl.size(); ++j) {
            if (lbl[j]) {
                bool flag = true;
                std::array<bool, 3> belong = {false, false, false};
                for (int k = 0; k < 3; ++k) {
                    if (!(lbl[j] & (7 << 3 * k)))
                        flag = false;
                    else belong[k] = true;
                }
                if (!flag) {
                    auto f = m_colMap[i].invmap[j];
                    partialy.push_back(m_colMap[i].invmap[j]);
                    int r = static_cast<int>(belong[0]) + belong[1] + belong[2];
                    int opposite = -1;
                    if (r == 1) {
                        for (int l = 0; l < 3; ++l)
                            if (belong[l])
                                opposite = l;
                    } else {
                        for (int l = 0; l < 3; ++l)
                            if (!belong[l])
                                opposite = l;
                    } 
                    auto vv = World3d::vert_around(shape.obj->m_mesh, f);
                    auto ee = World3d::edge_around(shape.obj->m_mesh, f);
                    E_ind e;
                    for (auto ei: ee) {
                        auto v = World3d::vert_around(shape.obj->m_mesh, ei);
                        if (v[0] != vv[opposite] && v[1] != vv[opposite])
                            e = ei;
                    }
                    auto ff = World3d::face_around(shape.obj->m_mesh, e);
                    std::pair<F_ind, bool> addf;
                    if (ff[0].first != f) addf = ff[0];
                    else addf = ff[1];
                    if (addf.second) partialy.push_back(addf.first);
                } else { //is face on central coaptation line?
                    char mask = (1 << ((i+1)%3)) | (1 << ((i+2)%3));
                    bool partialy_on_central_line = !((lbl[j] >> 3*0) & (lbl[j] >> 3*1) & (lbl[j] >> 3*2) & mask);
                    if (partialy_on_central_line)
                        partialy.push_back(m_colMap[i].invmap[j]);
                }
            }
        }
    }

    void computePartialyCollideTriags(int i, std::vector<F_ind>& partialy){
        partialy.resize(0);
        auto& lbl = m_colMap[i].face_lbl;
        auto& shape = m_colission_shapes[i];
        for (int j = 0; j < lbl.size(); ++j) {
            if (lbl[j]) {
                char mask = (1 << ((i+1)%3)) | (1 << ((i+2)%3));
                if (!((lbl[j] >> 3*0) & (lbl[j] >> 3*1) & (lbl[j] >> 3*2) & mask))
                    partialy.push_back(m_colMap[i].invmap[j]);
            }
        }
    }

    void refinePartiallyCollidingElems(int i, std::vector<F_ind>& partialy) {
        if (partialy.empty()) return;
        using TriangleSoup = Messurer::TriangleSoup;
        auto& lbl = m_colMap[i].face_lbl;
        auto& shape = m_colission_shapes[i];
        prev_npnts[i] = shape.mesh.points.size();
        //divide partially colliding triangles into 4 parts
        std::vector<std::pair<V_ind, Point>> newPnts;
        newPnts.reserve(partialy.size()*3);
        std::vector<std::pair<TriangleSoup::Triangle, ObjectShape::ExtraFace>> newTria;
        newTria.reserve(partialy.size()*4);

        for (int fi = 0; fi < partialy.size(); ++fi){
            const static std::array<Vector, 6> w = {Vector{1, 0, 0},
                                                    Vector{0, 1, 0},
                                                    Vector{0, 0, 1},
                                                    Vector{0, 0.5, 0.5},
                                                    Vector{0.5, 0, 0.5},
                                                    Vector{0.5, 0.5, 0}};
            const static std::array<std::array<int, 3>, 4> tri_nums = {std::array<int, 3>{0, 5, 4},
                                                                       {1, 3, 5},
                                                                       {2, 4, 3},
                                                                       {3, 4, 5}};
            auto f = partialy[fi];
            F_ind mf = f;
            auto v = shape.mesh.tria[f];
            auto off =  shape.mesh.points.size() + newPnts.size();
            auto& P = shape.mesh.points;
            newPnts.emplace_back(V_ind(off + 0), P[v[1]] + (P[v[2]] - P[v[1]]) / 2);
            newPnts.emplace_back(V_ind(off + 1), P[v[0]] + (P[v[2]] - P[v[0]]) / 2);
            newPnts.emplace_back(V_ind(off + 2), P[v[0]] + (P[v[1]] - P[v[0]]) / 2);
            std::array<V_ind, 6> vv{v[0], v[1], v[2], V_ind(off + 0), V_ind(off + 1), V_ind(off + 2)};
            std::array<Vector, 6> ww = w;
            //if we divide virtual triangle we need save real baricentric coords and reference to real triangle
            if (f >= shape.obj->m_mesh.num_faces()){
                mf = shape.extra_Finfo[f].f;
                auto wwP = shape.extra_Finfo[f].w;
                ww = {wwP[0], wwP[1], wwP[2], (wwP[1] + wwP[2])/2, (wwP[0] + wwP[2])/2, (wwP[0] + wwP[1])/2};
            }
            for (int ip = 0; ip < 4; ++ip){
                TriangleSoup::Triangle tr = {vv[tri_nums[ip][0]], vv[tri_nums[ip][1]], vv[tri_nums[ip][2]]};
                ObjectShape::ExtraFace ef;
                ef.f = mf;
                ef.w = {ww[tri_nums[ip][0]], ww[tri_nums[ip][1]], ww[tri_nums[ip][2]]};
                newTria.emplace_back(std::move(tr), std::move(ef));
            }
        }

        //remove repeated new Points
        std::sort(newPnts.begin(), newPnts.end(), [](const std::pair<V_ind, Point>& p1d, const std::pair<V_ind, Point>& p2d){
            const Point &p1 = p1d.second, &p2 = p2d.second;
            for (int i = 0; i < 3; ++i) {
                if (fabs(p1[i] - p2[i]) < 3 * DBL_EPSILON) continue;
                else return p1[i] < p2[i];
            }
            return p1d.first < p2d.first;
        });
        std::map<V_ind, V_ind> Pntmap;
        {
            auto norm = [](auto p) {
                std::array<double, 3> a;
                for (int i = 0; i < 3; ++i) a[i] = fabs(p[i]);
                return (a[0] < a[1]) ? (a[1] < a[2] ? a[2] : a[1]) : (a[0] < a[2] ? a[2] : a[0]);
            };
            int j = 0;
            Pntmap.insert({newPnts[0].first, V_ind(shape.mesh.points.size() + j)});
            Point p = newPnts[0].second;
            for (int i = 1; i < newPnts.size(); ++i) {
                if (norm(p - newPnts[i].second) < (norm(p) + norm(newPnts[i].second) + 1e-20) * DBL_EPSILON){
                    Pntmap.insert({newPnts[i].first, V_ind(shape.mesh.points.size() + j)});
                } else {
                    p = newPnts[i].second;
                    ++j;
                    Pntmap.insert({newPnts[i].first, V_ind(shape.mesh.points.size() + j)});
                    newPnts[j] = newPnts[i];
                }
            }
            newPnts.resize(j+1);
        }

        for (auto& tri: newTria)
            for (int l = 0; l < 3; ++l) {
                auto it = Pntmap.find(tri.first[l]);
                if (it != Pntmap.end())
                    tri.first[l] = it->second;
            }

        //сохраняем новые точки
        shape.mesh.points.reserve(shape.mesh.points.size() + newPnts.size());
        for (int ip = 0; ip < newPnts.size(); ++ip)
            shape.mesh.points.push_back(newPnts[ip].second);
        //удаление partially colliding triangles из shape.tria
        F_ind st = shape.mesh.tria.begin()->first;
        for (auto& f: shape.mesh.tria)
            st = std::max(st, f.first);
        for (auto& f: partialy) {
            shape.mesh.tria.erase(f);
            shape.extra_Finfo.erase(f);
        }
        for (auto& tri: newTria){
            ++st;
            shape.mesh.tria[st] = tri.first;
            shape.extra_Finfo[st] = tri.second;
        }

        m_colMap[i].clear();
    }

    static void refineTo(std::array<Object3D*, 3>& cusps, double margin, int level,
                         std::array<ObjectShape, 3>& colission_shapes, Messurer::CollisionInfo& colission_info, std::array<CollisionMap, 3>& colMap,
                         CBR::SpecialCase spc = nullptr){
        if (level <= 0) return;
        if (margin <= 0) return;
        if (level > 10) {
            std::cout << "Warning: level should be less or equal than 10, level = 10 will be used";
            level = 10;
        }
        CBR r(cusps, margin);
        if (spc) r.setSpecialCase(spc);
        r.prepareBoxes();
        r.computeCollisionInfo();
        r.makeCollisionMap();
        std::vector<F_ind> tmp;
        for (int i = 0; i < 3; ++i) {
            r.preparePartialyCollideTriags(i, tmp);
            r.refinePartiallyCollidingElems(i, tmp);
        }
        r.updateBoxes();
        r.computeCollisionInfo();
        r.makeCollisionMap();
        for (int l = 1; l < level; ++l){
            for (int i = 0; i < 3; ++i) {
                r.computePartialyCollideTriags(i, tmp);
                r.refinePartiallyCollidingElems(i, tmp);
            }
            r.updateBoxes();
            r.computeCollisionInfo();
            r.makeCollisionMap();
        }
        colission_shapes = std::move(r.m_colission_shapes);
        colission_info = std::move(r.m_colission_info.data);
        colMap = std::move(r.m_colMap);
    }
};

bool Messurer::saveObjectShape(const Messurer::ObjectShape &obj, std::string fname, std::string tag, std::function<bool (F_ind)> filter) {
    if (fname.substr(fname.size() - 4) != ".stl" && fname.substr(fname.size() - 4) != ".STL")
        fname += ".stl";

    std::ofstream ob(fname, std::ios::binary | std::ios::trunc);
    if (!ob)
        return false;
    auto m_x_it = obj.obj->m_mesh.property_map<V_ind, Point>(tag);
    if (!m_x_it.second)
        throw std::runtime_error("obj doesn't have tag = \"" + tag + "\"");
    auto& m_x = m_x_it.first;
    std::array<uint8_t, 80> zero = {0};
    ob.write(reinterpret_cast<const char *> ((zero.data())), 80);
    uint32_t nt = 0;
    for (auto &fi: obj.mesh.tria) nt += filter(fi.first);
    ob.write(reinterpret_cast<const char *> ((&nt)), sizeof(uint32_t));
    std::array<float, 3> v[4];
    for (auto &fi: obj.mesh.tria) {
        auto f = fi.first;
        if (!filter(f)) continue;
        auto evalP = [](const Vector& w, const Point& p1, const Point& p2, const Point& p3) -> Point{
            return Point(w[0]*p1.x()+w[1]*p2.x() + w[2]*p3.x(), w[0]*p1.y()+w[1]*p2.y() + w[2]*p3.y(), w[0]*p1.z()+w[1]*p2.z() + w[2]*p3.z());
        };
        std::array<Vector, 3> w = {Vector{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        auto fit = obj.extra_Finfo.find(fi.first);
        if (fit != obj.extra_Finfo.end()) {
            w = fit->second.w;
            f = fit->second.f;
        }
        auto vv = World3d::vert_around(obj.obj->m_mesh, f);
        Point PPT[3] = {m_x[vv[0]], m_x[vv[1]], m_x[vv[2]]};
        Point PP[3] = {evalP(w[0], PPT[0], PPT[1], PPT[2]), evalP(w[1], PPT[0], PPT[1], PPT[2]), evalP(w[2], PPT[0], PPT[1], PPT[2])};
        Vector n = CGAL::cross_product(PP[1] - PP[0], PP[2] - PP[0]);
        n /= sqrt(n.squared_length());
        for (int i = 0; i < 3; ++i)
            v[0][i] = n[i];
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i)
                v[j + 1][i] = PP[j][i];
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

void Messurer::computeCollidingBnd(int level, std::function<bool(const ObjectShape& from, const ObjectShape& to, std::vector<CollisionInfo::CPD>& info)> spec) {
    CBR::refineTo(m_cusps, m_margin, level, m_colission_shapes, m_colission_info, m_colMap, spec);
}

bool Messurer::checkCollisionShapes() {
    for (int i = 0; i < 3; ++i){
        auto& shape = m_colission_shapes[i];
        auto max_P_id =  shape.mesh.points.size();
        for (auto& triit: shape.mesh.tria){
            auto& tri = triit.second;
            for (int k = 0; k < 3; ++k){
                if (tri[k] >= max_P_id){
                    return false;
                }
            }
        }
    }
    return true;
}

std::vector<Messurer::CollisionPairData> &Messurer::CollisionInfo::getCollisionInfo(int from, int to) {
    assert(from != to && "from cusp number should be not equal to cusp number");
    int to_ind = (to - from + 3) % 3 - 1;
    return data[from][to_ind];
}

Messurer::TriangleSoup Messurer::TriangleSoup::makeTriangleSoup(const Object3D *obj) {
    TriangleSoup ts;
    ts.points.reserve(obj->m_mesh.num_vertices());
    for (auto f: obj->m_mesh.faces())
        ts.tria.insert({f, World3d::vert_around(obj->m_mesh, f)});
    for (auto v: obj->m_mesh.vertices())
        ts.points.push_back(obj->m_x[v]);
    return ts;
}

CGAL::Bbox_3 Messurer::TriangleSoup::bbox(Messurer::V_ind v, double size) const {
    const Point& p = points[v];
    double eps = size / 2;
    return CGAL::Bbox_3(p.x() - eps, p.y() - eps, p.z() - eps, p.x() + eps, p.y() + eps, p.z() + eps);
}

CGAL::Bbox_3 Messurer::TriangleSoup::bbox(const Point &p1, const Point &p2, const Point &p3, double size) {
    std::array<double, 3> p[2] = {{p1.x(), p1.y(), p1.z()}, {p2.x(), p2.y(), p2.z()}};
    for (int j = 0; j < 3; ++j) {
        p[0][j] = (p[0][j] > p2[j]) ? p2[j] : p[0][j];
        p[1][j] = (p[1][j] < p2[j]) ? p2[j] : p[1][j];
    }
    for (int j = 0; j < 3; ++j) {
        p[0][j] = (p[0][j] > p3[j]) ? p3[j] : p[0][j];
        p[1][j] = (p[1][j] < p3[j]) ? p3[j] : p[1][j];
    }
    double eps = size / 2;
    for (int j = 0; j < 3; ++j){
        p[0][j] -= eps;
        p[1][j] += eps;
    }
    return CGAL::Bbox_3(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]);
}

CGAL::Bbox_3 Messurer::TriangleSoup::bbox(const Point &p, double size) {
    double eps = size / 2;
    return CGAL::Bbox_3(p.x() - eps, p.y() - eps, p.z() - eps, p.x() + eps, p.y() + eps, p.z() + eps);
}

CGAL::Bbox_3 Messurer::TriangleSoup::bbox(const std::array<V_ind, 3> &v, double size) const {
    const Point& pp = points[v[0]];
    std::array<double, 3> p[2] = {{pp.x(), pp.y(), pp.z()}, {pp.x(), pp.y(), pp.z()}};
    for (int i = 1; i < 3; ++i) {
        const Point& q = points[v[i]];
        for (int j = 0; j < 3; ++j) {
            p[0][j] = (p[0][j] > q[j]) ? q[j] : p[0][j];
            p[1][j] = (p[1][j] < q[j]) ? q[j] : p[1][j];
        }
    }
    double eps = size / 2;
    for (int j = 0; j < 3; ++j){
        p[0][j] -= eps;
        p[1][j] += eps;
    }
    return CGAL::Bbox_3(p[0][0], p[0][1], p[0][2], p[1][0], p[1][1], p[1][2]);
}

CGAL::Bbox_3 Messurer::TriangleSoup::bbox(Messurer::F_ind f, double size) const {
    return bbox(tria.at(f), size);
}

World3d::Mesh Messurer::TriangleSoup::to_Mesh() const {
    std::vector<std::array<double, 3>> _points; _points.reserve(points.size());
    std::vector<std::vector<std::size_t>> _polygons; _polygons.reserve(tria.size());
    for (auto p: points) _points.push_back({p.x(), p.y(), p.z()});
    for (auto td: tria) { auto& t = td.second; _polygons.push_back({t[0], t[1], t[2]}); }
    return makeMesh(_points, _polygons);
}

World3d::Mesh Messurer::TriangleSoup::to_Mesh(function<bool(const TriangleSoup &, F_ind)> filter) {
    std::vector<std::array<double, 3>> _points; _points.reserve(points.size());
    std::vector<std::vector<std::size_t>> _polygons; _polygons.reserve(tria.size());
    for (auto p: points) _points.push_back({p.x(), p.y(), p.z()});
    for (auto td: tria) {
        if (!filter(*this, td.first)) continue;
        auto& t = td.second;
        _polygons.push_back({t[0], t[1], t[2]});
    }
    return makeMesh(_points, _polygons);
}

void Messurer::ObjectShape::Init(const Object3D* obj) {
    ObjectShape::obj = obj;
    mesh = TriangleSoup::makeTriangleSoup(obj);
}

using namespace std;
struct SewTwoFields{
    SewTwoFields(vector<Point_2>& p1, set<pair<int, int>>& link1,
                 vector<Point_2>& p2, set<pair<int, int>>& link2, map<int, int>& common_pnts_id):
            m_p1{p1}, m_p2{p2}, m_l1{link1}, m_l2{link2}, m_cn{common_pnts_id}
    {}
private:
    void set_data(gsl_vector *x);
    void set_new_points(gsl_vector* x);
    double compute(const gsl_vector *x);
    void compute_df(const gsl_vector *x, gsl_vector *df);
    double compute_fdf(const gsl_vector *x, gsl_vector *df);
    static void _fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df);
    static double _f4gsl (const gsl_vector *x, void *params);
    static void _df4gsl (const gsl_vector *x, void *params, gsl_vector *df);
public:
    int find_minimum_energy_df(int freq = 1, double step_sz = 1.e-4, double tol = 1.e-4, double epsabs = 1.e-2, int maxits = 50000, double time = 150);

private:
    vector<Point_2> &m_p1, &m_p2;
    set<pair<int, int>> &m_l1, &m_l2;
    map<int, int>& m_cn;
    const double m_wy = 10;
    map<int, int> m_map1, m_map2;
};

void SewTwoFields::set_data(gsl_vector *x){
    array<int, 4> off;
    off[0] = 0;
    off[1] = off[0] + m_p1.size() - m_cn.size();
    off[2] = off[1] + m_cn.size();
    off[3] = off[2] + m_p2.size() - m_cn.size();
    set<int> overlap1, overlap2;
    for (const auto& i: m_cn){
        overlap1.insert(i.first);
        overlap2.insert(i.second);
    }
    int id = 0;
    for (int i = 0; i < m_p1.size(); ++i){
        if (overlap1.count(i)) continue;
        m_map1.insert({i, id / 2});
        gsl_vector_set(x, id++, m_p1[i][0]);
        gsl_vector_set(x, id++, m_p1[i][1]);
    }
    for (auto i: overlap1){
        m_map1.insert({i, id / 2});
        m_map2.insert({m_cn[i], id / 2});
        gsl_vector_set(x, id++, m_p1[i][0]);
        gsl_vector_set(x, id++, m_p1[i][1]);
    }
    for (int i = 0; i < m_p2.size(); ++i){
        if (overlap2.count(i)) continue;
        m_map2.insert({i, id / 2});
        gsl_vector_set(x, id++, m_p2[i][0]);
        gsl_vector_set(x, id++, m_p2[i][1]);
    }
}

void SewTwoFields::set_new_points(gsl_vector* x){
    for(int i = 0, cnt = m_p1.size(); i < cnt; ++i) {
        int j = m_map1[i];
        m_p1[i] = Point_2(gsl_vector_get(x, 2 * j), gsl_vector_get(x, 2 * j + 1));
    }
    for(int i = 0, cnt = m_p2.size(); i < cnt; ++i) {
        int j = m_map2[i];
        m_p2[i] = Point_2(gsl_vector_get(x, 2 * j), gsl_vector_get(x, 2 * j + 1));
    }
}

double SewTwoFields::compute(const gsl_vector *x){
    double f = 0;
    auto adder = [&f = f, &m_wy = m_wy, &x = x](auto& container, auto& from, auto& _conv) {
        for (const auto &i: container) {
            double cur[4] = {from[i.first][0], from[i.first][1], from[i.second][0], from[i.second][1]};
            int idd[2] = {2*_conv[i.first],2*_conv[i.second]};
            double next[4] = {gsl_vector_get(x, idd[0]), gsl_vector_get(x, idd[0] + 1),
                              gsl_vector_get(x, idd[1]), gsl_vector_get(x, idd[1] + 1)};
            double weights[] = {1, m_wy};
            for (int k = 0; k < 2; ++k) {
                double loc = next[k] - next[k + 2] - (cur[k] - cur[k + 2]);
                f += loc * loc * weights[k];
                //cout << loc * loc * weights[k] << " ";
            }
        }
    };
    adder(m_l1, m_p1, m_map1);
    adder(m_l2, m_p2, m_map2);
    return f;
}

void SewTwoFields::compute_df(const gsl_vector *x, gsl_vector *df){
    auto adder = [&](auto& container, auto& from, auto& conv) {
        for (const auto &i: container) {
            int ids[2] = {i.first, i.second};
            int idd[2] = {2*conv[i.first], 2*conv[i.second]};
            double cur[4] = {from[i.first][0], from[i.first][1], from[i.second][0], from[i.second][1]};
            double next[4] = {gsl_vector_get(x, idd[0]), gsl_vector_get(x, idd[0] + 1),
                              gsl_vector_get(x, idd[1]), gsl_vector_get(x, idd[1] + 1)};
            double weights[] = {1, m_wy};
            for (int k = 0; k < 2; ++k) {
                double loc = next[k] - next[k + 2] - (cur[k] - cur[k + 2]);
                double ins = loc * weights[k];
                gsl_vector_set(df, idd[0] + k ,gsl_vector_get(df, idd[0] + k) + ins);
                gsl_vector_set(df, idd[1] + k,gsl_vector_get(df, idd[1] + k) - ins);
            }
        }
    };
    adder(m_l1, m_p1, m_map1);
    adder(m_l2, m_p2, m_map2);
}

double SewTwoFields::compute_fdf(const gsl_vector *x, gsl_vector *df){
    double f = 0;
    auto adder = [&](auto& container, auto& from, auto& conv) {
        for (const auto &i: container) {
            int ids[2] = {i.first, i.second};
            int idd[2] = {2*conv[i.first], 2*conv[i.second]};
            double cur[4] = {from[i.first][0], from[i.first][1], from[i.second][0], from[i.second][1]};
            double next[4] = {gsl_vector_get(x, idd[0]), gsl_vector_get(x, idd[0] + 1),
                              gsl_vector_get(x, idd[1]), gsl_vector_get(x, idd[1] + 1)};
            double weights[] = {1, m_wy};
            for (int k = 0; k < 2; ++k) {
                double loc = next[k] - next[k + 2] - (cur[k] - cur[k + 2]);
                double ins = loc * weights[k];
                f += ins * loc;
                gsl_vector_set(df, idd[0] + k ,gsl_vector_get(df, idd[0] + k) + ins);
                gsl_vector_set(df, idd[1] + k,gsl_vector_get(df, idd[1] + k) - ins);
            }
        }
    };
    adder(m_l1, m_p1, m_map1);
    adder(m_l2, m_p2, m_map2);
    return f;
}

void SewTwoFields::_fdf4gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    gsl_vector_set_all(df, 0);
    SewTwoFields* m = (SewTwoFields*)params;
    *f = m->compute_fdf(x, df);
}
double SewTwoFields::_f4gsl (const gsl_vector *x, void *params)
{
    SewTwoFields* m = (SewTwoFields*)params;

    return m->compute(x);
}
void SewTwoFields::_df4gsl (const gsl_vector *x, void *params, gsl_vector *df)
{
    gsl_vector_set_all(df, 0);
    SewTwoFields* m = (SewTwoFields*)params;
    m->compute_df(x, df);
}

int SewTwoFields::find_minimum_energy_df(int freq, double step_sz, double tol, double epsabs, int maxits, double time){
    class _Timer
    {
    private:
        using clock_t = std::chrono::high_resolution_clock;
        using second_t = std::chrono::duration<double, std::ratio<1> >;
        std::chrono::time_point<clock_t> m_beg;
    public:
        _Timer() : m_beg(clock_t::now()){}
        void reset() { m_beg = clock_t::now(); }

        double elapsed() const
        {
            return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
        }
    };

    size_t iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    void* par = this;

    int g_N = (m_p1.size() + m_p2.size() - m_cn.size()) * 2;
    gsl_vector *x = gsl_vector_alloc (g_N);
    set_data(x);

    gsl_multimin_function_fdf my_func;

    my_func.n = g_N;
    my_func.f = _f4gsl;
    my_func.df = _df4gsl;
    my_func.fdf = _fdf4gsl;
    my_func.params = par;

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, g_N);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_sz, tol);

    _Timer t;
    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, epsabs);

        if (status == GSL_SUCCESS && freq >= 0)
            printf ("Minimum found at: %5lu: f() = %10.5f\n", iter, s->f);

        if ((freq > 0) && !(iter % freq))
            printf ("%5lu: f() = %10.5f\n", iter, s->f);

    }
    while (status == GSL_CONTINUE && iter < maxits && t.elapsed() < time);
    if (status != GSL_SUCCESS && freq >= 0)
        printf ("Success didn't reached status = %d: %5lu: f() = %10.5f\n", status, iter, s->f);

    set_new_points(s->x);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);

    return status;
}

void Messurer::uniteMultiIntersectCuspHalves(int cusp_n, set<V_ind> &common, const std::array<std::set<V_ind>, 2>& coapt) {
//        std::vector<Point> points;
    //compute projections
    auto& shape = m_colission_shapes[cusp_n];
    std::array<std::map<V_ind, Point_2>, 2> projection;
    Point origin = shape.mesh.points[*common.begin()];
    for (int i = 0; i < 2; ++i)
        projection[i] = getHalfCoaptPointCloudProjection(cusp_n, (cusp_n + i + 1) % 3, coapt[i], origin);
    //compute existing links
    std::array<std::set<std::pair<V_ind, V_ind>>, 2> edges;
    for (int i = 0; i < 2; ++i) edges[i] = getHalfCoaptEdgesCloud(cusp_n, (cusp_n + 1 + i) % 3);

    //prepare data for using in energy minimizer
    std::array<std::vector<V_ind>, 2> remap;
    std::array<std::set<std::pair<int, int>>, 2> links;
    std::array<std::vector<Point_2>, 2> pnts;
    std::map<int, int> com;
    for (int i = 0; i < 2; ++i){
        std::map<V_ind, int> invremap;
        remap[i].reserve(projection[i].size());
        pnts[i].reserve(projection[i].size());
        int cnt = 0;
        for (auto d: projection[i]){
            remap[i].push_back(d.first);
            invremap[d.first] = remap[i].size() - 1;
            pnts[i].push_back(d.second);
        }
        for (auto e: edges[i])
            links[i].insert(std::pair<int, int>{invremap[e.first], invremap[e.second]});
        if (i == 0) for (auto c: common)
                com[invremap[c]] = c.idx();
        if (i == 1) for (auto& c: com)
                c.second = invremap[V_ind(c.second)];
    }

    SewTwoFields sw(pnts[0], links[0], pnts[1], links[1], com);
    sw.find_minimum_energy_df(-1, 1.e-4, 1.e-2, 5.e-1, 2000, 2);

    CoaptScan& scan = m_coaptScan[cusp_n];
    scan.remap.reserve(pnts[0].size() + pnts[1].size() - common.size());
    scan.pnts.reserve(pnts[0].size() + pnts[1].size() - common.size());
    std::map<V_ind, int> lremap;
    for (int i = 0; i < 2; ++i) {
        for (int p = 0; p < pnts[i].size(); ++p) {
            auto v = remap[i][p];
            if (i == 1 && common.count(v) > 0) continue;
            scan.pnts.push_back(pnts[i][p]);
            scan.remap.push_back(remap[i][p]);
            lremap[remap[i][p]] = scan.pnts.size() - 1;
        }
    }
    //reuse edges[0]
    for (auto q: edges[1])
        if (edges[0].count({q.second, q.first}) == 0)
            edges[0].insert(q);
    scan.edges.reserve(edges[0].size());
    for (auto e: edges[0])
        scan.edges.push_back({lremap[e.first], lremap[e.second]});
    //compute existing faces
    auto& lbl = m_colMap[cusp_n];
    auto& colFaces = scan.faces; colFaces.reserve(shape.mesh.tria.size());
    for (int j = 0; j < lbl.face_lbl.size(); ++j){
        int l = lbl.face_lbl[j];
        bool res = true;
        for (int k = 0; k < 3; ++k) {
            if (!((l & ((1 << (cusp_n + 2) % 3) << 3 * k)) || (l & ((1 << (cusp_n + 1) % 3) << 3 * k)))) res = false;
        }
        if (res) {
            auto& t = shape.mesh.tria[lbl.invmap[j]];
            colFaces.push_back({lremap[t[0]], lremap[t[1]], lremap[t[2]]});
        }
    }
}

double Messurer::computeHc() {
    std::array<std::set<V_ind>, 3> three_coapt;
    for (int i = 0; i < m_colission_shapes.size(); ++i){
        auto& shape = m_colission_shapes[i];
        auto& lbl = m_colMap[i];
        //find points coapting with both other cusps
        for (int j = 0; j < lbl.face_lbl.size(); ++j){
            int l = lbl.face_lbl[j];
            for (int k = 0; k < 3; ++k) {
                V_ind v = shape.mesh.tria[lbl.invmap[j]][k];
                if ((l & ((1 << (i + 1) % 3) << 3 * k)) && (l & ((1 << (i + 2) % 3) << 3 * k)))
                    three_coapt[i].insert(v);
            }
        }
    }
    bool empty = true;
    double bot = 0, up = 0;
    for (int i = 0; i < 3; ++i){
        if (!three_coapt[i].empty()){
            empty = false;
            auto& p = m_colission_shapes[i].mesh.points[*three_coapt[i].begin()];
            bot = up = (p - CGAL::ORIGIN) * m_direct;
            break;
        }
    }
    if (empty) return 0;
    for (int i = 0; i < 3; ++i){
        for (auto v: three_coapt[i]){
            auto& p = m_colission_shapes[i].mesh.points[v];
            double r = (p - CGAL::ORIGIN) * m_direct;
            if (r < bot) bot = r;
            if (r > up) up = r;
        }
    }
    return up - bot;
}

Messurer::CoaptationDistrib Messurer::CoaptScan::computeCoaptDistribution(int N, bool reject){
    CoaptationDistrib distr;
    CoaptScan& scan = *this;
    if (N <= 0) return distr;
    distr.down.resize(N+1, NAN);
    distr.up.resize(N+1, NAN);
    if (scan.edges.empty()) return distr;

    double l = scan.pnts[scan.edges.begin()->first].x(), r = scan.pnts[scan.edges.begin()->first].x();
    for (auto e: scan.edges){
        std::array<int, 2> v = {e.first, e.second};
        for (int k = 0; k < 2; ++k){
            double x = scan.pnts[v[k]].x();
            if (x < l) l = x;
            if (x > r) r = x;
        }
    }
    distr.width = r - l;
    r -= DBL_EPSILON*distr.width, l += DBL_EPSILON*distr.width;
    for (auto e: scan.edges){
        std::array<int, 2> v = {e.first, e.second};
        if (scan.pnts[v[0]].x() > scan.pnts[v[1]].x()) swap(v[0], v[1]);
        double l1 = scan.pnts[v[0]].x(), r1 = scan.pnts[v[1]].x();
        int kl = static_cast<int>(ceil((l1 - l) / (r - l) * N) + DBL_EPSILON),
                kr = static_cast<int>(floor((r1 - l) / (r - l) * N) + DBL_EPSILON);
        for (int k = kl; k <= kr; ++k){
            double w = (l + k * (r - l) / N - l1) / (r1 - l1);
            double y = (1 - w) * scan.pnts[v[0]].y() + w * scan.pnts[v[1]].y();
            if (isnan(distr.down[k]) || distr.down[k] > y) distr.down[k] = y;
            if (isnan(distr.up[k]) || distr.up[k] < y) distr.up[k] = y;
        }
    }
    auto remove_nan = [&](int k) {
        if (std::isnan(distr.down[k]) && std::isnan(distr.up[k]))distr.down[k] = distr.up[k] = 0;
        else if (std::isnan(distr.down[k])) distr.down[k] = distr.up[k];
        else if (std::isnan(distr.down[k])) distr.up[k] = distr.down[k];
    };
    remove_nan(0); remove_nan(N);

    if (reject)
        for (int i = 0; i < (N+1) / 2; ++i)
            swap(distr.down[i], distr.down[N - i]), swap(distr.up[i], distr.up[N - i]);

    return distr;
}

Mesh Messurer::CoaptScan::to_Mesh() const {
    std::vector<std::array<double, 3>> _points; _points.reserve(pnts.size());
    std::vector<std::vector<std::size_t>> _polygons; _polygons.reserve(faces.size());
    for (auto p: pnts) _points.push_back({p.x(), p.y(), 0});
    for (auto t: faces) {
        std::vector<std::size_t> pol(3);
        for (int i = 0; i < 3; ++i) pol[i] = t[i];
        _polygons.push_back(std::move(pol));
    }
    return makeMesh(_points, _polygons);
}

std::array<Messurer::CoaptationDistrib, 3> Messurer::computeCoaptDistribution(int N){
    std::array<Messurer::CoaptationDistrib, 3> res;
    for (int i = 0; i < 3; ++i){
        res[i] = m_coaptScan[i].computeCoaptDistribution(N);
    }
    return res;
}

std::array<Messurer::CoaptationDistrib, 3> Messurer::computeCoaptDistribution(std::array<int, 3> N){
    std::array<Messurer::CoaptationDistrib, 3> res;
    for (int i = 0; i < 3; ++i){
        res[i] = m_coaptScan[i].computeCoaptDistribution(N[i]);
    }
    return res;
}

Messurer::HalfCoaptDistrib Messurer::computeHalfCoaptDistrib(int Nparts){
    HalfCoaptDistrib res;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j){
        if (i == j) continue;
        res(i, j) = m_hcScan(i, j).computeCoaptDistribution(Nparts, j != (i+1)%3);
    }
    return res;
}

void Messurer::saveHalfCoaptDistribCSV(HalfCoaptDistrib& dist, std::string fname, bool by_row){
    std::ofstream csv(fname);
    if (!by_row) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                if (i == j) continue;
                for (int k = 0; k < 4; ++k)
                    csv << i << "->" << j << ", ";
            }
        csv << "\n";
        std::size_t max_r = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                if (i == j) continue;
                csv << "width, h, up, down,";
                max_r = std::max(max_r, dist(i, j).up.size());
            }
        csv << "\n";

        for (int l = 0; l < max_r; ++l) {
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    if (i == j) continue;
                    auto &up = dist(i, j).up;
                    auto &down = dist(i, j).down;
                    auto &width = dist(i, j).width;
                    int N = up.size();
                    if (l < N)
                        csv << l * width / (N - 1) << ", "
                            << (std::isnan(up[l] - down[l]) ? 0 : up[l] - down[l]) << ", "
                            << up[l] << ", "
                            << down[l] << ", ";
                    else
                        csv << ", , , , ";
                }
            csv << "\n";
        }
    } else {
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j){
            if (i == j) continue;
            auto& up = dist(i, j).up;
            auto& down = dist(i, j).down;
            auto& width = dist(i, j).width;
            int N = up.size();
            if (N == 0) {
                for (int k = 0; k < 4; ++k)
                csv << i << "->" << j << "\n";
            } else {
                csv << i << "->" << j << ", width";
                for (int i = 0; i < N; ++i)
                    csv << ", " << i * width / (N - 1);
                csv << "\n" << i << "->" << j <<  " ,h";
                for (int i = 0; i < N; ++i)
                    csv << ", " << (std::isnan(up[i] - down[i]) ? 0 : up[i] - down[i]);
                csv << "\n" << i << "->" << j << " ,up";
                for (int i = 0; i < N; ++i)
                    csv << ", " << up[i];
                csv << "\n" << i << "->" << j << " ,down";
                for (int i = 0; i < N; ++i)
                    csv << ", " << down[i];
                csv << "\n";
            }
        }
    }
    csv.close();
}

void Messurer::saveCoaptDistribCSV(array<CoaptationDistrib, 3> &dist, std::string fname, bool by_row) {
    std::ofstream csv(fname);
    if (by_row){
        for (int d = 0; d < 3; ++d) {
            auto &up = dist[d].up;
            auto &down = dist[d].down;
            auto &width = dist[d].width;
            int N = up.size();
            string dc = "l" + std::to_string(d);
            if (N == 0) {
                csv << dc << "\n\n\n";
            };
            csv << dc << ", width";
            for (int i = 0; i < N; ++i)
                csv << ", " << i * width / (N - 1);
            csv << "\n" << dc << " ,h";
            for (int i = 0; i < N; ++i)
                csv << ", " << up[i] - down[i];
            csv << "\n" << dc << " ,up";
            for (int i = 0; i < N; ++i)
                csv << ", " << up[i];
            csv << "\n" << dc << " ,down";
            for (int i = 0; i < N; ++i)
                csv << ", " << down[i];
            csv << "\n";
        }
    } else {
        for (int i = 0; i < 3; ++i)
                for (int k = 0; k < 4; ++k)
                    csv << "l" << i << ", ";
        csv << "\n";
        std::size_t max_r = 0;
        for (int i = 0; i < 3; ++i) {
            csv << "width, h, up, down,";
            max_r = std::max(max_r, dist[i].up.size());
        }
        csv << "\n";

        for (int l = 0; l < max_r; ++l) {
            for (int i = 0; i < 3; ++i){
                    auto &up = dist[i].up;
                    auto &down = dist[i].down;
                    auto &width = dist[i].width;
                    int N = up.size();
                    if (l < N)
                        csv << l * width / (N - 1) << ", "
                            << (std::isnan(up[l] - down[l]) ? 0 : up[l] - down[l]) << ", "
                            << up[l] << ", "
                            << down[l] << ", ";
                    else
                        csv << ", , , , ";
                }
            csv << "\n";
        }
    }
    csv.close();
}

void Messurer::computeMidPlanes() {
    using Triangle = CGAL::Simple_cartesian<double>::Triangle_3;
    std::vector<Triangle> triags;
    std::vector<Triangle> helper;
    std::array<Vector, 3> mid_normal;
    auto save_to_vector = [&m_colission_shapes = m_colission_shapes, &m_colMap = m_colMap, &mid_normal](std::vector<Triangle>& triags, int i, int col_id) {
        short mask = ((1 | 8 | 64) << col_id);//(1 << (col_id + 3 * 0)) | (1 << (col_id + 3 * 1)) | (1 << (col_id + 3 * 2));
        auto& P = m_colission_shapes[i].mesh.points;
        for (auto &fi: m_colission_shapes[i].mesh.tria) {
            auto t = fi.second;
            if (m_colMap[i].face_lbl[m_colMap[i].remap[fi.first]] & mask) {
                triags.emplace_back(P[t[0]], P[t[1]], P[t[2]]);
                mid_normal[i] += CGAL::cross_product(P[t[1]] - P[t[0]], P[t[2]] - P[t[0]]);
            }
        }
    };
    Point C[2]; int sz[2] = {-1, -1};
    for (int i = 0; i < 3; ++i)
        for (int j = i+1; j < 3; ++j){
            triags.clear(); helper.clear(); std::fill(mid_normal.begin(), mid_normal.end(), CGAL::NULL_VECTOR);
            save_to_vector(triags, i, j);
            if ((sz[0] = triags.size()) > 1)
                C[0] = CGAL::centroid(triags.begin(), triags.end(), CGAL::Dimension_tag<2>());
            save_to_vector(helper, j, i);
            //debug saving
//            write_STL(triags, "../result/MeassureDebug/plane_field_"+std::to_string(i)+"_"+std::to_string(j)+".stl");
//            write_STL(helper, "../result/MeassureDebug/plane_field_"+std::to_string(j)+"_"+std::to_string(i)+".stl");
            if ((sz[1] = helper.size()) > 1)
                C[1] = CGAL::centroid(helper.begin(), helper.end(), CGAL::Dimension_tag<2>());
            for (auto& t: helper) triags.push_back(std::move(t));
            if (triags.size() > 2) {
                Plane_3 plane;
                CGAL::linear_least_squares_fitting_3(triags.begin(), triags.end(), plane,
                                                     CGAL::Dimension_tag<2>());
                if (sz[0] > 1) {
                    m_planes(i, j).second = true;
                    Vector normal = plane.orthogonal_vector();
                    if (normal * mid_normal[i] < 0) normal *= -1;
                    m_planes(i, j).first = Plane_3(C[0], normal);
                } else m_planes(i, j).second = false;
                if (sz[1] > 1) {
                    m_planes(j, i).second = true;
                    Vector normal = plane.orthogonal_vector();
                    if (normal * mid_normal[j] < 0) normal *= -1;
                    m_planes(j, i).first = Plane_3(C[1], normal);
                } else m_planes(j, i).second = false;

            } else m_planes(j, i).second = m_planes(i, j).second = false;
        }
}

bool Messurer::isCoaptWith(int on_id, int collide_id){
    if (on_id == collide_id) return false;
    if (on_id < 0 || collide_id < 0 || on_id > 2 || collide_id > 2) return false;
    auto& shape = m_colission_shapes[on_id];
    auto& lbl = m_colMap[on_id];
    for (int r = 0; r < lbl.face_lbl.size(); ++r){
        int l = lbl.face_lbl[r];
        for (int k = 0; k < 3; ++k) {
            V_ind v = shape.mesh.tria[lbl.invmap[r]][k];
            if (l & ((1 << collide_id) << 3 * k)) return true;
        }
    }
    return false;
}

bool Messurer::isTripleCollide(int on_id){
    if (on_id < 0 || on_id > 2) return false;
    auto& shape = m_colission_shapes[on_id];
    auto& lbl = m_colMap[on_id];
    //find points coapting with both other cusps
    for (int j = 0; j < lbl.face_lbl.size(); ++j){
        int l = lbl.face_lbl[j];
        for (int k = 0; k < 3; ++k) {
            V_ind v = shape.mesh.tria[lbl.invmap[j]][k];
            if ((l & ((1 << (on_id + 1) % 3) << 3 * k)) && (l & ((1 << (on_id + 2) % 3) << 3 * k)))
                return true;
        }
    }
    return false;
}

std::set<V_ind> Messurer::getHalfCoaptShapePointCloud(int on_id, int collide_id){
    std::set<V_ind> coapt;
    if (on_id == collide_id) return coapt;
    if (on_id < 0 || collide_id < 0 || on_id > 2 || collide_id > 2) return coapt;

    auto& shape = m_colission_shapes[on_id];
    auto& lbl = m_colMap[on_id];
    for (int r = 0; r < lbl.face_lbl.size(); ++r){
        int l = lbl.face_lbl[r];
        for (int k = 0; k < 3; ++k) {
            V_ind v = shape.mesh.tria[lbl.invmap[r]][k];
            if (l & ((1 << collide_id) << 3 * k)) coapt.insert(v);
        }
    }
    return coapt;
}

std::map<V_ind, Point_2> Messurer::getHalfCoaptPointCloudProjection(int on_id, int collide_id, const std::set<V_ind>& cloud){
    //compute projections
    const auto& coapt = cloud;
    std::map<V_ind, Point_2> projection;
    if (coapt.empty()) return projection;

    return getHalfCoaptPointCloudProjection(on_id, collide_id, cloud, m_colission_shapes[on_id].mesh.points[*coapt.begin()]);
}

std::map<V_ind, Point_2> Messurer::getHalfCoaptPointCloudProjection(int on_id, int collide_id, const std::set<V_ind>& cloud, Point origin){
    //compute projections
    const auto& coapt = cloud;
    std::map<V_ind, Point_2> projection;
    if (coapt.empty()) return projection;

    auto& shape = m_colission_shapes[on_id];
    auto& lbl = m_colMap[on_id];
    auto plane = m_planes(on_id, collide_id).first;
    Vector e2 = plane.projection(CGAL::ORIGIN + m_direct) - plane.projection(CGAL::ORIGIN);
    e2 /= sqrt(e2.squared_length());
    Vector e1 = CGAL::cross_product(plane.orthogonal_vector(), e2);
    //e1 directed from center of valve to border if j == (i + 1)%3 else inverse direction
    e1 /= sqrt(e1.squared_length());
    for (auto v: coapt)
        projection[v] = Point_2(e1 * (shape.mesh.points[v] - origin), e2 * (shape.mesh.points[v] - origin));
    return projection;
}

std::set<std::pair<V_ind, V_ind>> Messurer::getHalfCoaptEdgesCloud(int on_id, int collide_id){
    std::set<std::pair<V_ind, V_ind>> edges;
    auto& shape = m_colission_shapes[on_id];
    auto& lbl = m_colMap[on_id];
    for (int r = 0; r < lbl.face_lbl.size(); ++r){
        int l = lbl.face_lbl[r];
        for (int k1 = 0; k1 < 3; ++k1)
            for (int k2 = k1 + 1; k2 < 3; ++k2) {
                V_ind v1 = shape.mesh.tria[lbl.invmap[r]][k1], v2 = shape.mesh.tria[lbl.invmap[r]][k2];
                if ((l & ((1 << collide_id) << 3 * k1))  && (l & ((1 << collide_id) << 3 * k2)))
                    if (edges.find({v2, v1}) == edges.end())
                        edges.insert({v1, v2});
            }
    }
    return edges;
}

void Messurer::alignOyPointCloud(std::map<V_ind, Point_2>& cloud, bool right){
    Vector_2 dp;
    if (right) {
        auto p = std::min_element(cloud.begin(), cloud.end(),
                                  [](const auto &a, const auto &b) { return a.second.x() < b.second.x(); });
        dp = Vector_2(p->second.x(), p->second.y());
    } else {
        auto p = std::max_element(cloud.begin(), cloud.end(),
                                  [](const auto &a, const auto &b) { return a.second.x() < b.second.x(); });
        dp = Vector_2(p->second.x(), p->second.y());
    }
    for(auto& pr: cloud) pr.second -= dp;
}

Messurer::CoaptScan Messurer::getHalfCoaptScan(int i, int j, const std::set<V_ind>& coapt){
    CoaptScan res;
    if (i == j) return res;
    if (coapt.empty()) return res;
    auto& shape = m_colission_shapes[i];
    auto& lbl = m_colMap[i];
    std::map<V_ind, Point_2> projection = getHalfCoaptPointCloudProjection(i, j, coapt);
    alignOyPointCloud(projection, j == (i+1)%3);
    std::set<std::pair<V_ind, V_ind>> edges = getHalfCoaptEdgesCloud(i, j);
    std::vector<V_ind> remap;
    std::vector<std::pair<int, int>> links;
    std::vector<std::array<int, 3>> faces;
    std::vector<Point_2> pnts;
    std::map<V_ind, int> invremap;
    remap.reserve(projection.size());
    pnts.reserve(projection.size());
    links.reserve(edges.size());
    faces.reserve(2 * edges.size() / 3);
    for (auto d: projection){
        remap.push_back(d.first);
        invremap[d.first] = remap.size() - 1;
        pnts.push_back(d.second);
    }
    for (auto e: edges)
        links.emplace_back(invremap[e.first], invremap[e.second]);
    for (int r = 0; r < lbl.face_lbl.size(); ++r){
        int l = lbl.face_lbl[r];
        bool resf = true;
        for (int k = 0; k < 3; ++k) { resf = resf && (l & ((1 << j) << 3 * k)); }
        if (resf) {
            auto& t = shape.mesh.tria[lbl.invmap[r]];
            faces.push_back({invremap[t[0]], invremap[t[1]], invremap[t[2]]});
        }
    }
    res = CoaptScan{std::move(pnts), std::move(links), std::move(remap), std::move(faces)};
    return res;
}

Messurer::CoaptScan Messurer::getHalfCoaptScan(int i, int j){
    if (i == j) return CoaptScan();
    std::set<V_ind> coapt = getHalfCoaptShapePointCloud(i, j);
    return getHalfCoaptScan(i, j, coapt);
}

void Messurer::computeHalfCoaptScans() {
    for (int i = 0; i < m_colission_shapes.size(); ++i){
        auto& shape = m_colission_shapes[i];
        auto& lbl = m_colMap[i];
        for (int j = 0; j < 3; ++j){
            //compute coapting points
            if (j == i) continue;
            std::set<V_ind> coapt = getHalfCoaptShapePointCloud(i, j);
            if (coapt.empty()) continue;
            //compute projections
            std::map<V_ind, Point_2> projection = getHalfCoaptPointCloudProjection(i, j, coapt);
            alignOyPointCloud(projection, j == (i+1)%3);
            std::set<std::pair<V_ind, V_ind>> edges = getHalfCoaptEdgesCloud(i, j);
            //save data to CoaptScan
            std::vector<V_ind> remap;
            std::vector<std::pair<int, int>> links;
            std::vector<std::array<int, 3>> faces;
            std::vector<Point_2> pnts;
            std::map<V_ind, int> invremap;
            remap.reserve(projection.size());
            pnts.reserve(projection.size());
            links.reserve(edges.size());
            faces.reserve(2 * edges.size() / 3);
            for (auto d: projection){
                remap.push_back(d.first);
                invremap[d.first] = remap.size() - 1;
                pnts.push_back(d.second);
            }
            for (auto e: edges)
                links.emplace_back(invremap[e.first], invremap[e.second]);
            for (int r = 0; r < lbl.face_lbl.size(); ++r){
                int l = lbl.face_lbl[r];
                bool res = true;
                for (int k = 0; k < 3; ++k) { res = res && (l & ((1 << j) << 3 * k)); }
                if (res) {
                    auto& t = shape.mesh.tria[lbl.invmap[r]];
                    faces.push_back({invremap[t[0]], invremap[t[1]], invremap[t[2]]});
                }
            }
            m_hcScan(i, j) = CoaptScan{std::move(pnts), std::move(links), std::move(remap), std::move(faces)};
            //debug saving
//            Object3D(m_hcScan(i, j).to_Mesh()).save("../result/MeassureDebug/coapt_proj_"+std::to_string(i)+"_"+std::to_string(j)+".stl");
        }
    }
}

Messurer::Billowing Messurer::computeBillowing() {
    V_ind av_min[3] = {V_ind(-1), V_ind(-1), V_ind(-1)};
    for (int i = 0; i < 3; ++i){
        auto& m = *m_cusps[i];
        V_ind& v_min = av_min[i];
        for (auto v: m.m_mesh.vertices()){
            if (!m.is_movable(v)) {
                v_min = v;
                break;
            }
        }
        if (v_min == V_ind(-1)) continue;
        double _min = (m.m_x[v_min] - CGAL::ORIGIN) * m_direct;
        for (auto v: m.m_mesh.vertices()){
            if (!m.is_movable(v)) {
                double _lmin = (m.m_x[v] - CGAL::ORIGIN) * m_direct;
                if (_lmin < _min){
                    _min = _lmin;
                    v_min = v;
                }
            }
        }
        m_bill.plane[i] = m.m_x[v_min];
    }
    if (av_min[0] == V_ind(-1) || av_min[1] == V_ind(-1) || av_min[2] == V_ind(-1)) {
        std::fill(m_bill.plane.begin(), m_bill.plane.end(), CGAL::ORIGIN);
        return m_bill;
    }
    Vector n = CGAL::cross_product(m_bill.plane[1] - m_bill.plane[0], m_bill.plane[2] - m_bill.plane[0]);
    n /= sqrt(n.squared_length());
    V_ind bv_min[3] = {av_min[0], av_min[1], av_min[2]};
    for (int i = 0; i < 3; ++i){
        double min = 0;
        auto& m = *m_cusps[i];
        for (auto v: m.m_mesh.vertices()){
            double cur = (m.m_x[v] - CGAL::ORIGIN) * n;
            if (cur < min) {
                min = cur;
                bv_min[i] = v;
            }
        }
        m_bill.bilPnt[i] = m.m_x[bv_min[i]];
    }
    return m_bill;
}

bool Messurer::isValveClosed() {
    for (int i = 0; i < m_colission_shapes.size(); ++i){
        auto& shape = m_colission_shapes[i];
        auto& lbl = m_colMap[i];
        //find points coapting with both other cusps
        for (int j = 0; j < lbl.face_lbl.size(); ++j){
            int l = lbl.face_lbl[j];
            for (int k = 0; k < 3; ++k) {
                V_ind v = shape.mesh.tria[lbl.invmap[j]][k];
                if ((l & ((1 << (i + 1) % 3) << 3 * k)) && (l & ((1 << (i + 2) % 3) << 3 * k)))
                    return true;
            }
        }
    }
    return false;
}

void Messurer::uniteSingleIntersectCuspHalves(int cusp_n, V_ind common, const array<std::set<V_ind>, 2> &coapt) {
    //compute projections
    auto& shape = m_colission_shapes[cusp_n];
    std::array<std::map<V_ind, Point_2>, 2> projection;
    Point origin = shape.mesh.points[common];
    for (int i = 0; i < 2; ++i)
        projection[i] = getHalfCoaptPointCloudProjection(cusp_n, (cusp_n + i + 1) % 3, coapt[i], origin);
    //compute existing links
    std::array<std::set<std::pair<V_ind, V_ind>>, 2> edges;
    for (int i = 0; i < 2; ++i) edges[i] = getHalfCoaptEdgesCloud(cusp_n, (cusp_n + 1 + i) % 3);

    CoaptScan& scan = m_coaptScan[cusp_n];
    scan.remap.reserve(coapt[0].size() + coapt[1].size() - 1);
    scan.pnts.reserve(coapt[0].size() + coapt[1].size() - 1);
    std::map<V_ind, int> lremap;
    for (int i = 0; i < 2; ++i) {
        for (auto v: coapt[i]){
            if (i == 1 && common == v) continue;
            scan.pnts.push_back(projection[i][v]);
            scan.remap.push_back(v);
            lremap[v] = scan.pnts.size() - 1;
        }
    }
    //reuse edges[0]
    for (auto q: edges[1])
        if (edges[0].count({q.second, q.first}) == 0)
            edges[0].insert(q);
    scan.edges.reserve(edges[0].size());
    for (auto e: edges[0])
        scan.edges.emplace_back(lremap[e.first], lremap[e.second]);
    //compute existing faces
    auto& lbl = m_colMap[cusp_n];
    auto& colFaces = scan.faces; colFaces.reserve(shape.mesh.tria.size());
    for (int j = 0; j < lbl.face_lbl.size(); ++j){
        int l = lbl.face_lbl[j];
        bool res = true;
        for (int k = 0; k < 3; ++k) {
            if (!((l & ((1 << (cusp_n + 2) % 3) << 3 * k)) || (l & ((1 << (cusp_n + 1) % 3) << 3 * k)))) res = false;
        }
        if (res) {
            auto& t = shape.mesh.tria[lbl.invmap[j]];
            colFaces.push_back({lremap[t[0]], lremap[t[1]], lremap[t[2]]});
        }
    }
}

int Messurer::computeCoaptStatus() {
    int res =   (isCoaptWith(0, 1) << 0) |
                (isCoaptWith(0, 2) << 1) |
                (isCoaptWith(1, 0) << 2) |
                (isCoaptWith(1, 2) << 3) |
                (isCoaptWith(2, 0) << 4) |
                (isCoaptWith(2, 1) << 5) |
                (isValveClosed() << 6);
    return res;
}

void Messurer::uniteNonIntersectCuspHalves(int cusp_n, array<std::set<V_ind>, 2> &coapts, bool with_gap) {
    if (coapts[0].empty() && coapts[1].empty()) return;
    auto& shape = m_colission_shapes[cusp_n];
    auto& lbl = m_colMap[cusp_n];
    for (int k = 0; k < 2; ++k)
        if (coapts[k].empty()){
            int j = (cusp_n + 1 + (k + 1)%2)%3;
            auto& coapt = coapts[(k+1)%2];
            m_coaptScan[cusp_n] = getHalfCoaptScan(cusp_n, j, coapt);
            return;
        }
    std::array<CoaptScan, 2> scans = {getHalfCoaptScan(cusp_n, (cusp_n + 1)%3, coapts[0]),
                                      getHalfCoaptScan(cusp_n, (cusp_n + 2)%3, coapts[1])};
    if (with_gap) {
        const double eps = 0.0025;
        auto compInterval = [](const CoaptScan &c, int axis = 0) -> std::pair<double, double> {
            auto it = std::minmax_element(c.pnts.begin(), c.pnts.end(),
                                          [axis](const Point_2 &p1, const Point_2 &p2) {
                                              return p1[axis] < p2[axis];
                                          });
            return {it.first->operator[](axis), it.second->operator[](axis)};
        };
        std::array<std::pair<double, double>, 2> X = {compInterval(scans[0]), compInterval(scans[1])},
                Y = {compInterval(scans[0], 1), compInterval(scans[1], 1)};
        double minXshift =
                (X[1].second - X[0].first) + eps * ((X[0].second - X[0].first) + (X[1].second - X[1].first));

        double Xshift = minXshift;
        double Yshift = 0;
        auto preq = m_colission_shapes[cusp_n].makePointRequier("v:point");
        {
            auto it1 = std::minmax_element(scans[0].remap.begin(), scans[0].remap.end(),
                                           [&preq, axis = 0](const V_ind &v1, const V_ind &v2) -> bool {
                                               auto p1 = preq(v1), p2 = preq(v2);
                                               return p1[axis] < p2[axis];
                                           });
            auto it2 = std::minmax_element(scans[1].remap.begin(), scans[1].remap.end(),
                                           [&preq, axis = 0](const V_ind &v1, const V_ind &v2) -> bool {
                                               auto p1 = preq(v1), p2 = preq(v2);
                                               return p1[axis] < p2[axis];
                                           });
            double ix[4] = {preq(*it1.first)[0], preq(*it1.second)[0], preq(*it2.first)[0], preq(*it2.second)[0]};
            if (ix[0] > ix[3]) {
                long vid1 = it1.first.base() - scans[0].remap.data(), vid2 =
                        it2.second.base() - scans[1].remap.data();
                if (ix[0] - ix[3] > minXshift) Xshift = ix[0] - ix[3];
                Yshift = (preq(*it1.first)[1] - preq(*it2.second)[1]) -
                         (scans[0].pnts[vid1][1] - scans[1].pnts[vid2][1]);
            } else if (ix[2] > ix[1]) {
                long vid1 = it1.second.base() - scans[0].remap.data(), vid2 =
                        it2.first.base() - scans[1].remap.data();
                if (ix[2] > ix[1] > minXshift) Xshift = ix[2] > ix[1];
                Yshift = (preq(*it1.second)[1] - preq(*it2.first)[1]) -
                         (scans[0].pnts[vid1][1] - scans[1].pnts[vid2][1]);
            }
        }
        Vector_2 dp{Xshift / 2, Yshift / 2};
        for (auto &p: scans[0].pnts) p += dp;
        for (auto &p: scans[1].pnts) p -= dp;
    }
    auto& scan = m_coaptScan[cusp_n];
    scan = std::move(scans[0]);
    int poff = scan.pnts.size();
    scan.pnts.resize(scan.pnts.size() + scans[1].pnts.size());
    std::copy(scans[1].pnts.begin(), scans[1].pnts.end(), scan.pnts.begin() + poff);
    scan.remap.resize(scan.remap.size() + scans[1].remap.size());
    std::copy(scans[1].remap.begin(), scans[1].remap.end(), scan.remap.begin() + poff);
    scan.edges.reserve(scan.edges.size() + scans[1].edges.size());
    for (auto e: scans[1].edges) scan.edges.emplace_back(e.first + poff, e.second + poff);
    scan.faces.reserve(scan.faces.size() + scans[1].faces.size());
    for (auto& t: scans[1].faces) scan.faces.push_back({t[0] + poff, t[1] + poff, t[2] + poff});
    return;
}

void Messurer::uniteCollidingCuspHalves() {
    for (int i = 0; i < m_colission_shapes.size(); ++i){
        auto& shape = m_colission_shapes[i];
        auto& lbl = m_colMap[i];
        //find common points
        std::set<V_ind> common;
        std::array<std::set<V_ind>, 2> coapt;
        for (int j = 0; j < lbl.face_lbl.size(); ++j){
            int l = lbl.face_lbl[j];
            for (int k = 0; k < 3; ++k) {
                V_ind v = shape.mesh.tria[lbl.invmap[j]][k];
                if (l & ((1 << (i + 1) % 3) << 3 * k)) coapt[0].insert(v);
                if (l & ((1 << (i + 2) % 3) << 3 * k)) coapt[1].insert(v);
                if ((l & ((1 << (i + 1) % 3) << 3 * k)) && (l & ((1 << (i + 2) % 3) << 3 * k)))
                    common.insert(v);
            }
        }
        if (common.size() == 0)
            uniteNonIntersectCuspHalves(i, coapt);
        else if (common.size() == 1)
            uniteSingleIntersectCuspHalves(i, *common.begin(), coapt);
        else if (common.size() > 1)
            uniteMultiIntersectCuspHalves(i, common, coapt);
        //debug save
//            Object3D(m_coaptScan[i].to_Mesh()).save("../result/MeassureDebug/coapt_unite_distr_"+std::to_string(i)+".stl");
    }
}

array<Messurer::ColArea, 3> Messurer::computeColArea() {
    std::array<ColArea, 3> res;
    for (int i = 0; i < m_colission_shapes.size(); ++i) {
        auto &shape = m_colission_shapes[i];
        auto &lbl = m_colMap[i];

        for (int r = 0; r < lbl.face_lbl.size(); ++r){
            int l = lbl.face_lbl[r];
            bool iscoapt = true;
            bool isbicol = true;
            bool iscoapt1 = true, iscoapt2 = true;
            for (int k = 0; k < 3 && iscoapt; ++k) {
                bool coapt1 = l & (1 << ((i+1)%3 + 3 * k)),
                        coapt2 = l & (1 << ((i+2)%3 + 3 * k) );
                if ( !coapt1 && !coapt2 ) { iscoapt = false;  isbicol = false;}
                else if ( !coapt1 || !coapt2 ) isbicol = false;
                if (!coapt1) iscoapt1 = false;
                if (!coapt2) iscoapt2 = false;
            }
            auto v = shape.mesh.tria[lbl.invmap[r]];
            auto& pp = shape.mesh.points;
            if (iscoapt) {
                double area = sqrt(CGAL::cross_product(pp[v[0]] - pp[v[1]], pp[v[0]] - pp[v[2]]).squared_length()) / 2;
                res[i].colArea += area;
                if (isbicol) res[i].biColArea += area;
                if (iscoapt1) res[i].partColArea[0] += area;
                if (iscoapt2) res[i].partColArea[1] += area;
            }
        }
    }

    return res;
}

std::ostream &operator<<(ostream &out, const Messurer::ColArea &a) {
    double bicol = a.biColArea;
    return out << "Full collide area = " << a.colArea << " (" << a.partColArea[0] - bicol << " + " << a.partColArea[1] - bicol << " + " << bicol << " + " << a.halfColArea() <<")";
}
