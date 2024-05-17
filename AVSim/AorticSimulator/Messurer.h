//
// Created by alex on 29.11.2020.
//

#ifndef AORTIC_VALVE_MESSURER_H
#define AORTIC_VALVE_MESSURER_H

#include "AVSim/Core/Object3D.h"
#include "AVSim/Core/TriangularMeshHelpers.h"

#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <fstream>

class Messurer {
public:
    using Object3D = World3d::Object3D;
    using Point = World3d::Point;
    using Point_2 = World3d::Point_2;
    using Vector = World3d::Vector;
    using Vector_2 = World3d::Vector_2;
    using V_ind = World3d::V_ind;
    using E_ind = World3d::E_ind;
    using F_ind = World3d::F_ind;

    struct TriangleSoup{
        typedef std::array<V_ind, 3> Triangle;
        std::vector<Point> points;
        std::map<F_ind, Triangle> tria;

        static TriangleSoup makeTriangleSoup(const Object3D* obj);
        static CGAL::Bbox_3 bbox(const Point &p1, const Point &p2, const Point &p3, double size);
        static CGAL::Bbox_3 bbox(const Point &p, double size);
        CGAL::Bbox_3 bbox(V_ind v, double size) const;
        CGAL::Bbox_3 bbox(const std::array<V_ind, 3>& v, double size) const;
        CGAL::Bbox_3 bbox(F_ind f, double size) const;
        World3d::Mesh to_Mesh() const;
        World3d::Mesh to_Mesh(std::function<bool(const TriangleSoup& ts, F_ind f)> filter);
    };

    struct ObjectShape{
        struct ExtraFace{
            F_ind f;    //face index of main face
            std::array<Vector, 3> w;//baricentric coords of owing points
        };
        const Object3D*  obj = nullptr;  //here initial mesh
        TriangleSoup mesh;  //here refined mesh
        std::map<F_ind, ExtraFace> extra_Finfo; //saved info to place point on initial mesh

        ObjectShape() {};
        ObjectShape(const Object3D* obj) { Init(obj); };
        void Init(const Object3D* obj);

        struct PointRequier{
            const ObjectShape& shape;
            std::map<V_ind, std::pair<F_ind, int>> shape_links;
            Mesh::Property_map<V_ind, Point> m_x;
            PointRequier(ObjectShape& _shape, std::string tag):
                shape{_shape}, m_x{_shape.obj->m_mesh.property_map<V_ind, Point>(tag).first}{
                for (auto t: shape.mesh.tria){
                    auto f = t.first;
                    for (int k = 0; k < 3; ++k)
                        if (shape_links.find(t.second[k]) == shape_links.end()) shape_links[t.second[k]] = {f, k};
                }
            }
            Point operator()(V_ind v){
                auto fi = shape_links[v];
                auto f = fi.first;
                int k = fi.second;
                std::array<Vector, 3> w = {Vector{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                auto fit = shape.extra_Finfo.find(f);
                if (fit != shape.extra_Finfo.end()) {
                    w = fit->second.w;
                    f = fit->second.f;
                }
                auto vv = World3d::vert_around(shape.obj->m_mesh, f);
                Point PPT[3] = {m_x[vv[0]], m_x[vv[1]], m_x[vv[2]]};
                auto evalP = [](const Vector& w, const Point& p1, const Point& p2, const Point& p3) -> Point{
                    return Point(w[0]*p1.x()+w[1]*p2.x() + w[2]*p3.x(), w[0]*p1.y()+w[1]*p2.y() + w[2]*p3.y(), w[0]*p1.z()+w[1]*p2.z() + w[2]*p3.z());
                };
                Point PP = evalP(w[k], PPT[0], PPT[1], PPT[2]);
                return PP;
            }
        };
        PointRequier makePointRequier(std::string tag){ return PointRequier(*this, tag); }
    };

    struct CollisionPairData{
        V_ind v; //colliding point
        F_ind f; //collided face
        Point fp; //closest point on f
        double dist2; //signed distance
    };

    struct CollisionInfo{
        typedef CollisionPairData CPD;
        std::array<std::array<std::vector<CollisionPairData>, 2>, 3> data;

        std::vector<CollisionPairData>& getCollisionInfo(int from, int to);
        std::vector<CollisionPairData>& operator()(int from, int to){ return getCollisionInfo(from, to); }
    };
    template <typename T>
    struct PairSharedData{
        using Handle = T;
        std::array<std::array<Handle, 2>, 3> data;
        Handle& operator()(int from, int to){
            assert(from != to);
            return data[from][(to - from + 3) % 3 - 1];
        }
    };

    using Plane_3 = CGAL::Simple_cartesian<double>::Plane_3;
    using PlaneData = PairSharedData<std::pair<Plane_3, bool>>;

    struct CollisionMap{
        std::map<F_ind, std::size_t> remap;
        std::vector<F_ind> invmap;
        std::vector<unsigned short> face_lbl; //F_ind -> lbl
        void clear() {
            remap.clear();
            invmap.clear();
            face_lbl.clear();
        }
    };
    struct CoaptationDistrib{
        std::vector<double> up;
        std::vector<double> down;
        double width = 0;

        void saveCSV(std::string fname){
            std::ofstream csv(fname);
            int N = up.size();
            if (N == 0) return;
            csv << "0";
            for (int i = 1; i < N; ++i)
                csv << ", " << i * width / (N - 1);
            csv << "\n" << up[0];
            for (int i = 1; i < N; ++i)
                csv << ", " << up[i];
            csv << "\n" << down[0];
            for (int i = 1; i < N; ++i)
                csv << ", " << down[i];
            csv << "\n";
            csv.close();
        }
    };
    struct CoaptScan{
        std::vector<Point_2> pnts;
        std::vector<std::pair<int, int>> edges;
        std::vector<V_ind> remap; //just for debug, it isn't used to compute distribution
        std::vector<std::array<int, 3>> faces;  //just for debug, it isn't used to compute distribution

        CoaptationDistrib computeCoaptDistribution(int N, bool reject=false); //reject will reject result relative Oy axis
        //debug function to save current Scan
        [[nodiscard]] Mesh to_Mesh() const;
    };
    using HalfCoaptScan = PairSharedData<CoaptScan>;
    struct HalfCoaptDistrib: public PairSharedData<CoaptationDistrib>{
        double getCoaptH(int on_id, int to_id){
            if (on_id == to_id || on_id < 0 || to_id < 0 || on_id > 2 || to_id > 2) return -1;
            double h = 0;
            auto& distr = operator()(on_id, to_id);
            for (int n = 0; n < distr.up.size(); ++n){
                double dh = distr.up[n] - distr.down[n];
                if (!std::isnan(dh) && dh > h) h = dh;
            }
            return h;
        }
        double getCoaptH(){
            double h = 0;
            for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j){
                if (i == j) continue;
                auto& distr = operator()(i, j);
                for (int n = 0; n < distr.up.size(); ++n){
                    double dh = distr.up[n] - distr.down[n];
                    if (!std::isnan(dh) && dh > h) h = dh;
                }
            }
            return h;
        }
    };
    struct Billowing{
        std::array<Point, 3> plane = {CGAL::ORIGIN, CGAL::ORIGIN, CGAL::ORIGIN};
        std::array<Point, 3> bilPnt;
        double getBillowing(int cusp_id){
            if (plane[0] == CGAL::ORIGIN && plane[1] == CGAL::ORIGIN && plane[2] == CGAL::ORIGIN) return NAN;
            Plane_3 p(plane[0], plane[1], plane[2]);
            double res = sqrt((bilPnt[cusp_id] - p.projection(bilPnt[cusp_id])).squared_length());
            return res;
        }
        std::array<double, 3> getBillowing(){
            return std::array<double, 3>{getBillowing(0), getBillowing(1), getBillowing(2)};
        }
        void print(std::ostream& out){
            out << "billow plate: {" << plane[0] << ", " << plane[1] << ", " << plane[2] << "}\n";
            out << "billow pnts : {" << bilPnt[0] << ", " << bilPnt[1] << ", " << bilPnt[2] << "}\n";
            auto bil = getBillowing();
            out << "billowing   : {" << bil[0] << ", " << bil[1] << ", " << bil[2] << "}\n";
        }
        void print() { print(std::cout); }
    };
    struct ColArea{
        double colArea = 0;
        double partColArea[2] = {0, 0};
        double biColArea = 0;
        double halfColArea() const { return colArea - (partColArea[0] + partColArea[1] - biColArea); }
    };
private:
    Object3D& m_aorta;
    std::array<Object3D*, 3> m_cusps;
    Vector m_direct;
    double m_margin = 0.0;

public:
    std::array<ObjectShape, 3> m_colission_shapes;
    CollisionInfo m_colission_info;
    std::array<CollisionMap, 3> m_colMap;
    PlaneData m_planes;
    std::array<CoaptScan, 3> m_coaptScan;
    HalfCoaptScan m_hcScan;
    Billowing m_bill;

    Messurer(Object3D& aorta, Object3D& NCC, Object3D& RCC, Object3D& LCC, Vector orientation = {0, 0, 1}):
            m_aorta{aorta}, m_cusps{&NCC, &RCC, &LCC}, m_direct{orientation} {
        for (int i = 0; i < 3; ++i)
            m_colission_shapes[i] = ObjectShape(m_cusps[i]);
    }
    Messurer& setMargin(double margin) { m_margin = margin; return *this; }
    void computeCollidingBnd(int level, std::function<bool(const ObjectShape& from, const ObjectShape& to, std::vector<CollisionInfo::CPD>& info)> spec = nullptr);
    [[maybe_unused]] static bool saveObjectShape(const ObjectShape& obj, std::string fname, std::string tag = "v:x",
                                                 std::function<bool (F_ind)> filter = [](F_ind) { return true; });

    void computeMidPlanes();
    void computeHalfCoaptScans();
    Billowing computeBillowing();
    std::array<ColArea, 3> computeColArea();
    bool isValveClosed();
    std::array<Object3D*, 3> getCusps() { return m_cusps; }

    void uniteCollidingCuspHalves();
    //coapt[i] - set of vertices colliding with (cusp_n + i + 1)%3 cusp
    void uniteNonIntersectCuspHalves(int cusp_n, std::array<std::set<V_ind>, 2>& coapts, bool with_gap = true);
    //common - common vertex of the cusp with id = cusp_n
    //coapt[i] - set of vertices colliding with (cusp_n + i + 1)%3 cusp including common vertex
    void uniteSingleIntersectCuspHalves(int cusp_n, V_ind common, const std::array<std::set<V_ind>, 2>& coapt);
    //common - set of common vertices of the cusp with id = cusp_n
    //coapt[i] - set of vertices colliding with (cusp_n + i + 1)%3 cusp including common vertices
    void uniteMultiIntersectCuspHalves(int cusp_n, std::set<V_ind>& common, const std::array<std::set<V_ind>, 2>& coapt);

    double computeHc();
    bool isCoaptWith(int on_id, int collide_id);
    bool isTripleCollide(int on_id);
    int computeCoaptStatus();
    HalfCoaptDistrib computeHalfCoaptDistrib(int Nparts);
    std::array<CoaptationDistrib, 3> computeCoaptDistribution(int Nparts);
    std::array<CoaptationDistrib, 3> computeCoaptDistribution(std::array<int, 3> Nparts);
    static void saveCoaptDistribCSV(std::array<CoaptationDistrib, 3>& dist, std::string fname, bool by_row = false);
    static void saveHalfCoaptDistribCSV(HalfCoaptDistrib& dist, std::string fname, bool by_row = false);

private:
    bool checkCollisionShapes();
    std::set<V_ind> getHalfCoaptShapePointCloud(int on_id, int collide_id);
    std::map<V_ind, Point_2> getHalfCoaptPointCloudProjection(int on_id, int collide_id, const std::set<V_ind>& cloud);
    std::map<V_ind, Point_2> getHalfCoaptPointCloudProjection(int on_id, int collide_id, const std::set<V_ind>& cloud, Point origin);
    std::set<std::pair<V_ind, V_ind>> getHalfCoaptEdgesCloud(int on_id, int collide_id);
    static void alignOyPointCloud(std::map<V_ind, Point_2>& cloud, bool right = true);

    CoaptScan getHalfCoaptScan(int on_id, int collide_id);
    CoaptScan getHalfCoaptScan(int on_id, int collide_id, const std::set<V_ind>& cloud);
};

std::ostream& operator<<(std::ostream& out, const Messurer::ColArea& a);



#endif //AORTIC_VALVE_MESSURER_H
