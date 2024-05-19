//
// Created by alex on 17.11.2020.
//

#include "AVSimulator.h"
#include "AVSim/LeafTemplates/TemplatesCollection.h"
#include "AVSim/Core/Collision/BulletCollisionManager.h"
#include "AVSim/Core/ForceAppliers/ForceAppliers.h"
#include <fstream>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Circle_3.h>
#include <CGAL/Plane_3.h>
#include "Sewer/MinEnergyDeformator.h"
#include "AVSim/Core/Renderers/Renderer.h"

#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/AABB_segment_primitive.h>
#include "AVSim/Core/NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverKinsol.h"
#include "AVSim/Core/NonLinearSolverCustom.h"

using namespace World3d;

AorticSurface fromMesh(const World3d::Mesh& m){
    AorticSurface as;
    as.nT = m.num_faces();
    as.tri = new uint32_t[as.nT*3];
    as.nP = m.num_vertices();
    as.pnt = new float[as.nP*3];
    for (auto f: m.faces()) {
        auto v = World3d::vert_around(m, f);
        for (int n = 0; n < 3; ++n)
            as.tri[3*f.idx()+n] = v[n].idx();
    }
    for (auto v: m.vertices())
        for(int x = 0; x < 3; ++x)
            as.pnt[3*v.idx() + x] = m.point(v)[x];
    return as;
}

void freeAorticSurface(AorticSurface* as){
    delete[] as->tri;
    delete[] as->pnt;
}

void InterConnector::update(IntermediateResults& ir){
    auto it = _interData.begin();
    //пропускаем первый элемент, потому что в нём лежит расчётная сетка поверхности аорты
    if (it != _interData.end())
        ++it;
    int cnt = _interData.size() - 1;
    if (cnt > ir.nL){
        if (ir.leaflets) delete[] ir.leaflets;
        ir.leaflets = new IntermediateResults::Leaflet[cnt];
    }
    ir.nL = cnt;
    int i = 0;
    while (it != _interData.end()){
        ir.leaflets[i].state = static_cast<uint8_t>(it->second.second);
        ir.leaflets[i].nm_sz = it->second.first.name.size();
        ir.leaflets[i].name = const_cast<char*>(it->second.first.name.data());
        if (it->second.first.value.second) {
            ir.leaflets[i].val_sz = it->second.first.value.first.size();
            ir.leaflets[i].value = it->second.first.value.first.data();
        }
        if (it->second.first.vertex.second) {
            ir.leaflets[i].v_sz = it->second.first.vertex.first.size();
            ir.leaflets[i].vertex = reinterpret_cast<float *>(it->second.first.vertex.first.data());
        }
        if (it->second.first.face.second){
            ir.leaflets[i].f_sz = it->second.first.face.first.size();
            ir.leaflets[i].faces = it->second.first.face.first.data();
        }
        if (it->second.first.vertex_lbl.second){
            int sz = ir.leaflets[i].v_sz;
            if (ir.leaflets[i].vlb_sz < sz){
                ir.leaflets[i].vlb_sz = sz;
                if (ir.leaflets[i].v_label) delete[] ir.leaflets[i].v_label;
                ir.leaflets[i].v_label = new int32_t[sz];
            }
            std::fill(ir.leaflets[i].v_label, ir.leaflets[i].v_label + sz, 0);
            for (auto& j: it->second.first.vertex_lbl.first)
                ir.leaflets[i].v_label[j.first] = j.second;
        }
        //в будущем я добавлю поле availableValues в GuiBackInterconnection
        ir.leaflets[i].av_sz = 0;
        ir.leaflets[i].availableValues = nullptr;
        ++it; ++i;
    }
}

void InterConnector::destructIntermediateResults(IntermediateResults& ir) {
    if (ir.leaflets) {
        for (int i = 0; i < ir.nL; ++i)
            if (ir.leaflets[i].v_label)
                delete[] ir.leaflets[i].v_label;
        delete[] ir.leaflets;
    }
}

SewLines readSewLines(std::string file){
    ifstream f(file);
    if (! f.good()) {
        throw std::runtime_error("Unsuccessful attempt to open file \"" + file + "\"");
    }
    SewLines sl;
    int cnt = 0;
    f >> cnt;
    sl.nL = cnt;
    sl.lines = new SewLines::Line[sl.nL];
    for (int i = 0; i < sl.nL; ++i){
        f >> sl.lines[i].nP;
        sl.lines[i].pnt = new float[3*sl.lines[i].nP];
        for (int n = 0; n < sl.lines[i].nP; ++n)
            f >> sl.lines[i].pnt[3*n + 0] >> sl.lines[i].pnt[3*n + 1] >> sl.lines[i].pnt[3*n + 2];
    }
    return sl;
}

void freeSewLines(SewLines* sl){
    for (int i = 0; i < sl->nL; ++i){
        delete[] sl->lines[i].pnt;
    }
    delete[] sl->lines;
}

void AVSimulator::setAorticSurface(const AorticSurface &as) {
    auto& points = asd.points;
    auto& polygons = asd.polygons;
    points.resize(as.nP);
    polygons.resize(as.nT);
    for (int i = 0; i < as.nP; ++i) points[i] = {as.pnt[3*i + 0], as.pnt[3*i + 1], as.pnt[3*i + 2]};
    for (int i = 0; i < as.nT; ++i) polygons[i] = {as.tri[3*i+0], as.tri[3*i+1], as.tri[3*i+2]};
}

void AVSimulator::setSewLines(const SewLines &sw) {
    sld.lines.resize(sw.nL);
    for (int i = 0; i < sw.nL; ++i){
        sld.lines[i].resize(sw.lines[i].nP);
        auto& line = sld.lines[i];
        for (int j = 0; j < sw.lines[i].nP; ++j) {
            line[j] = SewLineData::Node{sw.lines[i].pnt[3 * j], sw.lines[i].pnt[3 * j + 1], sw.lines[i].pnt[3 * j + 2]};
        }
    }
    auto& com_p = sld.commissure_points;
    for (int i = 0; i < 3; ++i)
        com_p[i] = CGAL::ORIGIN + (sld.lines[i][0] + sld.lines[(i + 2) % 3][sld.lines[(i + 2) % 3].size() - 1]) / 2;
    sld.orientation = CGAL::cross_product(com_p[1] - com_p[0], com_p[2] - com_p[0]);
    sld.orientation /= sqrt(sld.orientation.squared_length());
}

//эта функция вырезает часть аорты, которая попала в цилиндр, который описывает линии пришивания
//входные параметры масштабируют размеры цилиндра (lecu и lecd отвечают за изменение высоты), а rec - радиуса
Object3D AVSimulator::make_aorta_object(double lecu, double lecd, double rec) {
    assert(("SewLine is not setted or has less then 3 lines!",
            sld.lines.size() != 0
            && sld.lines[0].size() != 0 && sld.lines[1].size() != 0 && sld.lines[2].size() != 0));
    typedef CGAL::Simple_cartesian<DReal>::Circle_3 Circle_3;
    Circle_3 c(sld.commissure_points[0], sld.commissure_points[1], sld.commissure_points[2]);
    Point p0 = c.center();
    double min = DBL_MAX, max = DBL_MIN, maxR2 = c.squared_radius();
    for (auto& line: sld.lines)
        for (auto& pnt: line){
            double proj = ((CGAL::ORIGIN + pnt) - p0) * sld.orientation;
            min  = (min > proj) ? proj : min;
            max  = (max < proj) ? proj : max;
            double R2 = (((CGAL::ORIGIN + pnt) - p0) - proj * sld.orientation).squared_length();
            maxR2 = (maxR2 < R2) ? R2 : maxR2;
        }
    double mid = (max + min) / 2, half = (max - min) / 2;

    double actR2 = maxR2 * rec * rec, actmin = mid  - half * lecd, actmax = mid + half * lecu;
    double actmid = (actmin + actmax) / 2, acthalf = (actmax - actmin) / 2;

    set<int> inserted_p, inserted_f;
    for (int ip = 0; ip < asd.polygons.size(); ++ip){
        auto& polygon = asd.polygons[ip];
        std::array<Point, 3> tri;
        for (int i = 0; i < 3; ++i) {
            int id = polygon[i];
            auto& p = asd.points[id];
            tri[i] = Point(p[0], p[1], p[2]);
        }
        bool insert = false;
        for (int i = 0; i < 3 && !insert; ++i){
            double proj = (tri[i] - p0) * sld.orientation;
            double R2 = ((tri[i] - p0) - proj * sld.orientation).squared_length();
            if (abs(proj - actmid) < acthalf && R2 < actR2) insert = true;
        }
        if (insert){
            inserted_f.insert(ip);
            for (int i = 0; i < 3; ++i)
                inserted_p.insert(polygon[i]);
        }
    }
    std::vector<array<float, 3>> points(inserted_p.size());
    map<uint32_t, uint32_t> map_p;
    int new_id = 0;
    for (auto it = inserted_p.begin(); it != inserted_p.end(); ++it){
        points[new_id] = asd.points[*it];
        map_p[*it] = new_id++;
    }
    std::vector<AorticSurfaceData::CGAL_Polygon> polygons;
    polygons.reserve(inserted_f.size());
    for (auto ip: inserted_f){
        auto& pp = asd.polygons[ip];
        polygons.push_back({map_p[pp[0]], map_p[pp[1]], map_p[pp[2]]});
    }

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
    PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(Array_traits()));
    PMP::orient_polygon_soup(points, polygons);

    Mesh m;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, m);
    Object3D res(m);
    res.name = "aorta";

    return res;
}

void AVSimulator::transformLeafsIntoAorta(const vector<Object3D*>& leafs, std::array<Vector, 3> orients){
    assert(leafs.size() == 3 && sld.lines.size() == 3 && "This function valid only for 3-flap valve");
    for (auto& o: orients) o /= sqrt(o.squared_length());
    vector<double> heights = {0, 0, 0};
    for (int l = 0; l < sld.lines.size(); ++l){
        auto& line = sld.lines[l];
        for (auto& node: line){
            double h = (line[0] - node) * orients[l];
            if (heights[l] < h) heights[l] = h;
        }
    }

    for (int l = 0; l < 3; ++l){
        auto& leaf = *(leafs[l]);
        Point p[3];
        for (auto v: leaf.m_mesh.vertices()){
            if (leaf.m_boundary[v] == (2 | 1) )
                p[0] = leaf.m_x[v];
            else if (leaf.m_boundary[v] == (2 | 4) )
                p[1] = leaf.m_x[v];
        }
        Vector e0 = p[1] - p[0];
        e0 /= sqrt(e0.squared_length());
        //here we invert mesh orientation
        Vector e1(e0[1], -e0[0], 0);
        for (auto v: leaf.m_mesh.vertices()){
            leaf.m_x[v] = Point((leaf.m_x[v]-CGAL::ORIGIN)*e0, (leaf.m_x[v]-CGAL::ORIGIN)*e1, 0);
        }

        p[2] = *leaf.m_x.begin();
        for (auto v: leaf.m_mesh.vertices()){
            if (leaf.m_boundary[v] == (2 | 1) )
                p[0] = leaf.m_x[v];
            else if (leaf.m_boundary[v] == (2 | 4) )
                p[1] = leaf.m_x[v];
            if (p[2].y() > leaf.m_x[v].y())
                p[2] = leaf.m_x[v];
        }

        double h = (p[0].y() + p[1].y()) / 2 - p[2].y();
        double wy = heights[l] / h;
        double wx = sqrt((sld.commissure_points[(l+1)%3] - sld.commissure_points[l]).squared_length() / (p[1] - p[0]).squared_length());
        e0 = sld.commissure_points[(l+1)%3] - sld.commissure_points[l];
        e0 /= sqrt(e0.squared_length());
        e1 = orients[l];
        Point origin = sld.commissure_points[l];
        for (auto v: leaf.m_mesh.vertices()){
            leaf.m_x[v] = origin + (leaf.m_x[v][0] - p[0][0])* wx * e0 + (leaf.m_x[v][1] - p[0][1])* wy * e1;
        }
    }
}

//scale and rotate flat templates to insert their into aorta
//it is guaranteed that the inserted flaps do not intersect with each other
void AVSimulator::transformLeafsIntoAorta(const Object3D& aorta, const vector<Object3D*>& leafs){
    transformLeafsIntoAorta(leafs, {sld.orientation, sld.orientation, sld.orientation});
}

void AVSimulator::sewLeafBndToSewBnd(const vector<Object3D*>& leafs){
    for (int l = 0; l < 3 && l < leafs.size(); ++l){
        if (!leafs[l]) continue;
        auto& leaf = *(leafs[l]);
        auto& line = sld.lines[l];
        double l_s = 0, l_l = 0;
        std::vector<V_ind> bnd;
        //compute sequential vertices of boundary on template
        {
            std::map<V_ind, std::vector<E_ind>> bnd_map;
            E_ind bnd_start;
            for (auto e: leaf.m_mesh.edges()) {
                auto v = vert_around(leaf.m_mesh, e);
                if ((leaf.m_boundary[v[0]] & 2) && (leaf.m_boundary[v[1]] & 2)) {
                    l_l += sqrt((leaf.m_x0[v[1]] - leaf.m_x0[v[0]]).squared_length());
                    if (leaf.m_boundary[v[0]] == (2 | 1) || leaf.m_boundary[v[1]] == (2 | 1))
                        bnd_start = e;
                    bnd_map[v[0]].push_back(e);
                    bnd_map[v[1]].push_back(e);
                }
            }
            bnd.reserve(bnd_map.size());
            auto bnd_v = vert_around(leaf.m_mesh, bnd_start);
            V_ind v = bnd_v[0];
            bnd.push_back(bnd_v[1]);
            if (bnd_map[bnd_v[0]].size() < 2) {
                v = bnd_v[1];
                bnd[0] = bnd_v[0];
            }
            bnd.push_back(v);
            while (bnd.size() < bnd_map.size()){
                auto& bb = bnd_map[v];
                if (bb.size() != 2)
                    throw std::runtime_error("Template " + to_string(l) +" has non connected boundary or non standart boundary label: "
                                             + std::to_string(leaf.m_boundary[v]) + " v = "  + std::to_string(v));
                if (bb[0] != bnd_start)
                    bnd_start = bb[0];
                else
                    bnd_start = bb[1];
                bnd_v = vert_around(leaf.m_mesh, bnd_start);
                if (bnd_v[0] != v)
                    v = bnd_v[0];
                else
                    v = bnd_v[1];
                bnd.push_back(v);
            }
        }
        for (int i = 0; i < line.size() - 1; ++i)
            l_s += sqrt((line[i + 1] - line[i]).squared_length());
        double J0 = l_s / l_l;
        //linear mapping template boundary -> sew lines
        if (sld.cdfs.size() < l){
            double s1 = 0, s0 = 0;
            double d1 = 0, d0 = 0;
            for (int i = 1, j = 0; i < bnd.size(); ++i) {
                s0 = s1;
                s1 += J0 * sqrt((leaf.m_x0[bnd[i]] - leaf.m_x0[bnd[i - 1]]).squared_length());
                while (d1 <= s0 && j < line.size() - 1){
                    d0 = d1;
                    d1 += sqrt((line[j+1] - line[j]).squared_length());
                    j++;
                }

                double ww = (s0 - d0) / (d1 - d0);
                leaf.m_x[bnd[i-1]] = CGAL::ORIGIN + line[j - 1] * (1 - ww) + line[j] * ww;
            }
            leaf.m_x[bnd.back()] = CGAL::ORIGIN + line.back();
        } else {
            std::vector<double> s; s.reserve(bnd.size());
            s.push_back(0.0);
            for (int i = 1; i < bnd.size()-1; ++i)
                s.push_back(s.back() + (sqrt((leaf.m_x0[bnd[i]] - leaf.m_x0[bnd[i - 1]]).squared_length()))/l_l);
            s.push_back(1.0);
            if (!sld.cdfs[l](s, J0, l_l)) throw std::runtime_error("Error in CDF function [" + to_string(l) + "] for " + leafs[l]->name);

            leaf.m_x[bnd.front()] = CGAL::ORIGIN + line.front();
            double d1 = 0, d0 = 0;
            for (int i = 1, j = 0; i < s.size()-1; ++i){
                while (d1 <= s[i] && j < line.size() - 1){
                    d0 = d1;
                    j++;
                    d1 += sqrt((line[j] - line[j-1]).squared_length()) / l_s;
                }
                double ww = (s[i] - d0) / (d1 - d0);
                leaf.m_x[bnd[i]] = CGAL::ORIGIN + line[j - 1] * (1 - ww) + line[j] * ww;
            }
            leaf.m_x[bnd.back()] = CGAL::ORIGIN + line.back();
        }
    }

    for (int l = 0; l < 3 && l < leafs.size(); ++l)
        for (auto v: leafs[l]->m_mesh.vertices())
            leafs[l]->m_next_x[v] = leafs[l]->m_x[v];
}

//linear mapping template boundary -> sew lines
void AVSimulator::sewLeafletsCoarse(const Object3D& aorta, const vector<Object3D*>& leafs){
    transformLeafsIntoAorta(aorta, leafs);
    sewLeafBndToSewBnd(leafs);
}

//надо бы создать реализацию которая будет отыскивать разбиение
void AVSimulator::sewLeafletsMid(const Object3D& aorta, const vector<Object3D*>& leafs){
    assert(leafs.size() == 3 && sld.lines.size() == 3 && "This function valid only for 3-flap valve");
    typedef CGAL::Simple_cartesian<DReal>::Circle_3 Circle_3;
    typedef CGAL::Simple_cartesian<DReal>::Plane_3 Plane_3;
    Circle_3 c(sld.commissure_points[0], sld.commissure_points[1], sld.commissure_points[2]);
    Point p0 = c.center();
    std::vector<Plane_3> planes(3);
    for (int i = 0; i < planes.size(); ++i){
        planes[i] = Plane_3(p0 + sld.orientation, p0, sld.commissure_points[i]);
    }
    bool separate = true;
    for (int l = 0; l < 3; ++l){
        auto h0 = planes[l].orthogonal_vector(), h1 = -planes[(l+1)%3].orthogonal_vector();
        h0 /= sqrt(h0.squared_length()), h1 /= sqrt(h1.squared_length());
        for (int xi = 0; xi < sld.lines[l].size(); ++xi){
            auto& x = sld.lines[l][xi];
            auto p = x - (p0 - CGAL::ORIGIN);
            double q0 = p * h0, q1 = p * h1;
            if (q0 > 0 || q1 > 0)
                separate = false;
        }
    }
    if (!separate)
        std::cout << "Warning: This boundaries is not separable by this ad-hoc algorithm" << std::endl;

    for (int l = 0; l < 3; ++l){
        EnergyDeformatorOb3D ed_obj(*(leafs[l]));
//        ed_obj.setRenderer(std::make_unique<World3d::DefaultRenderer>());
        MinEnergyDeformator m(ed_obj);
        Configuration::SewEnergyParams& sep = conf.sep;
        auto h0 = planes[l].orthogonal_vector(), h1 = -planes[(l+1)%3].orthogonal_vector();
        h0 /= sqrt(h0.squared_length()), h1 /= sqrt(h1.squared_length());
        set_plane_constr(m, sep.plane_w, {h0[0], h0[1], h0[2]}, (p0 - CGAL::ORIGIN)*h0);
//        std::cout << "  plane: n = " << h0 << ", x0 = " << p0 << "\n";
        set_plane_constr(m, sep.plane_w, {h1[0], h1[1], h1[2]}, (p0 - CGAL::ORIGIN)*h1);
//        std::cout << "  plane: n = " << h1 << ", x0 = " << p0 << "\n";
        set_default_length_constr(m, sep.sqr_length_w);
        set_default_digedral_angle_constr(m, sep.digedral_angle_w, sep.convexity_w);
        auto cgal_wind = (-0.3 * sld.orientation + 1*(0.5 * h0 + 0.5 * h1));
        Eigen::Vector3d wind {cgal_wind[0], cgal_wind[1], cgal_wind[2]};
        wind.normalize();
        set_isotrop_force(m, conf.sep.force_w, wind);
//        std::cout << "  plane: wind = " << wind.transpose() << std::endl;
        auto& emp = sep.sp;
        m.find_minimum_energy_df(emp.freq, emp.step_sz, emp.tol, emp.epsabs, emp.maxits, emp.time);
    }

    for (int l = 0; l < 3; ++l)
        for (auto v: leafs[l]->m_mesh.vertices())
            leafs[l]->m_next_x[v] = leafs[l]->m_x[v];
}

void AVSimulator::sewLeaflets(const Object3D &aorta, vector<Object3D*> &leafs) {
    sewLeafletsCoarse(aorta, leafs);
    sewLeafletsMid(aorta, leafs);
}

void AVSimulator::read_leafs() {
    vector<Object3D> leaf(sld.lines.size());
    for (int i = 0; i < leaf.size(); ++i){
        if (conf.l_ins[i].generative) {
            AniMesh am;
            auto& prm = conf.l_ins[i].model_params;
            if (conf.l_ins[i].model == "old_ozaki")
                am = generate_old_ozaki(prm.value("template", 17), prm.value("mesh_h", 1.2));
            else if (conf.l_ins[i].model == "ozaki") {
                OzakiTempl tpl(prm.value("template", 21.0));
                auto sw_it = prm.find("suture_width"); if (sw_it != prm.end()) tpl.s = sw_it->get<double>();
                auto mw_it = prm.find("margins_width"); if (mw_it != prm.end()) tpl.m = mw_it->get<double>();
                auto h_it = prm.find("h"); if (h_it != prm.end()) tpl.h = h_it->get<double>();
                auto beta_it = prm.find("beta"); if (beta_it != prm.end()) tpl.beta = beta_it->get<double>();

                am = generate_ozaki(tpl, prm.value("mesh_h", 1.2));
            } else if (conf.l_ins[i].model == "custom_ozaki"){
                OzakiTempl tpl;
                tpl.D = prm.value("D", 25.0);
                tpl.h = prm.value("h", 11.0 + erf(2*(tpl.D-16.0)));
                tpl.beta = prm.value("beta", 2.5 + 0.5*erf(2*(tpl.D-24.0)));
                tpl.alpha = prm.value("alpha", 1.0);
                tpl.s = prm.value("s", 0.0);
                tpl.m = prm.value("m", 0.0);
                tpl.w = prm.value("w", 0.0);

                am = generate_ozaki(tpl, prm.value("mesh_h", 1.2));
            } else
                throw std::runtime_error("Other cases of templates are not implemented");
            Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
            leaf[i] = Object3D{m};
            leaf[i].name = "leaf_" + to_string(i);
        } else {
            leaf[i].read(conf.bs_dir + conf.l_ins[i].leaflet_file);
        }
    }
    static const int FIXED = 2, FREE = 1|4|8;
    DirichletCond dc = [](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & FIXED) return 0;
        return 7;
    };
    for (int i = 0; i < 3; ++i) {
        m_lid[i] = m_w.addObject3D(move(leaf[i]), 1);
        m_w.obj(m_lid[i]).setDirichletCond(dc);
    }
}

void AVSimulator::_run_computations(GuiBackInterconnection *ic) {
    World &w = m_w;
    auto& aorta = m_w.obj(m_aid);

    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BulletCollisionManager>();
    w.setCollider(move(colMan));
    auto renderer = std::make_unique<World3d::DefaultRenderer>();
    renderer->set_interconnector(ic);
    m_w.setRenderer(std::move(renderer));
    std::array<World3d::ObjectID, 3>& l_id = m_lid;
    std::array<Force, 3>  bending;
    for (int i = 0; i < 3; ++i) {
        reinterpret_cast<BulletCollisionManager*>(w.getCollider())->set_margin(l_id[i],conf.ots[i].collision_margin);
        Force Pr = SimplePressureLoad(conf.ots[i].pressure);
        w.addForce(Pr, l_id[i]);

        string dir = "../generated";
        Force elastic;
        auto& prm = conf.ots[i].elastic_model_params;
        if (conf.ots[i].elastic_model_type == "NeoGookModel"){
            double E = prm.value("E", 1.0e6);
            double mu = E/3, H = 0.5;
            elastic = NeoGookModel(prm.value("mu", mu), prm.value("H", H), dir);
        } else if (conf.ots[i].elastic_model_type == "GentModel"){
            double E = prm.value("E", 1.0e6);
            double mu = E/3, Jm = 2.3, H = 0.5;
            elastic = GentModel(prm.value("mu", mu), prm.value("H", H), prm.value("Jm", Jm), dir);
        } else if (conf.ots[i].elastic_model_type == "SVKirchhoffModel"){
            double E = prm.value("E", 1.0e6), nu = prm.value("nu", 0.5);
            double lambda = E * nu / (1 - nu * nu), mu = E / (2 * (1 + nu)), H = 0.5;
            elastic = SVKirchhoffModel(prm.value("H", H), prm.value("lambda", lambda), prm.value("mu", mu), dir);
        } else
            throw std::runtime_error("Other cases of forces are not implemented");
        if (conf.ots[i].use_bending)
            bending[i] = BendingForce(elastic.target<HyperElasticForceBase>()->f, true);

        w.addForce(elastic, l_id[i]);
    }
    //add bending force
    for (int i = 0; i < 3; ++i) {
        if (conf.ots[i].use_bending) w.addForce(bending[i], l_id[i]);
    }
    if (conf.sp.method == "simple_relaxation") {
        auto& prm = conf.sp.method_vals;
        double delta = prm.value("delta", 5.5e-7)/10;
        w.setForceApplier(StaticForceApplierP(delta));
        int maxits = prm.value("maxits", 72000);
        double err = prm.value("eps", 1e-3);
        int iter1 = 399; //first iter1 iteration for relaxation of initial guess
        int freq = 200;
        auto relax_time = World3d::Timer();
        w.Simulation([&err, iter1, &freq, &maxits, &out = std::cout, &tm = relax_time, &delta, &prm](StepSimInfo& info, World* w)->bool{
            static double resid_init = 1;
            int it = info.it;
            if (it == iter1) {
                delta = prm.value("delta", 5.5e-7);
            }
//                else if (it == 2000) {
//                    delta = 1e-4;
//                }
            if (it == iter1 + 1) resid_init = w->getWorldShift();
            if ((it > iter1 && it % freq == 0) || (it >= maxits)){
                double resid = w->getWorldShift();
                double eps = resid / resid_init;
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << tm.elapsed() << "\n";

                if (eps < err){
                    out << "Algorithm is converged: \n";
                    return true;
                }
                else if (it >= maxits){
                    out << "Algorithm have reached the maximum of iteration: " << maxits << "\n";
                    return true;
                }
                else if (eps > 300 || std::isnan(eps)){
                    out << "Algorithm is diverged: \n";
                    exit(-1);
                }
            }

            return false;
        });
    } else
        throw std::runtime_error("Other cases of solver method are not implemented");
}

void AVSimulator::_make_postsavings() {
    auto saveTag = [](const Object3D& obj, string filename, string tag){
        ofstream ob(filename, ios::binary | ios::trunc);
        if (!ob) return false;
        auto m_x = obj.m_mesh.property_map<V_ind, Point>(tag).first;
        for (auto v: obj.m_mesh.vertices()){
            std::array<double, 3> P = {m_x[v][0], m_x[v][1], m_x[v][2]};
            ob.write(reinterpret_cast<const char *> ((P.data())), 3 * sizeof(double));
        }
        ob.close();
        return true;
    };

    string& result_dir = conf.svr.result_dir;
    for (int i = 0; i < 3; ++i) {
        m_w.obj(m_lid[i]).save(conf.svr.result_dir + "leaf_res_" + to_string(i) + conf.svr.postfix + ".stl");
        m_w.obj(m_lid[i]).save(conf.svr.result_dir + "leaf_res_" + to_string(i) + conf.svr.postfix + ".txt");
        saveTag(m_w.obj(m_lid[i]), conf.svr.result_dir + "leaf_res" + to_string(i) + "_m_x0" + conf.svr.postfix + ".tag", "v:point");
    }
    std::cout << "Main direction: " << sld.orientation << std::endl;
}

void AVSimulator::OccludeAorticValve(int argc, const char **argv, const AorticSurface &as, const SewLines &sw,
                                     GuiBackInterconnection *ic) {
    conf.ConfigurationInit(argc, argv);
    setSewLines(sw);       //read sewed lines from input
    setAorticSurface(as);  //read aorta from input
    set_aorta_object();    //cut off excess parts of the aorta
    read_leafs();          //read/generate leaflets
    auto& aorta = getAorta();
    aorta.save(conf.svr.result_dir + "aorta_cut" + conf.svr.postfix + ".stl");
    auto leaf = getLeafs();
    for (int i = 0; i < 3; ++i) leaf[i]->save(conf.svr.result_dir + "leaf_init_"+to_string(i)+conf.svr.postfix + ".stl");
    sewLeafletsCoarse(aorta, leaf); //sew leafs to aorta
    for (int i = 0; i < 3; ++i) leaf[i]->save(conf.svr.result_dir + "leaf_sew_"+to_string(i)+conf.svr.postfix + ".stl");

    _run_computations(ic);

    _make_postsavings();
}

void AVSimulator::fillCoaptationMessures(CoaptationResult &cr) {
    Messurer mes(m_w.obj(m_aid), m_w.obj(m_lid[0]), m_w.obj(m_lid[1]), m_w.obj(m_lid[2]), sld.orientation);
    mes.setMargin(1.8*2*conf.ots[0].collision_margin);
    mes.computeCollidingBnd(4);
    mes.computeMidPlanes();
    string to_save = "../result/";
    for (int i = 0; i < 3; ++i) {
        auto colfilter = [&mes, &col = mes.m_colMap[i]](F_ind f) {
            auto lbl = col.face_lbl[col.remap[f]];
            bool res = true;
            for (int k = 0; k < 3; ++k)
                res &= ((lbl & (7 << 3 * k)) > 0);
            return res;
        };
        auto nocolfilter = [&colfilter](F_ind f) { return !colfilter(f); };
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + m_w.obj(m_lid[i]).name + "_mes_init.stl", "v:point");
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + m_w.obj(m_lid[i]).name + "_mes_init_col.stl", "v:point",
                            colfilter);
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + m_w.obj(m_lid[i]).name + "_mes_init_nocol.stl", "v:point",
                            nocolfilter);
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + m_w.obj(m_lid[i]).name + "_mes_real.stl", "v:x");
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + m_w.obj(m_lid[i]).name + "_mes_real_col.stl", "v:x", colfilter);
    }
    mes.computeHalfCoaptScans();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            Object3D(mes.m_hcScan(i, j).to_Mesh()).save(
                    to_save + "half_scan_" + std::to_string(i) + "_" + std::to_string(j) + ".stl");
        }
    auto distr = mes.computeHalfCoaptDistrib(500);
    mes.uniteCollidingCuspHalves();
    for (int i = 0; i < 3; ++i)
        Object3D(mes.m_coaptScan[i].to_Mesh()).save(to_save + "unite_scan_" + std::to_string(i) + ".stl");
    auto fdistr = mes.computeCoaptDistribution({(int)cr.cd[0].N - 1, (int)cr.cd[1].N - 1, (int)cr.cd[2].N - 1});
    Messurer::saveCoaptDistribCSV(fdistr, to_save + "full_distrib.csv");
    cr.coaptStatus = mes.computeCoaptStatus();
    cr.centralCoaptation = mes.computeHc();
    cr.maxCoaptationHeight = distr.getCoaptH();
    std::cout << "Hc = " << cr.centralCoaptation << " H = " << cr.maxCoaptationHeight << " Closed = " << ((cr.coaptStatus & (1 << 6)) > 0) <<   std::endl;
    auto bill = mes.computeBillowing(); bill.print();
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            cr.bill.billowPlate[i][k] = bill.plane[i][k];
            cr.bill.billowPnt[i][k] = bill.bilPnt[i][k];
        }
        cr.bill.billowing[i] = bill.getBillowing(i);
    }
    for (int i = 0; i < 3; ++i){
        std::copy(fdistr[i].up.begin(), fdistr[i].up.end(), cr.cd[i].top);
        std::copy(fdistr[i].down.begin(), fdistr[i].down.end(), cr.cd[i].bottom);
        cr.cd[i].width = fdistr[i].width;
    }
    Messurer::saveHalfCoaptDistribCSV(distr, to_save + "half_distrib.csv");
}

CoaptationResult createCoaptationResult(int N) {
    CoaptationResult res;
    for (int i = 0; i < 3; ++i){
        res.cd[i].N = N;
        res.cd[i].bottom = new float [N];
        res.cd[i].top = new float [N];
    }
    return res;
}

void freeCoaptationResult(CoaptationResult *cr) {
    for (int i = 0; i < 3; ++i){
        cr->cd[i].N = -1;
        if (cr->cd[i].bottom) delete[] cr->cd[i].bottom;
        if (cr->cd[i].top) delete[] cr->cd[i].top;
        cr->cd[i].bottom = nullptr;
        cr->cd[i].top = nullptr;
    }
}

int MeasureTest() {
    Object3D aorticSurface("../result/aorta_cut.stl");
    std::array<Object3D, 3> cusps;
    auto readTag = [](Object3D& obj, string filename, string tag){
        ifstream ob(filename, std::ios::binary);
        if (!ob) return false;
        auto m_x = obj.m_mesh.add_property_map<V_ind, Point>(tag);
        for (auto v: obj.m_mesh.vertices()){
            std::array<double, 3> P;
            ob.read(reinterpret_cast<char *> ((P.data())), 3 * sizeof(double));
            m_x.first[v] = Point(P[0], P[1], P[2]);
        }
        ob.close();
        return true;
    };
    for (int i = 0; i < 3; ++i){
        cusps[i].read("../result/leaf_res" + to_string(i) + ".txt");
        cusps[i].name = "leaf" + std::to_string(i);
        readTag(cusps[i], "../result/leaf_res" + to_string(i) + "_m_x0.tag", "v:point");
    }

    string to_save = "../result/";
    Messurer mes(aorticSurface, cusps[0], cusps[1], cusps[2], Vector{-0.706728, -0.271291, 0.653404});
    mes.setMargin(0.18);
    mes.computeCollidingBnd(4);
    std::cout << "Hc = " << mes.computeHc() << std::endl;
    mes.computeMidPlanes();

    for (int i = 0; i < 3; ++i) {
        auto colfilter = [&mes, &col = mes.m_colMap[i]](F_ind f) {
            auto lbl = col.face_lbl[col.remap[f]];
            bool res = true;
            for (int k = 0; k < 3; ++k)
                res &= ((lbl & (7 << 3 * k)) > 0);
            return res;
        };
        auto nocolfilter = [&colfilter](F_ind f) { return !colfilter(f); };
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_init.stl", "v:point");
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_init_col.stl", "v:point",
                            colfilter);
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_init_nocol.stl", "v:point",
                            nocolfilter);
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_real.stl", "v:x");
        mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_real_col.stl", "v:x", colfilter);
    }
    mes.computeHalfCoaptScans();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            Object3D(mes.m_hcScan(i, j).to_Mesh()).save(
                    to_save + "half_scan_" + std::to_string(i) + "_" + std::to_string(j) + ".stl");
        }
    auto distr = mes.computeHalfCoaptDistrib(500);
    mes.uniteCollidingCuspHalves();
    for (int i = 0; i < 3; ++i)
        Object3D(mes.m_coaptScan[i].to_Mesh()).save(to_save + "unite_scan_" + std::to_string(i) + ".stl");
    auto fdistr = mes.computeCoaptDistribution(1000);
    Messurer::saveCoaptDistribCSV(fdistr, to_save + "full_distrib.csv");
    auto coaptStatus = mes.computeCoaptStatus();
    auto centralCoaptation = mes.computeHc();
    auto maxCoaptationHeight = distr.getCoaptH();
    std::cout << "Hc = " << centralCoaptation << " H = " << maxCoaptationHeight << " Closed = " << ((coaptStatus & (1 << 6)) > 0) <<   std::endl;
    auto bill = mes.computeBillowing(); bill.print();
    Messurer::saveHalfCoaptDistribCSV(distr, to_save + "half_distrib.csv");

    return 0;
}

AVSimulator::SewLinesMeta AVSimulator::SewLinesMeta::Init(AVSimulator& m){
    SewLinesMeta s;
    s.commissure_points = m.sld.commissure_points;
    s.orientation = m.sld.orientation;

    typedef CGAL::Simple_cartesian<DReal>::Circle_3 Circle_3;
    Circle_3 c(s.commissure_points[0], s.commissure_points[1], s.commissure_points[2]);
    s.circle_center = c.center(); s.circle_R = sqrt(c.squared_radius());

    for (int l = 0; l < 3; ++l){
        auto& mes = s.meassures[l];
        mes.comissPnts[0] = CGAL::ORIGIN+m.sld.lines[l][0], mes.comissPnts[1] = CGAL::ORIGIN+m.sld.lines[l][m.sld.lines[l].size()-1];
        mes.linDist = sqrt((mes.comissPnts[1] - mes.comissPnts[0]).squared_length());
        auto low_it = std::min_element(m.sld.lines[l].begin(), m.sld.lines[l].end(), [n = s.orientation](auto& a, auto& b){ return n*a < n*b; });
        mes.lowPoint = CGAL::ORIGIN + (*low_it);
        mes.dist_low_plane = s.orientation * (s.commissure_points[0] - mes.lowPoint);
        mes.len = 0;
        for (int i = 1; i < m.sld.lines[l].size(); ++i)
            mes.len += sqrt((m.sld.lines[l][i] - m.sld.lines[l][i-1]).squared_length());
        //TODO: below simple not effective algorithm O(n^2), but exists algo O(n)
        typedef CGAL::Simple_cartesian<DReal>::Plane_3 Plane_3;
        auto findMostDistant = [](Vector n, std::vector<Vector>& crv){
            double dist = 0;
            std::array<Vector, 2> res;
            std::vector<Vector> intersect;
            for (int i = 0; i < crv.size(); ++i) {
                intersect.resize(0);
                for (int j = crv.size()-1; j > 0; --j){
                    if (((crv[j] - crv[i]) * n) * ((crv[j-1] - crv[i]) * n) < 0){
                        Vector v = crv[j-1] + (crv[j] - crv[j-1]) * ((crv[i] - crv[j-1])*n) / ((crv[j] - crv[j-1])*n);
                        intersect.push_back(v);
                    }
                }
                for (auto& v: intersect){
                    auto d = sqrt((crv[i] - v).squared_length());
                    if (d > dist)
                        dist = d, res[0] = crv[i], res[1] = v;
                }
            }
            return std::tuple<double, std::array<Vector, 2>>{dist, res};
        };
        Point p0 = mes.comissPnts[0], p1 = mes.comissPnts[1];
        auto maxDist1 = findMostDistant(Plane_3(p0,p1,s.commissure_points[(l+2)%3]).orthogonal_vector(),
                                        m.sld.lines[l]);
        mes.maxDist1 = get<0>(maxDist1);
        mes.pnts1 = {CGAL::ORIGIN + get<1>(maxDist1)[0], CGAL::ORIGIN + get<1>(maxDist1)[1]};

        Plane_3 low_plane(p0, p1, mes.lowPoint);
        Plane_3 plane2(p0, p1, p0 + low_plane.orthogonal_vector());
        auto maxDist2 = findMostDistant(plane2.orthogonal_vector(),m.sld.lines[l]);
        mes.maxDist2 = get<0>(maxDist2);
        mes.pnts2 = {CGAL::ORIGIN + get<1>(maxDist2)[0], CGAL::ORIGIN + get<1>(maxDist2)[1]};
    }

    return s;
}

bool AVSimulator::SewLinesMeta::prepareOzakiLine(const Object3D& aorta) {
    typedef CGAL::Simple_cartesian<DReal>::Plane_3 Plane_3;
    for (int l = 0; l < 3; ++l){
        auto& mes = meassures[l];
        Plane_3 ozaki_plane(mes.comissPnts[0], mes.comissPnts[1], mes.lowPoint);
        Vector orient = ozaki_plane.projection(mes.comissPnts[0] + orientation) - mes.comissPnts[0];
        orient /= sqrt(orient.squared_length());

        typedef CGAL::Simple_cartesian<DReal>                     K;
        typedef std::list<Polyline_type>                          Polylines;
        namespace PMP = CGAL::Polygon_mesh_processing;

        CGAL::Polygon_mesh_slicer<Mesh, K> slicer(aorta.m_mesh);
        Polylines polylines;
        slicer(ozaki_plane, std::back_inserter(polylines));
        auto absn = [](Vector a){ return abs(*std::max_element(a.cartesian_begin(), a.cartesian_end(), [](DReal a, DReal b){ return abs(a) < abs(b); })); };

        typedef CGAL::Simple_cartesian<DReal>::Segment_3 Segment;
        typedef std::vector<Segment>::iterator Iterator;
        typedef CGAL::AABB_segment_primitive<K, Iterator> Primitive;
        typedef CGAL::AABB_traits<K, Primitive> Traits;
        typedef CGAL::AABB_tree<Traits> Tree;

        std::pair<Point, int> p1, p2;
        double d1 = DBL_MAX, d2 = DBL_MAX;
        bool choosable = true;
        auto itl = polylines.end();
        for(auto it = polylines.begin(); it != polylines.end(); ++it){
            auto& line = *it;
            std::vector<Segment> lsegments;
            for (int i = 0; i < line.size()-1; ++i) lsegments.push_back(Segment{line[i], line[i+1]});
            Tree tree(lsegments.begin(),lsegments.end());
            auto lp1 = tree.closest_point_and_primitive(mes.comissPnts[0]);
            auto lp2 = tree.closest_point_and_primitive(mes.comissPnts[1]);
            double ld1 = CGAL::squared_distance(lp1.first, mes.comissPnts[0]);
            double ld2 = CGAL::squared_distance(lp2.first, mes.comissPnts[1]);
            if (ld1 < d1 && ld2 < d2){
                p1 = {lp1.first, lp1.second - lsegments.begin()}, p2 = {lp2.first, lp2.second - lsegments.begin()}, d1 = ld1, d2 = ld2, itl = it;
            } else if (ld1 < d1 || ld2 < d2){
                choosable = false;
                std::cout << "Warrning: Can't choose real slice\n";
            }
        }
        double err = 1e-5*(1 + absn(mes.comissPnts[0] - CGAL::ORIGIN));
        if (d1 > err || d2 > err || !choosable) continue;
        if (p1.second > p2.second) std::swap(p1, p2);

        auto& line = *itl;
        Vector mid = CGAL::NULL_VECTOR;
        double len = 0;
        Segment st(p1.first, line[p1.second+1]);
        double stl = sqrt(st.squared_length());
        mid += ((p1.first - CGAL::ORIGIN) + (line[p1.second+1] - CGAL::ORIGIN))/2 * stl;
        Segment end(line[p2.second], p2.first);
        double endl = sqrt(st.squared_length());
        mid += ((p2.first - CGAL::ORIGIN) + (line[p2.second] - CGAL::ORIGIN))/2 * endl;
        len +=  stl+endl;

        for (int i = p1.second+1; i < p2.second; ++i){
            double llen = sqrt((line[i+1] - line[i]).squared_length());
            len += llen;
            mid += ((line[i+1] - CGAL::ORIGIN) + (line[i] - CGAL::ORIGIN))/2 * llen;
        }
        mid /= len;
        if ((mes.comissPnts[0] - (CGAL::ORIGIN + mid)) * orient > 0) {
            mes.ozakiLineLen = len;
            mes.ozakiLine.reserve(p2.second - p1.second + 1);
            mes.ozakiLine.push_back(p1.first);
            for (int i = p1.second+1; i <= p2.second; ++i) mes.ozakiLine.push_back(line[i]);
            mes.ozakiLine.push_back(p2.first);
        } else {
            st = Segment(p2.first, line[(p2.second+1)%line.size()]);
            end = Segment(line[p1.second], p1.first);
            stl = sqrt(st.squared_length());
            endl = sqrt(st.squared_length());
            mid = CGAL::NULL_VECTOR;
            len = stl + endl;
            mid += ((p2.first - CGAL::ORIGIN) + (line[(p2.second+1)%line.size()] - CGAL::ORIGIN))/2 * stl;
            mid += ((line[p1.second] - CGAL::ORIGIN) + (p1.first - CGAL::ORIGIN))/2 * endl;
            for (int i = p2.second+1; i <= line.size() + p1.second; ++i){
                int j = i%line.size(), jp1 = (i+1)%line.size();
                double llen = sqrt((line[jp1] - line[j]).squared_length());
                len += llen;
                mid += ((line[jp1] - CGAL::ORIGIN) + (line[j] - CGAL::ORIGIN))/2 * llen;
            }
            mid /= len;
            if ((mes.comissPnts[0] - (CGAL::ORIGIN + mid)) * orient <= 0)
                std::cout << "Warrning: slicing has errors" << std::endl;
            mes.ozakiLineLen = len;
            mes.ozakiLine.reserve(line.size() - (p2.second - p1.second));
            mes.ozakiLine.push_back(p2.first);
            for (int i = p2.second+1; i < line.size() - 1; ++i) mes.ozakiLine.push_back(line[i]);
            for (int i = 0; i <= p1.second; ++i) mes.ozakiLine.push_back(line[i]);
            mes.ozakiLine.push_back(p1.first);
        }
        //save_polyline_csv(mes.ozakiLine, "../result/PigTest/NCC_poly6_oz" + to_string(l) + ".csv");
    }

    return (meassures[0].ozakiLine.size() != 0) && (meassures[1].ozakiLine.size() != 0) && (meassures[2].ozakiLine.size() != 0);
}

void AVSimulator::SewLinesMeta::save_polyline_csv(AVSimulator::SewLinesMeta::Polyline_type &line, string fname) {
    std::ofstream f(fname);
    int i = 0;
    for (auto& p: line)
        f << p[0] << ", " << p[1] << ", " << p[2] << ", " << i++ << "\n";
    f.close();
}



