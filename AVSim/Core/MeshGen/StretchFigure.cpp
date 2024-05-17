//
// Created by alex on 04.05.2021.
//

#include "StretchFigure.h"
#include "Helper.h"

Ani3dSurfDiscrWrap StretchFigure::prepareFigure() const {
    double st = m_tPnts.front(), end = m_tPnts.back();
    bool st_tw = m_dirLen(st) > 0, end_tw = m_dirLen(end) > 0;
    Ani3dSurfDiscrWrap asd;
    auto p = asd.addPoint(Point{-1000, -1000, -1000});
    auto p0 = p;
    vector<LineID> lines;
    for (int i = 1; i < m_tPnts.size(); ++i){
        auto p1 = asd.addPoint(Point{1000, 1000, 1000});
        AniSurfLine asl;
        asl.st = p, asl.end = p1;
        asl.use_func = true;
        asl.func_param = 1;
        asl.lbl = m_eLbls[i];
        asl.is_reverse = false;
        asl.st_t = m_tPnts[i-1], asl.end_t = m_tPnts[i];
        lines.push_back(asd.addLine(asl));
        p = p1;
    }
    if (end_tw){
        auto p1 = asd.addPoint();
        AniSurfLine asl;
        asl.st = p, asl.end = p1;
        asl.use_func = true;
        asl.func_param = 3;
        asl.lbl = m_eLbls[m_tPnts.size()];
        asl.st_t = 0, asl.end_t = 1;
        lines.push_back(asd.addLine(asl));
        p = p1;
    }
    for (int i = static_cast<int>(m_tPnts.size()) - 2; i > 0; --i){
        auto p1 = asd.addPoint();
        AniSurfLine asl;
        asl.st = p, asl.end = p1;
        asl.func_param = 2;
        asl.use_func = true;
        asl.lbl = m_eLbls[2*m_tPnts.size() - i - 1]; //m_eLbls[2*m_tPnts.size() - i]
        asl.is_reverse = false;
        asl.st_t = m_tPnts[i+1], asl.end_t = m_tPnts[i];
        lines.push_back(asd.addLine(asl));
        p = p1;
    }
    if (st_tw) {
        auto p1 = asd.addPoint();
        AniSurfLine asl;
        asl.st = p, asl.end = p1;
        asl.func_param = 2;
        asl.use_func = true;
        asl.lbl = m_eLbls[2*m_tPnts.size() - 1];
        asl.is_reverse = false;
        asl.st_t = m_tPnts[1], asl.end_t = m_tPnts[0];
        lines.push_back(asd.addLine(asl));
        p = p1;

        AniSurfLine asl1;
        asl1.st = p, asl1.end = p0;
        asl1.use_func = true;
        asl1.func_param = 4;
        asl1.lbl = m_eLbls[0];
        asl1.st_t = 1, asl1.end_t = 0;
        lines.push_back(asd.addLine(asl1));
        p = p0;
    }
    else {
        auto p1 = p0;
        AniSurfLine asl;
        asl.st = p, asl.end = p1;
        asl.func_param = 2;
        asl.use_func = true;
        asl.lbl = m_eLbls[2*m_tPnts.size() - 1];
        asl.is_reverse = false;
        asl.st_t = m_tPnts[1], asl.end_t = m_tPnts[0];
        lines.push_back(asd.addLine(asl));
        p = p1;
    }
    AniSurfSurf as;
    as.lines = std::move(lines);
    as.use_func = true;
    as.is_positive_orient = true;
    as.func_param = 0;
    as.st_u = st, as.end_u = end;
    as.st_v = 0; as.end_v = 1;
    asd.addSurface(as);
    asd.set_LineParamFcn([st, end](int param, double t, double* pu, double* pv)->int{
        switch(param){
            case 1: { *pu = t  , *pv = 0; break; }
            case 2: { *pu = t  , *pv = 1; break; }
            case 3: { *pu = end, *pv = t; break; }
            case 4: { *pu = st , *pv = t; break; }
        }
        return 1;
    });
    asd.set_SurfParamFcn([dir = m_dir, curve = m_curve, len = m_dirLen](int prm, double u, double v, double* px, double* py, double* pz) -> int{
        Point p = curve(u) + len(u) * v * dir;
        *px = p.x(), *py = p.y(), *pz = p.z();
        return 1;
    });


    return asd;
}

std::function<Point_2(double, double)> StretchFigure::makeParametricFlatter() const{
    Point p0 = m_curve(m_tPnts.back());
    using Plane_3 = CGAL::Simple_cartesian<double>::Plane_3;
    Plane_3 plane(p0, m_dir);
    auto projection = [plane = Plane_3(p0, m_dir), curve = m_curve](double t) { return plane.projection(curve(t)); };
    PiecewieceLineIntegral pli;
    pli.pieces.reserve(m_tPnts.size() - 1);
    for (int i = 1; i < m_tPnts.size(); ++i){
        QuickLineIntegral qli;
        qli.m_curve = projection;
        qli.Tst = m_tPnts[i-1], qli.Tend = m_tPnts[i];
        qli.rel = 1e-7;
        qli.Ninit = 20;
        pli.pieces.push_back(std::move(qli));
    }
    pli.setup();
    double L = pli.pieces.back().quick_index.back().l;
    auto flat = [p0, len = m_dirLen, curve = m_curve, dir = m_dir, pli = std::move(pli), L](double u, double v) -> Point_2{
        double y = (curve(u) - p0) * dir + len(u) * v;
        double x = pli(u) - L / 2;
        return Point_2(x, y);
    };
    return flat;
}

Object3D StretchFigure::generateMesh(std::function<double(double, double, double)> fsize){
    auto asd = prepareFigure();
    AniMesh am;
    Ani3dMeshOut amo(am);
    std::vector<double> v_u;
    amo.vu_map = &v_u;
    makeAft3dSurface(asd.getAniSurfDiscr(), fsize, amo);
    mark_vertices(am);
    invertFaceOrientation(am);
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj(m);
    auto flatter = makeParametricFlatter();
    for (auto v: m.vertices()){
        Point_2 p2 = flatter(v_u[2*v.idx()+0], v_u[2*v.idx()+1]);
        obj.m_x0[v] = Point(p2.x(), p2.y(), 0);
    }
    return obj;
}

StretchFigure &StretchFigure::setPntParams(const vector<double> &t_pnts) {
    m_eLbls.resize(t_pnts.size()*2);
    for (int i = 0; i <= t_pnts.size(); ++i) m_eLbls[i] = 2;
    for (unsigned i = 1+t_pnts.size(); i < 2*t_pnts.size(); ++i) m_eLbls[i] = 1;
    return m_tPnts = t_pnts, *this;
}

StretchFigure &StretchFigure::setPntParams(const vector<double> &t_pnts, const vector<int> &eLbls) {
    m_eLbls.resize(t_pnts.size()*2);
    for (int i = 0; i <= t_pnts.size(); ++i) m_eLbls[i] = 2;
    for (unsigned i = 1+t_pnts.size(); i < 2*t_pnts.size(); ++i) m_eLbls[i] = 1;
    for (unsigned i = 0, cnt = std::min(eLbls.size(), m_eLbls.size()); i < cnt; ++i) m_eLbls[i] = eLbls[i];
    return m_tPnts = t_pnts, *this;
}

void QuickLineIntegral::setup() {
    struct _IntegrData{
        Point p;
        double t;
        double dl;
        bool isRedusedEnough;
    };
    std::list<_IntegrData> integr;
    {
        double dt = (Tend - Tst) / Ninit, T = Tst;
        Point prev = m_curve(T);
        integr.push_back(_IntegrData{prev, T, 0, true});
        for (int i = 1; i <= Ninit; ++i) {
            T += dt;
            Point next = m_curve(T);
            Vector dlt = next - prev;
#if __cplusplus < 201703L
            auto hypot = [](auto x, auto y, auto z){ return sqrt(x*x+y*y+z*z); };
#endif
            double dl = hypot(dlt.x(), dlt.y(), dlt.z());
            integr.push_back(_IntegrData{next,T, dl, false});
            prev = next;
        }
    }
    bool is_changed = true;
    int ni = 0;
    for (ni = 0; ni < max_refine_depth && is_changed; ++ni){
        bool _is_changed = false;
        auto itl = integr.begin();
        auto end = integr.end();
        auto itr = itl;
        ++itr;
        while (itr != end) {
            if (!itr->isRedusedEnough) {
                double tm = (itl->t + itr->t) / 2;
                Point l = itl->p;
                Point r = itr->p;
                Point m = m_curve(tm);
#if __cplusplus < 201703L
                auto hypot = [](auto x, auto y, auto z){ return sqrt(x*x+y*y+z*z); };
#endif
                double lm = hypot(m.x() - l.x(), m.y() - l.y(), m.z() - l.z());
                double mr = hypot(r.x() - m.x(), r.y() - m.y(), r.z() - m.z());
                if (lm + mr - itr->dl > rel * itr->dl) {
                    integr.insert(itr, _IntegrData{m, tm, lm, false});
                    itr->dl = mr;
                    _is_changed = true;
                } else {
                    itr->isRedusedEnough = true;
                }
                itl = itr;
                itr++;
            } else {
                itl = itr;
                itr++;
            }
        }
        is_changed = _is_changed;
    }
    quick_index.reserve(integr.size());
    double len = 0;
    for (auto it: integr){
        len += it.dl;
        quick_index.push_back(QuickIntegrData{it.t, len, it.p});
    }
}

double QuickLineIntegral::operator()(double t) const {
    if (quick_index.size() < 2) return -1;
    if (t >= Tend) return quick_index.back().l;
    if (t <= Tst) return quick_index.front().l;
    int id;
    {
        int a = 0, b = quick_index.size();
        while (b - a > 1){
            int m = (b + a) / 2;
            if (quick_index[m].t > t){
                b = m;
            } else if (quick_index[m].t == t){
                a = b = m; b++;
            } else {
                a = m;
            }
        }
        id = a;
    }
    Vector dlt = m_curve(t) - quick_index[id].p;
#if __cplusplus < 201703L
    auto hypot = [](auto x, auto y, auto z){ return sqrt(x*x+y*y+z*z); };
#endif
    double L = quick_index[id].l + hypot(dlt.x(), dlt.y(), dlt.z());
    return L;
}

int StretchFigureTest(string dir) {
    auto obj = StretchFigure()
            .setStretchDirection(Vector{0, 1, 1})
            .setStretchLengthFunc([](double t) { return 1 - abs(t); })
            .setCurveFunc([](double t) { return Point(t, t*t, 0); })
            .setPntParams({-1, 0, 1})
            .generateMesh(Ani3dFSize(0.1).fsize);
    obj.save(dir + "test_x0.vtk", obj.m_x0);
    obj.save(dir + "test_x.vtk");

    return 0;
}

double PiecewieceLineIntegral::operator()(double t) const {
    if (t <= pieces.front().Tst) return pieces.front().quick_index.front().l;
    if (t >= pieces.back().Tend) return pieces.back().quick_index.back().l;
    int id = -1;
    {
        int a = 0, b = pieces.size();
        while (b - a >= 1){
            int m = (b + a) / 2;

            if (t > pieces[m].Tend){
                a = m+1;
            } else if (t < pieces[m].Tst){
                b = m;
            } else if (t <= pieces[m].Tend && t >= pieces[m].Tst){
                id = m; break;
            }
        }
        if (id < 0) throw std::runtime_error("Pieces is not contiguous");
    }
    return pieces[id](t);
}

void PiecewieceLineIntegral::setup() {
    for (auto& p: pieces) p.setup();
    for (int i = 1; i < pieces.size(); ++i){
        if (pieces[i].Tst != pieces[i-1].Tend) throw std::runtime_error("Pieces is not contiguous");
        double dl = pieces[i-1].quick_index.back().l;
        for (auto& q: pieces[i].quick_index) q.l += dl;
    }
}
