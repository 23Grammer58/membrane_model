//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"
#include "TemplatesCollection.h"
#include "Elliptic.h"
#include <Eigen/Dense>

struct SurfPrms{
    double a;   //elliptical big semiaxis
    double b;   //elliptical small semiaxis
    double d;   //shift of cutting through a-axis counting from center
    double dc;  //shift of cutting central plates
    double dz;
    double phi[3]; //angles of cutting
    double mesh_h;
};

static int _check_initial_intersection(SurfPrms p)
{
    auto non_inter_len = p.dz - (p.a + p.d) * tan (p.phi[0]);
    if (non_inter_len < p.mesh_h/4)
    {
        auto lost_len = non_inter_len;
        if (lost_len > 0){
            std::cout << "Warning: Length of non intersected part is " <<  lost_len
            << " will be ignored because too big mesh step = " << p.mesh_h
            << ". Please consider mesh step less " << lost_len * 4 << std::endl;
        }
        return 0;   // вся граница определяется пересечением
    }
    //если не так, то пересечения может не быть
    double r = p.b * sqrt(1 - (p.d / p.a) * (p.d / p.a));
    Eigen::Vector3d p0(- p.d - p.a, -r, -p.dz), p1(cos(p.phi[2]), 0, sin(p.phi[2]));
    auto n = p1.cross(p0);
    Eigen::AngleAxis<double> q(p.phi[0], Eigen::Vector3d(0, 1, 0));
    n = q*n;
    Eigen::Vector2d pp(n[1], -n[0]);
    double k1 = pp[1] / pp[0];
    double k2 = -(p.d/r) * (p.b / p.a) * (p.b / p.a) * cos (p.phi[0]);
    if (k1 < k2) return 2;                  //граница определяется начальной плоскостью
    return 1;                               //граница кусочная
}

struct IntersectionParam{
    double k, a, d;
    double l0, l1, l00;
    double phi[2];
    Eigen::Vector3d pa, pb, pdif_z;
    double inter_phi[2];
};

static IntersectionParam init_bound_params(SurfPrms p, int mode)
{
    IntersectionParam res;

    double k = sqrt(1 - (p.b / p.a) * (p.b / p.a));
    double a = p.a;
    double d = p.d;
    double l1 = p.a * (-gsl_sf_ellint_E(M_PI_2, k, GSL_PREC_DOUBLE));
    double* pphi = p.phi;

    res.k = k;
    res.a = a;
    res.d = d;
    res.l1 = l1;
    res.phi[0] = pphi[0];
    res.phi[1] = pphi[1];

    double r00 = asin(p.d / p.a);
    double l00 = p.a * (gsl_sf_ellint_E(r00, k, GSL_PREC_DOUBLE));
    res.l00 = l00;

    if (mode == 2)
    {
        res.l0 = l00;
        return res;
    }

    double r = p.b * sqrt(1 - (p.d / p.a) * (p.d / p.a));
    Eigen::Vector3d f0(- p.d - p.a, -r, -p.dz), f1(cos(p.phi[2]), 0, sin(p.phi[2]));
    Eigen::Vector3d n = f1.cross(f0).normalized();
    Eigen::Vector3d pa, pb;
    solve_plane_to_cylinder_ellips(p.a, p.b, n.data(), pa.data(), pb.data());

    Eigen::Vector3d p0(-p.d, -r, 0), p1;
    if (mode == 0)
    {
        p1 << p.a, 0, p.dz;
        res.l0 = l1;
    }
    else if (mode == 1)
    {
        Eigen::AngleAxis<double> q(p.phi[0], Eigen::Vector3d(0, 1, 0));
        n = q*n;
        double d = n[0] * p0[0] / cos(p.phi[0]) + n[1]*p0[1];
        double loca = p.a / cos (p.phi[0]), locb = p.b;
        Eigen::Vector2d intersect;
        solve_2d_line_to_ellips(n[0], n[1], d, loca, locb, 2, intersect.data());
        Eigen::Vector3d newp(intersect[0] + p.d / cos(p.phi[0]), intersect[1], 0);
        newp = q.inverse()*newp;
        n = q.inverse()*n;
        newp[0] -= p.d;
        p1 = newp;

        double r0 = asin(-p1[0] / p.a);
        double l0 = p.a * (gsl_sf_ellint_E(r0, k, GSL_PREC_DOUBLE));
        res.l0 = l0;
    }

    double difz = n.dot(p0) / n[2];
    Eigen::Vector3d pdif_z(0, 0, difz);
    p0[2] -= difz;
    p1[2] -= difz;
    double loca = pa.norm(), locb = pb.norm();
    Eigen::Vector2d p2[2] = { {p0.dot(pa)/loca, p0.dot(pb)/locb}, {p1.dot(pa)/loca, p1.dot(pb)/locb} };
    double phi[2];
    for (int i = 0; i < 2; ++i)
    {
        double sinp = p2[i][0] / loca, cosp = p2[i][1] / locb;
        phi[i] = acos(cosp);
        if (sinp < 0) phi[i] = 2 * M_PI - phi[i];
    }


    res.pa = pa;
    res.pb = pb;
    res.pdif_z = pdif_z;
    res.inter_phi[0] = phi[0];
    res.inter_phi[1] = phi[1];

    return res;
}

static SurfPrms computeSurfParams(ThubricarValvePrms vp, double size){
    double R1 = vp.R1;
    double kR = vp.R2 / vp.R1;
    double kH = vp.H / vp.R1;
    double kdh = vp.dh / vp.R1;
    double kr = cos(M_PI / vp.nleafs);
    double phi = atan2(kH, (kR - kr));
    double alpha = vp.alpha;
    double omega = phi - alpha;
    double teta = atan2(kdh, ((kR - kr) + (kH - kdh) * tan(alpha)));
    double psi = alpha - teta;
    double ka =(1 + (1 + kr) / (2 * kR - 1 - kr)) * cos (omega) / sin (phi) * kH / 2;
    double kb = sqrt((1 - kr) / (2 * kR - 1 - kr) ) * kR;
    if (ka < kb) ka = kb;
    double kd = cos(omega) / sin(phi) * (kR - 1 - kr) / (2*kR - 1 - kr) * kH;
    double kdc = kd + kr * (1 + (tan(phi) - tan(alpha)) / (1.0/tan(2 * phi) + tan(alpha))) / cos(alpha);
    double kdz = (kR * cos(2 * phi) + (kH + kr * tan(phi)) * sin(2*phi)) / sin(2 * phi - alpha) - kH * sin(phi - alpha) / sin (phi);
    double gamma = M_PI_2 + alpha - 2 * phi;

    if (ka < kb) std::cout <<"Warning: a = " << ka*R1 << " < " << kb*R1 <<" = b" << std::endl;
    SurfPrms p{ka * R1, kb * R1, kd * R1, kdc * R1, kdz * R1, {psi, -omega, gamma}, size};
    return p;
}

static int createThubrikarFlatScan(AniMesh& am, ThubricarValvePrms vp, double size){
    auto p = computeSurfParams(vp, size);
    double k = sqrt(1 - (p.b / p.a) * (p.b / p.a));
    inv_ellint_E_init(k);
    int mode = _check_initial_intersection(p);
//    switch (mode) {
//        case 0: std::cout << "All border defined by symmetry intersection" << std::endl; break;
//        case 1: std::cout << "The border has intersections, but not entirely"<< std::endl; break;
//        case 2: std::cout << "Borders doesn't intersects" << std::endl; break;
//    }
    if (mode < 0 || mode > 2) throw std::runtime_error("Error: Check_initial_intersection fails");
    auto inp = init_bound_params(p, mode);

    auto get_flatted_bound = [&p = inp](const double t, int id)  -> std::array<double, 2>{           //t in [-1, 1]
        double l = p.l1 + (((id==1) ? p.l00 : p.l0) - p.l1) * t;
        double z = (p.d - p.a * sin(inv_ellint_E(l / p.a, p.k))) * tan(p.phi[id]);
        return std::array<double, 2>{l - p.l1, z}; //центрированный интервал
    };
    auto get_intersected_bound = [&p = inp](const double t, int id) -> std::array<double, 2>{    //t in [0, 1] (id = 0) and [1, 0] (id = 1)
        auto phi = p.inter_phi;
        double angle = phi[0] + (phi[1] - phi[0]) * t;
        Eigen::Vector3d pp = sin(angle)*p.pa + cos(angle)*p.pb + p.pdif_z;
        double rel = -pp[0] / p.a;
        if (fabs(rel) > 1) rel = (rel > 0) ? 1 : -1;
        double glob_angle = asin(rel);

        double l = p.a * (gsl_sf_ellint_E(glob_angle, p.k, GSL_PREC_DOUBLE));
        return std::array<double, 2>{((id == 1) ? -1 : 1) * (l - p.l1), pp[2]}; //центрированный интервал
    };
    auto CRVParametricFunc = [get_flatted_bound, get_intersected_bound, &inp](int i, double t, double *pu, double *pv) -> int{
        std::array<double, 2> p = {0, 0};
        switch (i)
        {
            case 0: std::cout << "something wrong, i = 0" << std::endl;
            case 1:
            case 2: p = get_flatted_bound(t - 1, i - 1); break;
            case 3:
            case 4: p = get_intersected_bound(t, i - 3); break;
            default: return -1;
        }

        return *pu = p[0], *pv = p[1], 1;
    };
    auto FlatSurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1;};
    int status = -1;
    switch (mode) {
        case 0: {
            const int nVVert = 4, nLine = 5, nSurface = 2;
            std::array<std::array<double, 2>, nVVert> pp{get_flatted_bound(1, 1),
                                                        get_flatted_bound(0, 1),
                                                        get_flatted_bound(-1, 1),
                                                        get_intersected_bound(1, 0)};
            double VVert[nVVert*3] = {pp[0][0],pp[0][1],0,  pp[1][0],pp[1][1],0,  pp[2][0],pp[2][1],0,  pp[3][0],pp[3][0],0};
            int LineD[nLine*3] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 1, 1,  2, 4, 1};
            int LineP[nLine*2] = {1, 2,  1, 2,  1, 4,  1, 3,  0, 0};
            int exportCurves[nLine] = {2, 2, 1, 4, 0};
            double LineT[nLine*2] = {2, 1,  1, 0,   0, 1,  1, 0,  0, 0};
            int SurfL[nSurface*5] = {3, 0, 1, 0, 0,
                                     3, 0, 1, 0, 0};
            int SurfI[] = {1, 0,  5, 0,  4, 0,
                           2, 0,  3, 0,  5, 1};
            double y_max = std::max(fabs(pp[1][1]), fabs(pp[3][1]));
            double x_max = std::max(fabs(pp[0][0]), fabs(pp[2][0]));
            double SurfT[nSurface*4] = {-x_max*2,-y_max*2, +x_max*2,+y_max*2,
                                        -x_max*2,-y_max*2, +x_max*2,+y_max*2};
            Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
            asd.bounline = CRVParametricFunc; asd.bounsurf = FlatSurf;

            status = makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
            break;
        }
        case 1: {
            const int nVVert = 5, nLine = 5, nSurface = 1;
            std::array<std::array<double, 2>, nVVert> pp{get_flatted_bound(1, 1),
                                                         get_flatted_bound(0, 1),
                                                         get_flatted_bound(-1, 1),
                                                         get_intersected_bound(1, 1),
                                                         get_intersected_bound(1, 0)};
            double VVert[nVVert*3] = {pp[0][0],pp[0][1],0,  pp[1][0],pp[1][1],0,  pp[2][0],pp[2][1],0,  pp[3][0],pp[3][1],0,  pp[4][0],pp[4][1],0};
            int LineD[nLine*3] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 5, 1,  5, 1, 1};
            int LineP[nLine*2] = {1, 2,  1, 2,  1, 4,  1, 1,  1, 3};
            int exportCurves[nLine] = {2, 2, 1, 8, 4};
            double LineT[nLine*2] = {2, 1,  1, 0,   0, 1,  0, 2,  1, 0};
            int SurfL[nSurface*5] = {5, 0, 1, 0, 0};
            int SurfI[] = {1, 0,  2, 0,  3, 0,  4, 0,  5, 0};
            double y_max = std::max(fabs(pp[1][1]), fabs(get_flatted_bound(0, 0)[1]));
            double x_max = std::max(fabs(pp[0][0]), fabs(pp[2][0]));
            double SurfT[nSurface*4] = {-x_max*2,-y_max*2, +x_max*2,+y_max*2};
            Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
            asd.bounline = CRVParametricFunc; asd.bounsurf = FlatSurf;

            status = makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
            break;
        }
        case 2: {
            const int nVVert = 4, nLine = 5, nSurface = 2;
            std::array<std::array<double, 2>, nVVert> pp{get_flatted_bound(1, 1),
                                                         get_flatted_bound(0, 1),
                                                         get_flatted_bound(-1, 1),
                                                         get_flatted_bound(0, 0)};
            double VVert[nVVert*3] = {pp[0][0],pp[0][1],0,  pp[1][0],pp[1][1],0,  pp[2][0],pp[2][1],0,  pp[3][0],pp[3][0],0};
            int LineD[nLine*3] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 1, 1,  2, 4, 1};
            int LineP[nLine*2] = {1, 2,  1, 2,  1, 1,  1, 1,  0, 0};
            int exportCurves[nLine] = {2, 2, 1, 4, 0};
            double LineT[nLine*2] = {2, 1,  1, 0,   0, 1,  1, 2, 0, 0};
            int SurfL[nSurface*5] = {3, 0, 1, 0, 0,
                                     3, 0, 1, 0, 0};
            int SurfI[] = {1, 0,  5, 0,  4, 0,
                           2, 0,  3, 0,  5, 1};
            double y_max = std::max(fabs(pp[1][1]), fabs(pp[3][1]));
            double x_max = std::max(fabs(pp[0][0]), fabs(pp[2][0]));
            double SurfT[nSurface*4] = {-x_max*2,-y_max*2, +x_max*2,+y_max*2,
                                        -x_max*2,-y_max*2, +x_max*2,+y_max*2};
            Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
            asd.bounline = CRVParametricFunc; asd.bounsurf = FlatSurf;

            status = makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
            break;
        }
    }
    inv_ellint_E_free();

    return status;
}

struct ConvParams{
    double l1, k, a, b, d, r, alpha, R1, shift;
    ConvParams(ThubricarValvePrms vp){
        R1 = vp.R1;
        double kR = vp.R2 / vp.R1;
        double kH = vp.H / vp.R1;
        double kdh = vp.dh / vp.R1;
        double kr = cos(M_PI / vp.nleafs);
        double phi = atan2(kH, (kR - kr));
        alpha = vp.alpha;
        double omega = phi - alpha;
        double teta = atan2(kdh, ((kR - kr) + (kH - kdh) * tan(alpha)));
        double psi = alpha - teta;
        double ka =(1 + (1 + kr) / (2 * kR - 1 - kr)) * cos (omega) / sin (phi) * kH / 2;
        double kb = sqrt((1 - kr) / (2 * kR - 1 - kr) ) * kR;
        if (ka < kb) ka = kb;
        double kd = cos(omega) / sin(phi) * (kR - 1 - kr) / (2*kR - 1 - kr) * kH;
        a = ka * R1, b = kb * R1;
        k = sqrt(1 - (kb / ka) * (kb / ka));
        l1 = a * (-gsl_sf_ellint_E(M_PI_2, k, GSL_PREC_DOUBLE));
        r = kr * R1;
        d = kd * R1;
        shift = (vp.ht / R1 < 1.0e-7) ?  1e-7 * R1 : vp.ht / (2 * sin(M_PI / vp.nleafs));
    }
};

template<class PointType>
std::function<PointType(PointType)> make_converter(ThubricarValvePrms vp, int ileaf = 0){
    ConvParams cp(vp);
    Eigen::AngleAxis<double> q(M_PI - cp.alpha, Eigen::Vector3d(0, 1, 0));
    Eigen::AngleAxis<double> qleaf(2 * M_PI / vp.nleafs * ileaf, Eigen::Vector3d(0, 0, 1));
    return [cp, ql = qleaf, iei = InvEllInt(cp.k), q](PointType p) -> PointType{
        Eigen::Vector3d rv;
        rv[2] = p[1];
        double phi = iei.compute((p[0] + cp.l1) / cp.a);
        rv[0] = cp.a * sin(phi);
        rv[1] = cp.b * cos(phi);
        rv[0] -= cp.d;
        rv = q * rv;
        rv[0] += cp.r;
        rv[0] += cp.shift; //small shift to guarantee that leafs is not intersect
        rv = ql * rv;
        return PointType(rv[0], rv[1], rv[2]);
    };
}

static void convertFromFlatTo3D(Object3D& am, ThubricarValvePrms vp){
    ConvParams cp(vp);
}

static void convertFromFlatTo3D(AniMesh& am, ThubricarValvePrms vp){
    ConvParams cp(vp);

    Eigen::AngleAxis<double> q(M_PI - cp.alpha, Eigen::Vector3d(0, 1, 0));

    inv_ellint_E_init(cp.k);
    for (int i = 0; i < am.vertices.size()/3; ++i){
        auto p = am.vertices.data() + 3*i;
        p[2] = p[1];
        double phi = inv_ellint_E((p[0] + cp.l1) / cp.a, cp.k);
        p[0] = cp.a * sin(phi);
        p[1] = cp.b * cos(phi);
        p[0] -= cp.d;
        auto pv = Eigen::Map<Eigen::Vector3d>(p);
        pv = q * pv;
        p[0] += cp.r;
    }
    inv_ellint_E_free();

    return;
}

AniMesh generateThubricarFlatLeaf(ThubricarValvePrms tvp, double size){
    if (tvp.nleafs < 3) throw std::runtime_error("Thubricar valve defined only for at least 3 cusps");
    AniMesh am;
    createThubrikarFlatScan(am, tvp, size);
    mark_vertices(am);
    invertFaceOrientation(am);

    return am;
}

std::vector<Object3D> generateThubricarOpenValve(ThubricarValvePrms tvp, double size){
    if (tvp.nleafs < 3) throw std::runtime_error("Thubricar valve defined only for at least 3 cusps");
    AniMesh am = generateThubricarFlatLeaf(tvp, size);
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    std::vector<Object3D> valve(tvp.nleafs);
    for (int i = 0; i < tvp.nleafs; ++i){
        Object3D& obj = (valve[i] = Object3D(m));
        auto conv = make_converter<Point>(tvp, i);
        std::for_each(obj.m_mesh.vertices().begin(), obj.m_mesh.vertices().end(), [&obj, &conv](auto v){ obj.m_next_x[v] = obj.m_x[v] = conv(obj.m_x0[v]); });
    }

    return valve;
}

std::vector<Object3D> generateThubricarCloseValve(ThubricarValvePrms tvp, double size){
    auto res = generateThubricarOpenValve(tvp, size);
    using Plane_3 = CGAL::Simple_cartesian<double>::Plane_3;
    for (auto& obj: res){
        if (obj.m_mesh.num_vertices() < 3) continue;
        Point p0 = obj.m_x[V_ind(0)], p1 = obj.m_x[V_ind(1)], p2 = obj.m_x[V_ind(2)];
        Plane_3 pl(p0, p1, p2);
        for (auto v: obj.m_mesh.vertices()){ obj.m_next_x[v] = (obj.m_x[v] += 2 * (pl.projection(obj.m_x[v]) - obj.m_x[v])); }
    }
    return res;
}
