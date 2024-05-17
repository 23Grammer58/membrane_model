//
// Created by alex on 13.04.2021.
//

#include "VanLoonBench.h"

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

bool ObSaveTag(const Object3D& obj, string filename, string tag){
    ofstream ob(filename, ios::binary | ios::trunc);
    if (!ob) return false;
    auto m_x = obj.m_mesh.property_map<V_ind, Point>(tag).first;
    for (auto& v: obj.m_mesh.vertices()){
        std::array<double, 3> P = {m_x[v][0], m_x[v][1], m_x[v][2]};
        ob.write(reinterpret_cast<const char *> ((P.data())), 3 * sizeof(double));
    }
    ob.close();
    return true;
}

bool ObReadTag(Object3D& obj, string filename, string tag){
    ifstream ob(filename, std::ios::binary);
    if (!ob) return false;
    auto m_x = obj.m_mesh.add_property_map<V_ind, Point>(tag);
    for (auto& v: obj.m_mesh.vertices()){
        std::array<double, 3> P;
        ob.read(reinterpret_cast<char *> ((P.data())), 3 * sizeof(double));
        m_x.first[v] = Point(P[0], P[1], P[2]);
    }
    ob.close();
    return true;
}

StretchFigure prepareInitialLeaf(double R, double phi, double Hc, double Hb, double Omega){
    auto flat_bnd = [l = R*phi, Hc, Hb](double t) -> Point_2 {
        if (t < 0) return Point_2(-l/2, Hc);
        if (t < 1) return Point_2(-l/2, Hc*(1-t));
        if (t < 3) return Point_2(l/2 * sin(M_PI_2*(t - 2)), -Hb * cos(M_PI_2*(t - 2)));
        if (t < 4) return Point_2(l/2, Hc*(t - 3));
        return Point_2( l/2, Hc);
    };
    auto cil_bnd = [flat_bnd, R](double t) -> Point{
        Point_2 x_f = flat_bnd(t);
        return Point(R*cos(x_f.x()/R), R*sin(x_f.x()/R), x_f.y());
    };
    auto proj_line = [cil_bnd, Hc](double Omega){
        return [&cil_bnd, Hc, Omega](double t) -> Point{
            Point x0 = cil_bnd(t);
            double dz = -(x0.z() - Hc);
            return Point(x0.x() + dz * tan(Omega), x0.y(), Hc);
        };
    };
    auto get_free_length = [proj_line](double Omega) {
        auto projection = proj_line(Omega);
        PiecewieceLineIntegral pli;
        pli.pieces.reserve(4);
        for (int i = 0; i < 4; ++i) {
            QuickLineIntegral qli;
            qli.m_curve = projection;
            qli.Tst = i, qli.Tend = i + 1;
            qli.rel = 1e-7;
            qli.Ninit = 20;
            pli.pieces.push_back(std::move(qli));
        }
        pli.setup();
        double L = pli.pieces.back().quick_index.back().l;
        return L;
    };
    std::cout << "L(" << Omega << ") = " << get_free_length(Omega) << std::endl;
    auto len_func = [cil_bnd, Hc, Omega](double t)->double{
        Point x0 = cil_bnd(t);
        double dz = -(x0.z() - Hc);
        auto t_scl = [](double t) {
            double lt = t - 2;
            if (t < 1 || t > 3) return 1.0;
            double t0 = 1.8, t1 = 2.2;
            if (t < t0) {
                double a = 1, b = t0;
                double pt = (t-(a+b)/2);
                double m = (b - a) / 2;
                return 0.8 + 0.2 * pt * pt / (m * m);
            }
            if (t > t1) {
                double a = t1, b = 3;
                double pt = (t-(a+b)/2);
                double m = (b - a) / 2;
                return 0.8 + 0.2 * pt * pt / (m * m);
            }
            double pt = t - 2;
            double m = t1 - 2;
            return 1.0 + 0.02*(1.0 - pt*pt / (m*m));
        };
        double res = abs(dz) / cos(Omega);
        res *= t_scl(t);
        return res;
    };

    auto res = StretchFigure()
            .setStretchDirection(Vector{sin(Omega), 0, cos(Omega)})
            .setStretchLengthFunc(len_func)
            .setCurveFunc(cil_bnd)
            .setPntParams({0, 1, 2, 3, 4});

    return res;
}

Object3D computeInitialLeaf(double R, double phi, double Hc, double Hb, double Omega, double mesh_h){
    auto flat_bnd = [l = R*phi, Hc, Hb](double t) -> Point_2 {
        if (t < 0) return Point_2(-l/2, Hc);
        if (t < 1) return Point_2(-l/2, Hc*(1-t));
        if (t < 3) return Point_2(l/2 * sin(M_PI_2*(t - 2)), -Hb * cos(M_PI_2*(t - 2)));
        if (t < 4) return Point_2(l/2, Hc*(t - 3));
        return Point_2( l/2, Hc);
    };
    auto cil_bnd = [&flat_bnd, R](double t) -> Point{
        Point_2 x_f = flat_bnd(t);
        return Point(R*cos(x_f.x()/R), R*sin(x_f.x()/R), x_f.y());
    };
    auto proj_line = [&cil_bnd, Hc](double Omega){
        return [&cil_bnd, Hc, Omega](double t) -> Point{
            Point x0 = cil_bnd(t);
            double dz = -(x0.z() - Hc);
            return Point(x0.x() + dz * tan(Omega), x0.y(), Hc);
        };
    };
    auto get_free_length = [&proj_line](double Omega) {
        auto projection = proj_line(Omega);
        PiecewieceLineIntegral pli;
        pli.pieces.reserve(4);
        for (int i = 0; i < 4; ++i) {
            QuickLineIntegral qli;
            qli.m_curve = projection;
            qli.Tst = i, qli.Tend = i + 1;
            qli.rel = 1e-7;
            qli.Ninit = 20;
            pli.pieces.push_back(std::move(qli));
        }
        pli.setup();
        double L = pli.pieces.back().quick_index.back().l;
        return L;
    };
    std::cout << "L(" << Omega << ") = " << get_free_length(Omega) << std::endl;
    auto len_func = [&cil_bnd, Hc, Omega](double t)->double{
        Point x0 = cil_bnd(t);
        double dz = -(x0.z() - Hc);
        auto t_scl = [](double t) {
            double lt = t - 2;
//            double st = 0.9;
//            if(t < 1 || t > 3) {
//                return (st + (1-st)*abs(lt)/2);
//            }
//            return st;//*(0.9 + 0.1*abs(lt));
            if (t < 1 || t > 3) return 1.0;
            double t0 = 1.8, t1 = 2.2;
            if (t < t0) {
                double a = 1, b = t0;
                double pt = (t-(a+b)/2);
                double m = (b - a) / 2;
                return 0.8 + 0.2 * pt * pt / (m * m);
            }
            if (t > t1) {
                double a = t1, b = 3;
                double pt = (t-(a+b)/2);
                double m = (b - a) / 2;
                return 0.8 + 0.2 * pt * pt / (m * m);
            }
            double pt = t - 2;
            double m = t1 - 2;
            return 1.0 + 0.02*(1.0 - pt*pt / (m*m));
        };
        auto t_scl1 = [](double t){
            double sc = 0.8;
            double add_sc = 0.15;//0.15;//0.05
            double rm_sc = 0.0;//0.1;//0.05
            double lt = t - 2;
            double reg = 1;
            double w = 0.5;
            auto make_bell_func = [](double mid, double weight, double width){
                return [m = mid, w = weight, r = width](double x){
                    return (exp(- (x-m)*(x-m) / (w * w)) -  exp(- (r * r) / (w * w))) / (1 - exp(- (r * r) / (w * w) ));
                };
            };
//            double f = (exp(- lt*lt / (w * w)) -  exp(- (reg * reg) / (w * w))) / (1 - exp(- (reg * reg) / (w * w) ));
//            return sc + add_sc * f;
            double mid = 0.5;//0.7;
            return sc + add_sc * make_bell_func(0, w, reg)(lt)
                - rm_sc * (make_bell_func(-mid, 0.4, 0.6)(lt) + make_bell_func(+mid, 0.4, 0.6)(lt));
        };
        double res = abs(dz) / cos(Omega);
        res *= t_scl1(t);
        return res;
    };

    auto leaf = StretchFigure()
            .setStretchDirection(Vector{sin(Omega), 0, cos(Omega)})
            .setStretchLengthFunc(len_func)
            .setCurveFunc(cil_bnd)
            .setPntParams({0, 1, 2, 3, 4})
            .generateMesh(Ani3dFSize(mesh_h).fsize);
    leaf.name = "deformable";

//    std::cout << "Wanted L = " << 2*M_PI/3 * R << std::endl; //-28.5 градусов или +6.5
//    for (int i = -30; i <= 20; ++i){
//        double Omega = i / 180.0 * M_PI;
//        std::cout << "L(" << i << ") = " << get_free_length(Omega) << std::endl;
//    }
    return leaf;
}

//auto t_scl = [](double t) {
//    if (t < 1 || t > 3) return 1.0;
//    double t0 = 1.8, t1 = 2.2;
//    if (t < t0) {
//        double a = 1, b = t0;
//        double pt = (t-(a+b)/2);
//        double m = (b - a) / 2;
//        return 0.8 + 0.2 * pt * pt / (m * m);
//    }
//    if (t > t1) {
//        double a = t1, b = 3;
//        double pt = (t-(a+b)/2);
//        double m = (b - a) / 2;
//        return 0.8 + 0.2 * pt * pt / (m * m);
//    }
//    double pt = t - 2;
//    double m = t1 - 2;
//    return 1.0 + 0.05*(1.0 - pt*pt / (m*m));
//};
//double R = 15, phi = 2*M_PI/3 - 0.4, Hc = 1.5, Hb = 15, Omega = -40.0/180 * M_PI, mesh_h = 0.5;
//хороший шаблон

//auto t_scl = [](double t) {
//    double lt = t - 2;
//    if (t < 1 || t > 3) return 1.0;
//    double t0 = 1.8, t1 = 2.2;
//    if (t < t0) {
//        double a = 1, b = t0;
//        double pt = (t-(a+b)/2);
//        double m = (b - a) / 2;
//        return 0.8 + 0.2 * pt * pt / (m * m);
//    }
//    if (t > t1) {
//        double a = t1, b = 3;
//        double pt = (t-(a+b)/2);
//        double m = (b - a) / 2;
//        return 0.8 + 0.2 * pt * pt / (m * m);
//    }
//    double pt = t - 2;
//    double m = t1 - 2;
//    return 1.0 + 0.02*(1.0 - pt*pt / (m*m));
//};
//double R = 15, phi = 2*M_PI/3 - 0.4, Hc = 1.5, Hb = 15, Omega = 25.0/180.0 * M_PI;
//неплохой шаблон

map<E_ind, BendingForce::BoundaryData::Entry>
    computeBndDataOnCilAligned(Object3D& leaf,
         std::function<bool(Object3D&, E_ind, V_ind*/*[2]*/, V_ind)> choose_edge,
         double angle_align)
 {
    map<E_ind, BendingForce::BoundaryData::Entry> bdata;
    Eigen::AngleAxis<double> q(-angle_align, Eigen::Vector3d(0, 0, 1));

    for (auto e: leaf.m_mesh.edges()) {
        auto f = face_around(leaf.m_mesh, e);
        if (f[0].second && f[1].second)
            continue;
        auto v = vert_around(leaf.m_mesh, e);
        auto vop = vert_opposite(leaf.m_mesh, e);
        if (choose_edge(leaf, e, v.data(), vop.first[0])) {
            Point _p[3] = {leaf.m_x[v[0]], leaf.m_x[v[1]], leaf.m_x[vop.first[0]]};
            Eigen::Vector3d p[3] = {{_p[0][0], _p[0][1], _p[0][2]}, {_p[1][0], _p[1][1], _p[1][2]}, {_p[2][0], _p[2][1], _p[2][2]}};
            for (int k = 0; k < 3; ++k) p[k] = q * p[k];
            double lphis[2] = {atan2(p[0].y(), p[0].x()), atan2(p[1].y(), p[1].x())};
            double lphi = (lphis[0] + lphis[1]) / 2;
            Eigen::Vector3d ephi(-sin(lphi), cos(lphi), 0);
            Eigen::Vector3d el(0, 0, 1);
            Eigen::Vector3d er = el.cross(ephi);
            Eigen::Vector3d lt = p[1] - p[0];
            lt.normalize();
            if (lt.dot(ephi) > 1e-12) lt *= -1;
            if (lphis[0] * lphis[1] > 1e-12 && lt.dot(el) * lphi > 0) lt *= -1;
            Eigen::Vector3d eb = er.cross(lt);
            eb.normalize();
            eb = q.inverse() * eb;
            eb.normalize();
            std::array<double, 3> ext_normal = {eb[0], eb[1], eb[2]};
            bdata[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, ext_normal);
        } else {
            bdata[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
        }
    }
    return bdata;
}

Mesh::Property_map<V_ind, Vector> addBdataTag(Object3D& leaf, const map<E_ind, BendingForce::BoundaryData::Entry>& bdata){
    auto bdit = leaf.m_mesh.add_property_map<V_ind, Vector>("v:bdata").first;
    for (auto v: leaf.m_mesh.vertices()) bdit[v] = CGAL::NULL_VECTOR;
    for (auto dat: bdata) {
        auto e = dat.first;
        auto v = vert_around(leaf.m_mesh, e);
        for (int i = 0; i < 2 && dat.second.btype == BendingForce::BoundaryData::CLAMPED; ++i)
            bdit[v[i]] += Vector(dat.second.normal[0], dat.second.normal[1], dat.second.normal[2]);
    }
    for (auto v: leaf.m_mesh.vertices()) {
        double norm = sqrt(bdit[v].squared_length());
        if (norm > 1e-7) bdit[v] /= norm;
    }
    return bdit;
}

int computeInitialMesh(ComputeInitMesh p){
    int status = 0;
    double R = p.R;
    double phi = p.phi, Hc = p.Hc, Hb = p.Hbc*R, Omega = p.Omega;//-30.0; +8.0;
    double mesh_h = p.mesh_h;
//    double H1 = 10, H2 = 4, R = 10, phi = 2*M_PI/3 * 0.9;
//    double mesh_h = 0.3;
    double P = p.P, mu = p.mu, Ht = p.Ht; //Ht = 0.8
    double delta = p.c_delta / (mu), err = p.err, abs_err = p.abs_err, diverge_lim = p.diverge_lim;
    int maxits = p.maxits, pfreq = p.pfreq;
    bool use_clamped = p.use_clamped;
    bool use_true_bc = p.use_true_bc;
    bool use_bending = p.use_bending;
    string postfix = "_" + std::to_string(static_cast<int>(2*R))+"_" + p.case_name;
    if (!use_bending) use_true_bc = use_clamped = false;
    if (use_true_bc && use_clamped) throw std::runtime_error("It doesn't allowed to use clamped-bc-forces and true neumann bc");
    Object3D leaf = computeInitialLeaf(R, phi, Hc, Hb, Omega, mesh_h);
    leaf.name = p.leaf_name;
    string dir = p.dir;
    // leaf.save(dir + "init_mesh_x0" + postfix + ".vtk", leaf.m_x0);
    // leaf.save(dir + "init_mesh_x"+ postfix + ".vtk");
    // ObSaveTag(leaf, "../result/ModelLeafBench/template/leaf_presew_x0"+ postfix + ".txt", "v:x");

    static const int FIXED = 2, FREE = 1|4|8;
    DirichletCond dc = [](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & FIXED) return 0;
//        if (obj.m_boundary[v] & FREE) return 3;
        return 7;
    };
    leaf.setDirichletCond(dc);
    auto bdata = computeBndDataOnCilAligned(leaf,
        [](Object3D& leaf, E_ind e, V_ind* v, V_ind vop){
        return (leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED) && !(leaf.m_boundary[vop] & FREE);
    });
    auto bdit = addBdataTag(leaf, bdata);
    auto evaluate_bnd_quality = [bdata](Object3D& obj){
        double q = 0;
        for (auto d: bdata){
            auto e = d.first;
            auto f = face_around(obj.m_mesh, e)[0].first;
            auto b = Eigen::Map<Eigen::Vector3d>(d.second.normal);
            Eigen::Vector3d n(obj.m_normal[f][0], obj.m_normal[f][1], obj.m_normal[f][2]);
            double q_loc = abs(b.dot(n));
            if (q_loc > q) q = q_loc;
        }
        return q;
    };

    Force elastic = NeoGookModel(mu, Ht, "../../../generated", false);
    elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
    Force bending = BendingForce(elastic.target<HyperElasticForceBase>()->f, false);
    bending.target<BendingForce>()->prepareJacobianFunction(false);
    Force pressure = SimplePressureLoad(P);
    Force bottom_constraint = SpurnPlane(p.bpln_Pc*P, Hb*p.bpln_Htc,  SpurnPlane::Plane_3(Point(0, 0, Hb*p.bpln_bilc), Vector(0, 0, 1)));
    double d_pln = p.fpln_posc*R, d_lr = p.spln_posc*R;//
    Force left_constraint = SpurnPlane(p.spln_Pc*P, p.spln_Htc*R, SpurnPlane::Plane_3(Point(d_lr, 0, 0), Vector(-sin(-p.spln_phi), cos(-p.spln_phi), 0.0)));
    Force right_constraint = SpurnPlane(p.spln_Pc*P, p.spln_Htc*R, SpurnPlane::Plane_3(Point(d_lr, 0, 0), -Vector(-sin(p.spln_phi), cos(p.spln_phi), 0.0)));
    Force face_constraint = SpurnPlane(p.fpln_Pc*P, p.fpln_Htc*R, SpurnPlane::Plane_3(Point(d_pln, 0, 0), Vector(1, 0, 0.2)));
    std::map<E_ind, Vector> cdat; for(const auto& i: bdata) if (i.second.btype == 2) cdat.insert({i.first, Vector{i.second.normal[0], i.second.normal[1], i.second.normal[2]}});
    Force clamped = ClampedBndPenalty(1.0, cdat);

    World w;
    auto oid = w.addObject3D(std::move(leaf), 1);
    if (use_true_bc) bending.target<BendingForce>()->set_boundary_data(w.obj(oid), bdata);
    w.addForce(std::move(pressure));
    w.addForce(std::move(elastic));
    w.addForce(bottom_constraint);
    w.addForce(left_constraint);
    w.addForce(right_constraint);
    w.addForce(face_constraint);
    Force_ID bfid;
    if (use_bending) bfid = w.addForce(std::move(bending), oid);
    Force_ID clfid;
    if (use_clamped) clfid = w.addForce(std::move(clamped), oid);

    // w.obj(oid).save("../result/VanLoonBench/temp/leaf2"+ postfix + ".vtk");//.m_mesh.add_property_map<F_ind, std::array<double, 3>>("")

    w.UpdateForces();
    double abs_err_init = w.getWorldNorm();
    thread* t;
    if (p.with_view){
        if (p.argc < 1 || p.argv == nullptr) throw std::runtime_error("View window can't be open without env information from argv[]");
        t = new thread(start_debug_gui, p.argc, p.argv);
    }
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());

    World3d::Timer sym_time;
    StaticForceApplierP sfa(delta);
    sfa.addDamper(Damper(0.01*Ht, 0.8, 100, 1000000));
    double temp_err = err;
    err = 0;
    std::vector<double> cmaxits = {};//8000, 4000};//, 4000, 4000, 4000, 5000};
    std::vector<double> stiffness = { 1,    3,   10,   30,  100,  400};
    auto stopCond = World3d::StaticStopCondition(err, diverge_lim, maxits, &pfreq, &sym_time, nullptr, true).setAbsErr(abs_err);
    for (int i = 0; i < cmaxits.size(); ++i) {
        if (use_clamped) w.obj(oid).m_forces[clfid].target<ClampedBndPenalty>()->setStiffnes(stiffness[i]);
        w.setForceApplier(sfa);
        maxits = cmaxits[i];
        status = w.Simulation(stopCond);
    }
    // w.obj(oid).save("../result/VanLoonBench/temp/init_mesh2"+ postfix + ".stl");

    err = temp_err;

    NSWorldWrapper nsww(w);
    NLProblem nlp = NSWorldWrapper::makeProblem(nsww);
    NonLinearSolverKinsol nlsp(nlp);
    LinearSolver ls("inner_mptiluc");
    ls.setVerbosityLevel(1);
    ls.SetParameterReal("relative_tolerance", 1e-20);
    double droptol = 8e-3;
    ls.SetParameterReal("drop_tolerance", droptol);
    ls.SetParameterReal("reuse_tolerance", 0.1 * droptol);
    nlsp.SetLinearSolver(&ls);
    static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
    nlsp.SetInfoHandlerFn([](const char *module, const char *function, char *msg){
        std::cout << "[" << module << "] " << function << "\n   " << msg << std::endl;
    });
    nlsp.SetVerbosityLevel(1);
    nlsp.SetInitialGuess(nsww.getCurrentX());
    nlsp.setInterIterationMonitorFunc([&nsww, &time = sym_time](const double * x){
        nsww.setX(x);
        std::cout << "time = " << time.elapsed() << std::endl;
        nsww.RenderIteration();
        return 0;
    });
    nlsp.SetFuncNormTol(1e-7 * abs_err_init);
    nlsp.SetNumMaxIters(50);
    nlsp.SetMaxNewtonStep(R / 10 * sqrt(nlp.m_dofs));
    bool slvFlag = false;
    nlsp.SetParameterInt("MaxSetupCalls", 1);
    nlsp.SetParameterInt("MaxSubSetupCalls", 1);
    nlsp.SetScaledStepTol(1e-9);
    nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
    nlsp.SetMaxBetaFails(40);

    slvFlag = nlsp.Solve();

    if (use_clamped) {
        nlsp.SetFuncNormTol(1e-11 * abs_err_init);
        nlsp.SetScaledStepTol(1e-11);
        double q = evaluate_bnd_quality(w.obj(oid)),
        stiff = w.obj(oid).m_forces[clfid].target<ClampedBndPenalty>()->m_K;
        double q_eps = 1e-2;
        while (q > q_eps){
            w.obj(oid).m_forces[clfid].target<ClampedBndPenalty>()->setStiffnes((stiff*=2));
            slvFlag = nlsp.Solve();
            std::cout << "bnd_quality = " << (q = evaluate_bnd_quality(w.obj(oid))) << std::endl;
        }
    }
//    w.obj(oid).remove_force(clfid);
//    w.obj(oid).m_forces[bfid].target<BendingForce>()->set_boundary_data(w.obj(oid), bdata);
//    slvFlag = nlsp.Solve();

    std::cout <<"\t#NonLinIts = "<<nlsp.GetNumNolinSolvIters() << " Residual = "<< nlsp.GetResidualNorm()
              << "\n\t#linIts = " << nlsp.GetNumLinIters() << " #funcEvals = " << nlsp.GetNumLinFuncEvals() << " #jacEvals = " << nlsp.GetNumJacEvals()
              << "\n\t#convFails = " << nlsp.GetNumLinConvFails() << " #betaCondFails = " << nlsp.GetNumBetaCondFails() << " #backtrackOps = " << nlsp.GetNumBacktrackOps() << "\n";
    std::cout << "Reason: " << nlsp.GetReason() << std::endl;

    // w.obj(oid).save("../result/ModelLeafBench/template/leaf_init"+ postfix + ".txt");
    // ObSaveTag(w.obj(oid), "../result/ModelLeafBench/template/leaf_init_x0"+ postfix + ".txt", "v:point");

    //compute quality of mesh deformation
    if (p.evalute_quality_of_res){
        auto& obj = w.obj(oid);
        std::cout << "bnd_quality = " << evaluate_bnd_quality(obj) << std::endl;
        auto Lit = set_L(obj);
        bool flat = check_flat_initial_template(obj);
        LocalCartesian _S(obj, flat);
        auto def = obj.m_mesh.add_property_map<F_ind, std::array<double, 2>>("f:deform").first;
        double q = 0, s_q[2] = {1, 1};
        for (auto f: obj.m_mesh.faces()){
            auto v = vert_around(obj.m_mesh, f);
            auto L = Eigen::Map<Eigen::Matrix<double, 2, 3>>(Lit[f].data());
            auto S = _S(f);
            Eigen::Matrix3d Q;
            Q <<    obj.m_x[v[0]][0], obj.m_x[v[1]][0], obj.m_x[v[2]][0],
                    obj.m_x[v[0]][1], obj.m_x[v[1]][1], obj.m_x[v[2]][1],
                    obj.m_x[v[0]][2], obj.m_x[v[1]][2], obj.m_x[v[2]][2];
            auto F = Q * L.transpose();
            Eigen::JacobiSVD<Eigen::Matrix<double, 3, 2>> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
            const auto& s = svd.singularValues();
            double q_loc = (s[0] / s[1] + s[1] / s[0]) / 2 - 1;
            if (q_loc > q) q = q_loc;
            s_q[0] = std::min({s_q[0], s[0], s[1]});
            s_q[1] = std::max({s_q[1], s[0], s[1]});
            def[f][0] = s[0], def[f][1] = s[1];
        }
        std::cout << "q = " << q << " {" << s_q[0] << ", " << s_q[1] << "}" << std::endl;
        auto bit = obj.m_forces[bfid].target<BendingForce>()->add_watch(
        [](BendingForce::VarContainer<BendingForce::LocalVar> &loc_vars,
           BendingForce::VarContainer<BendingForce::Invariant> &invariants,
           casadi::SX u,
           bool is_bnd) -> casadi::SX {
                using SX = casadi::SX;
                SX S(2, 2);
                for (const auto& i: invariants){
                    S += SX::simplify((SX::jacobian(u, i.sym) * i.C2dderiv));
                    S = SX::simplify(SX::substitute(S, i.sym, i.C2dval));
                    u = SX::substitute(u, i.sym, i.C2dval);
                }
                S *= 2;
                auto H_it = loc_vars.find("H");
                if (H_it == loc_vars.end())
                    throw std::runtime_error("To use this force you should provide parameter \"H\" that means height");
                SX H = H_it->sym;
                SX Ap = loc_vars["Ap"];
                SX Aq = loc_vars["Aq"];
                SX lambda = Ap / Aq;

                SX z = SX::sym("z");
                SX k = loc_vars["k_int"];
                SX k0 = loc_vars["k0"];
                SX nf = 0;
                SX B[3] = {SX(3, 6), SX(3, 6), SX(3, 6)};
                if (is_bnd) {
                    k = loc_vars["k_bnd"];
                    SX Bb = loc_vars["BbMatrix_bnd"];
                    SX C2d = loc_vars["C2d"];
                    SX SS = SX::simplify(S / (Ap * H));
                    SX D = SX::vertcat({SS(0, 0), SS(1, 1), SS(0, 1), SS(1, 0)});
                    D = SX::simplify(SX::jacobian(D, SX::horzcat({C2d(0, 0), C2d(1, 1), C2d(0, 1), C2d(1, 0)})));
                    auto Ds = SX::horzsplit(D, 1);
                    D = SX::horzcat({Ds[0], Ds[1], Ds[2] + Ds[3]});
                    Ds = SX::vertsplit(D, 1);
                    D = SX::vertcat({Ds[0], Ds[1], Ds[2] + Ds[3]});

                    SX xi = SX::vertcat({(k - k0)(0, 0), (k - k0)(1, 0), (k - k0)(2, 0)});
                    nf = SX::if_else_zero(loc_vars["b_lbl0"] == 1, 1) +
                         SX::if_else_zero(loc_vars["b_lbl1"] == 1, 1) +
                         SX::if_else_zero(loc_vars["b_lbl2"] == 1, 1);
                    SX n = SX::if_else_zero(loc_vars["b_lbl0"] == 1, loc_vars["b_n0"]) +
                           SX::if_else_zero(loc_vars["b_lbl1"] == 1, loc_vars["b_n1"]) +
                           SX::if_else_zero(loc_vars["b_lbl2"] == 1, loc_vars["b_n2"]);
                    SX t = SX::vertcat({-n(1, 0), n(0, 0)});
                    SX NN = SX::vertcat({SX::sq(n(0, 0)), SX::sq(n(1, 0)), n(0, 0) * n(1, 0)});
                    SX Tn = SX::vertcat(
                            {n(0, 0) * t(0, 0), t(1, 0) * n(1, 0), (n(1, 0) * t(0, 0) + n(0, 0) * t(1, 0)) / 2});
                    SX T = SX::vertcat({SX::sq(t(0, 0)), SX::sq(t(1, 0)), t(0, 0) * t(1, 0)});

                    SX A = SX::vertcat({SX::mtimes(NN.T(), D), SX::mtimes(Tn.T(), D), T.T()});

                    SX R = SX::vertcat({0, 0, T(0, 0) * xi(0, 0) + T(1, 0) * xi(1, 0) + 2 * T(2, 0) * xi(2, 0)});
                    R = SX::if_else_zero(nf < 1.5 && nf > 0.5, R);
                    SX invA = SX::if_else_zero(nf < 1.5 && nf > 0.5, SX::inv(A));
                    xi = SX::horzsplit(invA, 1)[2] * R(2, 0);
                    k = SX::if_else(nf < 1.5 && nf > 0.5, xi + k0, k);
                }

                SX res = SX::horzcat({k});

                for (int i = loc_vars.size() - 1; i >= 0; --i) {
                    res = SX::simplify(SX::substitute(res, loc_vars.at(i).sym, loc_vars.at(i).val));
                }
                return res;
            });
        auto k = obj.m_mesh.add_property_map<F_ind, std::array<double, 3>>("f:k_comp").first;
        auto kg = obj.m_mesh.add_property_map<F_ind, std::array<double, 2>>("f:k_gauss").first;
        obj.apply_updaters();
        double K[2] = {0, 0};
        for (auto f: obj.m_mesh.faces()) {
            auto dat = bit.eval_watches(obj, f);
            std::copy(dat[0].data(), dat[0].data()+3, k[f].data());
            Eigen::Matrix2d kmat;
            kmat << k[f][0], k[f][2],
                    k[f][2], k[f][1];
            auto vals = kmat.selfadjointView<Eigen::Upper>().eigenvalues();
            if (vals[0] < vals[1]) swap(vals[0], vals[1]);
            K[0] = std::min({K[0], vals[1]});
            K[1] = std::max({K[1], vals[0]});
            kg[f][0] = vals[0];
            kg[f][1] = vals[1];
        }
        std::cout << "K_gauss in [" << K[0] << ", " << K[1] << "]" << std::endl;
    }

    w.obj(oid).save(p.dir + "/leaf" + postfix + ".vtk");

    if (p.with_view){
        t->join();
        delete t;
    }

    return status;
}

#ifdef USE_MATIO
#include <matio.h>
#endif
int save_ModelLeaf_to_mat(string postfix = "_4_3"){
    Object3D leaf;
    //read sewed state
    leaf.read("../result/ModelLeafBench/template/leaf_init" + postfix + ".txt");
    //read flat state
    string init_aniso_crd = "v:initial";
    auto init_x = leaf.m_mesh.add_property_map<V_ind, Point>(init_aniso_crd).first;
    ObReadTag(leaf, "../result/ModelLeafBench/template/leaf_init_x0" + postfix + ".txt",
              init_aniso_crd);

    //read generated noflat state
    string presew_crd = "v:presew";
    auto presew_x = leaf.m_mesh.add_property_map<V_ind, Point>(presew_crd).first;
    ObReadTag(leaf, "../result/ModelLeafBench/template/leaf_presew_x0" + postfix + ".txt",
              presew_crd);

    Eigen::MatrixXd Xf(leaf.m_mesh.num_vertices(), 2), Xg(leaf.m_mesh.num_vertices(), 3), Xr(leaf.m_mesh.num_vertices(), 3);
    std::vector<unsigned> bnd(leaf.m_mesh.num_vertices());
    std::vector<unsigned> TRI(leaf.m_mesh.num_faces()*3);
    for (auto f: leaf.m_mesh.faces()) {
        auto v = vert_around(leaf.m_mesh, f);
        for (int k = 0; k < 3; ++k) TRI[f.idx() + k*leaf.m_mesh.num_faces()] = v[k].idx() + 1;
    }
    for (auto v: leaf.m_mesh.vertices()){
        for (int k = 0; k < 3; ++k) {
            Xr(v.idx(), k) = leaf.m_x[v][k];
            if (k != 2) Xf(v.idx(), k) = init_x[v][k];
            Xg(v.idx(), k) = presew_x[v][k];
        }
        bnd[v.idx()] = static_cast<unsigned>(leaf.m_boundary[v]);
    }
#ifdef USE_MATIO
    string path = "../result/ModelLeafBench/template/data" + postfix + ".mat";
    mat_t *matfp;
    matfp = Mat_CreateVer(path.c_str(), nullptr, MAT_FT_MAT5);
    if (!matfp ) throw std::runtime_error("Error while creating MAT file \"" + path + "\"");
    size_t dims[2];
    matvar_t *matvar;
    dims[0] = leaf.m_mesh.num_faces(), dims[1] = 3;
    matvar = Mat_VarCreate("TRI", MAT_C_UINT32, MAT_T_UINT32, 2, dims, TRI.data(), 0);
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);

    dims[0] = leaf.m_mesh.num_vertices(), dims[1] = 1;
    matvar = Mat_VarCreate("LBL", MAT_C_UINT32, MAT_T_UINT32, 2, dims, bnd.data(), 0);
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);

    dims[0] = Xr.rows(), dims[1] = Xr.cols();
    matvar = Mat_VarCreate("X", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, Xr.data(), 0);
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);

    dims[0] = Xf.rows(), dims[1] = Xf.cols();
    matvar = Mat_VarCreate("X_flat", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, Xf.data(), 0);
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);

    dims[0] = Xg.rows(), dims[1] = Xg.cols();
    matvar = Mat_VarCreate("X_generated", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, Xg.data(), 0);
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);

    Mat_Close(matfp);
    return 0;
#endif
    return -1;
}

int SaveDisplaceMap(){
    Object3D leaf;
    string init_config_prefix[3] = {"_4", "_4_2", "_4_3"};
    leaf.name = "ModelLeaf";
    leaf.read("../result/ModelLeafBench/template/leaf_init" + init_config_prefix[0] + ".txt");
    string init_aniso_crd = "v:initial";
    auto init_x = leaf.m_mesh.add_property_map<V_ind, Point>(init_aniso_crd).first;
    ObReadTag(leaf, "../result/ModelLeafBench/template/leaf_init_x0" + init_config_prefix[0] + ".txt",
              init_aniso_crd);//"v:point");//"v:initial");
    for (auto v: leaf.m_mesh.vertices()) leaf.m_next_x[v] = leaf.m_x[v], leaf.m_F[v] = CGAL::NULL_VECTOR;
    //set boundary conditions
    static const int FIXED = 2, FREE = 1|4|8;
    DirichletCond dc = [](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & FIXED) return 0;
        return 7;
    };
    leaf.setDirichletCond(dc);
    auto bdata = computeBndDataOnCilAligned(leaf,
                                            [](Object3D& leaf, E_ind e, V_ind* v, V_ind vop){
                                                return (leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED) && !(leaf.m_boundary[vop] & FREE);
                                            });
    {
        leaf.apply_updaters();
        for (auto &d: bdata) {
            auto e = d.first;
            auto f = face_around(leaf.m_mesh, e)[0].first;
            auto b = Eigen::Map<Eigen::Vector3d>(d.second.normal);
            Eigen::Vector3d n(leaf.m_normal[f][0], leaf.m_normal[f][1], leaf.m_normal[f][2]);
            b = b - b.dot(n) * n;
            b.normalize();
        }
    }
    auto bdit = addBdataTag(leaf, bdata);
    ObReadTag(leaf, "../result/ModelLeafBench/displ_map/21.txt","v:21");
    ObReadTag(leaf, "../result/ModelLeafBench/displ_map/45.txt","v:45");
    ObReadTag(leaf, "../result/ModelLeafBench/displ_map/69.txt","v:69");
    ObReadTag(leaf, "../result/ModelLeafBench/displ_map/93.txt","v:93");
    leaf.save("../result/ModelLeafBench/displ_map/together.vtk");
    return 0;
}

int ModelLeafBench1(int argc, char* argv[]){
    double Ht = 0.4;
    double mu = 9e2; //kPa
    double theta = M_PI_4, k1 = 1.1e3 /*kPa*/, k2 = 1e3, kappa = 0.15;//.29;
    double P = 90.0_mmHg / 1e3; //kPa
    int scenario = 0; // 0 - membrane, 1 - shape, 2 - shape + BC
    int templ = 0;
    double mesh_h[3] = {0.5, 0.25, 0.125};
    int nit = 0, maxnit = 1000;
    double f_err_lim = 1e-4;
    bool start_from_previous = true;
    bool do_solve = true;
    bool process_result = true;
    bool regenerate = false;
    bool view = false;
    double delta_relax = 1e-4;

    struct CustomSolverIteration{
        int relax_its = 1000, relax_pfreq = 50;
        int ns_type = 0;
        double delta_full_step = 0.1;
        int max_ns_its = 100;
    };
    std::vector<CustomSolverIteration> solverItData;

    string  main_save = "../result/ModelLeafBench/",
            to_save = "",
            extr_dir = "extra/",
            to_gen = "../generated/";

    int _case = -1;

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << "\n";
            std::cout << " Command line options: " << "\n";
            std::cout << "  -sc, --scenario     <Choosing of the running scenario: 0 - membrane, 1 - shape, 2 - shape + BC>" << "\n";
            std::cout << "  -ht, --thickness    <Leafs thickness>" << "\n";
            std::cout << "  -pr, --pressure     <Diastolic pressure>" << "\n";
            std::cout << "  -mu, --elastic_mod  <Shear modulus>" << "\n";
            std::cout << "  -an, --aniso        <Parameters of anisotropy: k1, k2, kappa, theta>" << "\n";
            std::cout << "  -tl, --template     <Which template: 0 - 0.5 mm, 1 - 0.25 mm, 2 - 0.125 mm>" << "\n";
            std::cout << "  -tr, --target       <Directory in common directory to save results>" << "\n";
            std::cout << "  -md, --maindir      <Common directory to save results from different simulations>" << "\n";
            std::cout << "  -gd, --gendir       <Path to save generated codes>" << "\n";
            std::cout << "  -sp, --slvstop      <Solver stop criterion parameters: f_abs_err, maxnits>" << "\n";
            std::cout << "  -cs, --case         <Set default parameters for the case>" << "\n";
            std::cout << "  -nc, --no_continue  <Do not start from previous calculations>" << "\n";
            std::cout << "  -ns, --no_solve     <Do not solve the problem>" << "\n";
            std::cout << "  -np, --no_process   <Do not process result of solution>" << "\n";
            std::cout << "  -rg, --regenerate   <Regenerate symbolic functions>" << "\n";
            std::cout << "  -dr, --delta_relax  <Set initial relaxation delta>" << "\n";
            std::cout << "  -vw, --view         <Create ImGui debug view window>" << "\n";
            std::cout << "  -ss, --spec_solver  <Add specific solver parameters for last iterations followed "
                         "after previous added specific parameters>\n"
                         "                       Syntax is follow: NONLINEAR_SLV_TYPE NONLINEAR_SLV_TYPE_PARAMETERS [RELAX_PARAMS]\n"
                         "                       NONLINEAR_SLV_TYPE: KINSOL or NEWTON\n"
                         "                       NONLINEAR_SLV_TYPE_PARAMETERS: for KINSOL <maxnits>,\n"
                         "                                                      for NEWTON <maxnits full_step>\n"
                         "                       RELAX_PARAMS: <relax_its relax_pfreq>" << "\n";
            std::cout << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-vw") == 0 || strcmp(argv[i], "--view") == 0){
            view = true;
            continue;
        }
        if (strcmp(argv[i], "-sc") == 0 || strcmp(argv[i], "--scenario") == 0) {
            try {
                if (i + 1 < argc) {
                    scenario = stoi(argv[++i]);
                    if (scenario < 0 || scenario > 3){
                        throw std::runtime_error("Wrong value of scenario, you should choose one of the following:\n\t0 - membrane\n\t1 - shape\n\t2 - shape + BC\n");
                    }
                }
            } catch (std::exception& e){
                std::cout << "Waited number of scenario but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-ht") == 0 || strcmp(argv[i], "--thickness") == 0) {
            try {
                if (i + 1 < argc)
                    Ht = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited thickness value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-pr") == 0 || strcmp(argv[i], "--pressure") == 0) {
            try {
                if (i + 1 < argc)
                    P = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited thickness value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-mu") == 0 || strcmp(argv[i], "--elastic_mod") == 0) {
            try {
                if (i + 1 < argc)
                    mu = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited mu value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-an") == 0 || strcmp(argv[i], "--aniso") == 0) {
            int j = i;
            const char* v[] = {"k1", "k2", "kappa", "theta"};
            try {
                if (i + 1 < argc) k1 = stof(argv[++i]);
                if (i + 1 < argc) k2 = stof(argv[++i]);
                if (i + 1 < argc) kappa = stof(argv[++i]);
                if (i + 1 < argc) theta = stof(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-tl") == 0 || strcmp(argv[i], "--template") == 0) {
            try {
                if (i + 1 < argc) {
                    templ = stoi(argv[++i]);
                    if (templ < 0 || templ > 3){
                        throw std::runtime_error("Wrong value of template, you should choose one of the following:\n\t0 - 0.5mm\n\t1 - 0.25mm\n\t2 - 0.125mm\n");
                    }
                }
            } catch (std::exception& e){
                std::cout << "Waited number of template but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-tr") == 0 || strcmp(argv[i], "--target") == 0) {
            if (i + 1 < argc) to_save = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-md") == 0 || strcmp(argv[i], "--maindir") == 0) {
            if (i + 1 < argc) main_save = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-gd") == 0 || strcmp(argv[i], "--gendir") == 0) {
            if (i + 1 < argc) to_gen = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-sp") == 0 || strcmp(argv[i], "--slvstop") == 0) {
            int j = i;
            const char* v[] = {"err", "maxints"};
            try {
                if (i + 1 < argc) f_err_lim = stof(argv[++i]);
                if (i + 1 < argc) maxnit = stoi(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-cs") == 0 || strcmp(argv[i], "--case") == 0) {
            try {
                if (i + 1 < argc) {
                    int cs = stoi(argv[++i]);
                    if (cs < 0 || cs >= 96+32){
                        throw std::runtime_error("Wrong value of case, case should be integer in [0, 95]");
                    }
                    _case = cs;
                    int tetacs, kcs, mucs;
                    if (cs < 96) {
                        int cluster = (cs / 48); cs %= 48;
                        tetacs = cs % 3; cs /= 3;
                        kcs = cs % 4; cs /= 4;
                        mucs = 2 * cluster + cs % 2; cs /= 2;
                        int sc_cs = cs % 2;
                        scenario = (sc_cs == 1) ? 2 : 0;
                    } else if (cs < 96 + 32){
                        cs -= 96;
                        tetacs = 3 + cs % 2; cs /= 2;
                        kcs = cs % 4; cs /= 4;
                        mucs = 1 + 2 * (cs % 2); cs /= 2;
                        scenario = ((cs % 2) == 1) ? 2 : 0;
                    }
                    switch (tetacs) {
                        case 0: theta = 0; break;
                        case 1: theta = M_PI_4; break;
                        case 2: theta = M_PI_2; break;
                        case 3: theta = M_PI_4 / 2; break;
                        case 4: theta = 3 * M_PI_4 / 2; break;
                    }
                    switch (kcs) {
                        case 0: kappa = 0; break;
                        case 1: kappa = 0.2; break;
                        case 2: kappa = 0.29; break;
                        case 3: kappa = 1.0/3; break;
                    }
                    switch (mucs) {
                        case 0: mu = 9e1; break;
                        case 1: mu = 9e2; break;
                        case 2: mu = 3e2; break;
                        case 3: mu = 3e3; break;
                    }
                }
            } catch (std::exception& e){
                std::cout << "Waited number of scenario but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-nc") == 0 || strcmp(argv[i], "--no_continue") == 0) {
            start_from_previous = false;
            continue;
        }
        if (strcmp(argv[i], "-ns") == 0 || strcmp(argv[i], "--no_solve") == 0) {
            do_solve = false;
            continue;
        }
        if (strcmp(argv[i], "-np") == 0 || strcmp(argv[i], "--no_process") == 0) {
            process_result = false;
            continue;
        }
        if (strcmp(argv[i], "-rg") == 0 || strcmp(argv[i], "--regenerate") == 0) {
            regenerate = true;
            continue;
        }
        if (strcmp(argv[i], "-dr") == 0 || strcmp(argv[i], "--delta_relax") == 0) {
            if (i + 1 < argc) {
                delta_relax = stof(argv[++i]);
            }
            continue;
        }
        if (strcmp(argv[i], "-ss") == 0 || strcmp(argv[i], "--spec_solver") == 0) {
            CustomSolverIteration csi;
            if (i + 1 < argc) {
                if (strcmp(argv[i+1], "KINSOL") == 0) csi.ns_type = 0;
                else if (strcmp(argv[i+1], "NEWTON") == 0) csi.ns_type = 1;
                else throw std::runtime_error("Wrong NONLINEAR_SLV_TYPE \"" + string(argv[i+1]) + "\". Waited \"KINSOL\" or \"NEWTON\"" );
                ++i;
            }
            if (i + 1 < argc){
                if (csi.ns_type == 0){
                    csi.max_ns_its = stoi(argv[++i]);
                } else if (csi.ns_type == 1){
                    csi.max_ns_its = stoi(argv[++i]);
                    if (i + 1 < argc){
                        csi.delta_full_step = stof(argv[++i]);
                    }
                }
            }

            if (i + 1 < argc && argv[i+1][0] != '-'){
                csi.relax_its = stoi(argv[++i]);
            }
            if (i + 1 < argc && argv[i+1][0] != '-'){
                csi.relax_pfreq = stoi(argv[++i]);
            }
            solverItData.push_back(csi);
            continue;
        }
        if (true){
            throw std::runtime_error("Faced unparsed value: \"" + std::string(argv[i]) + "\"");
        }
    }

    if (to_save.empty()){
        if (_case >= 0) to_save = "case" + to_string(_case) + "_" + std::to_string(templ) + "/";
        else {
            std::ostringstream oss;
            auto tt = std::time(nullptr);
            auto tm = *std::localtime(&tt);
            oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
            to_save = oss.str() + "/";
        }
    }
    if (to_save[to_save.size()-1] != '/') to_save += "/";
    if (solverItData.empty()){
        if (scenario == 0){
            CustomSolverIteration csi;
            csi.ns_type = 0;
            csi.max_ns_its = 200;
            csi.relax_its = 10000;
            csi.relax_pfreq = 200;
            solverItData.push_back(csi);
        } else {
            CustomSolverIteration csi;
            csi.ns_type = 0;
            csi.max_ns_its = 200;
            csi.relax_its = 5000;
            csi.relax_pfreq = 100;
            solverItData.push_back(csi);
            csi.relax_its = 1000;
            csi.relax_pfreq = 50;
            solverItData.push_back(csi);
        }
    }

    auto print_input = [&](std::ostream& out) -> std::ostream& {
        out << "Input parameters: " <<
                   "\n\tscenario = " << scenario <<
                   "\n\tHt = " << Ht <<
                   "\n\tP = "  << P <<
                   "\n\tmu = " << mu <<
                   "\n\tk1 = " << k1 << ", k2 = " << k2 << ", kappa = " << kappa << ", theta = " << theta <<
                   "\n\tmesh_h = " << mesh_h[templ] <<
                   "\n\tf_abs_err = " << f_err_lim << ", maxnits = " << maxnit <<
                   "\n\tmain_dir = \"" << main_save << "\", save_dir = \"" << to_save <<
                   "\", gen_dir = \"" << to_gen << "\"" << ", extr_dir = \"" << extr_dir << "\"" <<
                   "\n\tinitial_relaxation_delta = " << delta_relax <<
                   "\n\tstart_from_previous = " << start_from_previous << ", do_solve = " << do_solve <<
                        ", process_result = " << process_result << ", regenerate = " << regenerate << "\n";
        out << "\tCustomSolverIteration{\n";
        string NLS_TYPES[2] = {"KINSOL", "NEWTON"};
        for (int i = 0; i < solverItData.size(); ++i){
            out << "\t  " << i << ":{" << NLS_TYPES[solverItData[i].ns_type] << " [" << solverItData[i].max_ns_its << " ";
            if (solverItData[i].ns_type == 1) out << solverItData[i].delta_full_step << "] ";
            else out << "] ";
            out << "RELAX [" << solverItData[i].relax_its << " " << solverItData[i].relax_pfreq << "] }\n";
        }
        out << "\t}" << std::endl;
        return out;
    };

    print_input(std::cout);
    {
        auto make_dir = [](std::string save) {
            path path = save;
            auto s = status(path);
            if (!status_known(s) || s.type() != file_type::directory)
                create_directory(path);
        };
        make_dir(main_save);
        make_dir(main_save + to_save);
        make_dir(main_save + to_save + extr_dir);
        make_dir(to_gen);
    }
    to_save = main_save + to_save;
    std::string extr_save = to_save + extr_dir;
    std::ofstream lgfile;
    if (start_from_previous) lgfile.open(to_save + "log.txt", std::ios_base::app);
    else lgfile.open(to_save + "log.txt");
    print_input(lgfile);

    Object3D leaf;
    string init_config_prefix[3] = {"_4", "_4_2", "_4_3"};
    leaf.name = "ModelLeaf";
    leaf.read("../result/ModelLeafBench/template/leaf_init" + init_config_prefix[templ] + ".txt");
    string init_aniso_crd = "v:initial";
    auto init_x = leaf.m_mesh.add_property_map<V_ind, Point>(init_aniso_crd).first;
    ObReadTag(leaf, "../result/ModelLeafBench/template/leaf_init_x0" + init_config_prefix[templ] + ".txt",
              init_aniso_crd);//"v:point");//"v:initial");
    for (auto v: leaf.m_mesh.vertices()) leaf.m_next_x[v] = leaf.m_x[v], leaf.m_F[v] = CGAL::NULL_VECTOR;
    //set anisotropy directions
    auto aniso_f = leaf.m_mesh.add_property_map<F_ind, std::array<DReal , 3>>("f:aniso_f").first;
    auto aniso_s = leaf.m_mesh.add_property_map<F_ind, std::array<DReal , 3>>("f:aniso_s").first;
    {
        auto _init_x = leaf.m_mesh.property_map<V_ind, Point>(init_aniso_crd);
        if (!_init_x.second) throw std::runtime_error("Tag \"" + init_aniso_crd + "\" not found");
        auto& init_x = _init_x.first;
        for (auto f: leaf.m_mesh.faces()) {
            auto& o = leaf;
            auto v = vert_around(leaf.m_mesh, f);
            Eigen::Matrix3d Q, P;
            Q <<    o.m_x0[v[0]][0], o.m_x0[v[1]][0], o.m_x0[v[2]][0],
                    o.m_x0[v[0]][1], o.m_x0[v[1]][1], o.m_x0[v[2]][1],
                    o.m_x0[v[0]][2], o.m_x0[v[1]][2], o.m_x0[v[2]][2];
            P <<    init_x[v[0]][0], init_x[v[1]][0], init_x[v[2]][0],
                    init_x[v[0]][1], init_x[v[1]][1], init_x[v[2]][1],
                    init_x[v[0]][2], init_x[v[1]][2], init_x[v[2]][2];
            Eigen::Vector3d n0 = (P.col(1) - P.col(0)).cross(P.col(2) - P.col(0));
            double Ap = n0.stableNorm() / 2;
            n0 /= (2*Ap);
            Eigen::Matrix3d Dp;
            for (int i = 0; i < 3; ++i) {
                Dp.col(i) = n0.cross(P.col((i + 2) % 3) - P.col((i + 1) % 3));
            }
            Dp /= (2*Ap);
            Eigen::Vector3d Mfp(cos(theta),  sin(theta), 0), Msp(cos(theta),  -sin(theta), 0);
            auto F = Q * Dp.transpose();
            Eigen::Vector3d Mf = F * Mfp, Ms = F * Msp;
            Mf.normalize(), Ms.normalize();
            aniso_f[f] = std::array<DReal , 3>{Mf[0],  Mf[1], Mf[2]};
            aniso_s[f] = std::array<DReal , 3>{Ms[0],  Ms[1], Ms[2]};
        }
    }
    //tags for debug view
    {
        auto vaniso_f = leaf.m_mesh.add_property_map<V_ind, Vector>("v:aniso_f").first;
        auto vaniso_s = leaf.m_mesh.add_property_map<V_ind, Vector>("v:aniso_s").first;
        for (auto v: leaf.m_mesh.vertices()) {
            vaniso_f[v] = vaniso_s[v] = CGAL::NULL_VECTOR;
            auto fs = face_around(leaf.m_mesh, v);
            for (auto f: fs) {
                Vector ff(aniso_f[f][0], aniso_f[f][1], aniso_f[f][2]), ss(aniso_s[f][0], aniso_s[f][1],
                                                                           aniso_s[f][2]);
                vaniso_f[v] += ff / fs.size();
                vaniso_s[v] += ss / fs.size();
            }
        }
    }
    //set boundary conditions
    static const int FIXED = 2, FREE = 1|4|8;
    DirichletCond dc = [](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & FIXED) return 0;
        return 7;
    };
    leaf.setDirichletCond(dc);
    auto bdata = computeBndDataOnCilAligned(leaf,
                                            [](Object3D& leaf, E_ind e, V_ind* v, V_ind vop){
                                                return (leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED) && !(leaf.m_boundary[vop] & FREE);
                                            });
    {
        leaf.apply_updaters();
        for (auto &d: bdata) {
            auto e = d.first;
            auto f = face_around(leaf.m_mesh, e)[0].first;
            auto b = Eigen::Map<Eigen::Vector3d>(d.second.normal);
            Eigen::Vector3d n(leaf.m_normal[f][0], leaf.m_normal[f][1], leaf.m_normal[f][2]);
            b = b - b.dot(n) * n;
            b.normalize();
        }
    }
    auto bdit = addBdataTag(leaf, bdata);
    //read previous computed step if neccessary
    int ait = -1;
    path path = to_save;
    auto _s = status(path);
    if (_s.type() == file_type::directory){
        if (start_from_previous) {
            for (auto &p: directory_iterator(to_save)) {
                if (p.status().type() == file_type::regular) {
                    string name = p.path().filename();
                    auto pos = name.find_last_of('=');
                    if (pos != name.npos && name.size() - pos > 1) {
                        string nitstr = name.substr(name.find_last_of('=') + 1);
                        nitstr = nitstr.substr(0, nitstr.find_first_of('_'));
                        int nit = stoi(nitstr);
                        if (nit > ait) ait = nit;
                    }
                }
            }
        }
        if (ait >= 0){
            ObReadTag(leaf, to_save + "ModelLeaf_it="+ to_string(ait) +"_tag_x.txt","v:x");
        }
        ait++;
    }

    /*initialize forces*/
    SimplePressureLoad Pf(P);
    HGO_GSTModel Ef(Ht, mu, k1, k2, kappa, theta, to_gen, regenerate, 15);
    Ef.prepareJacobianFunction(regenerate); Ef.f.setVerbosity(2);
    BendingForce bf(Ef.f, regenerate);
    bf.prepareJacobianFunction(regenerate); bf.setVerbosity(2);
    double mpR = 15;
    SpurnPlane left_constraint(1.1*P, 0.01*mpR, SpurnPlane::Plane_3(CGAL::ORIGIN, Vector(-sin(-M_PI/3), cos(-M_PI/3), 0.0)));
    SpurnPlane right_constraint(1.1*P, 0.01*mpR, SpurnPlane::Plane_3(CGAL::ORIGIN, -Vector(-sin(M_PI/3), cos(M_PI/3), 0.0)));
    std::map<E_ind, Vector> cdat; for(const auto& i: bdata)
        if (i.second.btype == BendingForce::BoundaryData::CLAMPED)
            cdat.insert({i.first, Vector{i.second.normal[0], i.second.normal[1], i.second.normal[2]}});
    ClampedBndPenalty clamped(0, cdat);
    /*set world*/
    World w;
    auto oid = w.addObject3D(std::move(leaf), 1);
    /*set symmetrical leafs*/
    Object3D leaf1 = w.obj(oid), leaf2 = w.obj(oid);
    leaf1.name = "ModelLeaf1"; leaf2.name = "ModelLeaf2";
    auto oid1 = w.addObject3D(std::move(leaf1), 1);
    auto oid2 = w.addObject3D(std::move(leaf2), 1);
    std::array<ObjectID, 3> oids {oid, oid1, oid2};
    auto copy_update = [&w, oids](){
        auto& mob = w.obj(oids[0]);
        for (int i = 1; i < 3; ++i) {
            Object3D& obj = w.obj(oids[i]);
            Eigen::AngleAxis<double> q(2.0 * M_PI / 3.0 * i, Eigen::Vector3d(0, 0, 1));
            for (auto v: obj.m_mesh.vertices()){
                Eigen::Vector3d p(mob.m_x[v][0], mob.m_x[v][1], mob.m_x[v][2]);
                Eigen::Vector3d pc = q * p;
                obj.m_next_x[v] = obj.m_x[v] = Point(pc[0], pc[1], pc[2]);
            }
        }
    };
    copy_update();
    w.objs()[oid1].second = 1;
    w.objs()[oid2].second = 1;

    bool slvFlag = false;
    std::thread* view_window;
    if (do_solve){
        /*set forces*/
        Force_ID fpi, fei, fbti, flti, frti, fbi, fcli;
        fpi = w.addForce(std::move(Pf), oid);
        fei = w.addForce(std::move(Ef), oid);
        flti = w.addForce(left_constraint, oid);
        frti = w.addForce(right_constraint, oid);
        if (scenario == 1 || scenario == 2) {
            fbi = w.addForce(std::move(bf), oid);
        }
        if (scenario == 2) {
                w.obj(oid).m_forces[fbi].target<BendingForce>()->set_boundary_data(w.obj(oid), bdata);
        }
        //initial saving
        if (ait == 0)
            w.obj(oids[0]).save(to_save + w.obj(oids[0]).name + "_init_state" + ".vtk", w.obj(oids[0]).m_x0);

        /*init solvers*/
        NSWorldWrapper nsww(w);
        NLProblem nlp = NSWorldWrapper::makeProblem(nsww);
        LinearSolver ls("inner_mptiluc");
        ls.setVerbosityLevel(1);
        ls.SetParameterReal("relative_tolerance", 1e-20);
        double droptol = 1e-4;
        if (templ != 0) droptol = 8e-3;
        ls.SetParameterReal("drop_tolerance", droptol);
        ls.SetParameterReal("reuse_tolerance", 0.1 * droptol);
        auto monitor_fcn = [&nsww, &copy_update](const double *x) {
            nsww.setX(x);
            copy_update();
            nsww.RenderIteration();
            return 0;
        };

        //initialize KINSOL solver
        NonLinearSolverKinsol nlsp(nlp);
        nlsp.SetLinearSolver(&ls);
        static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
        nlsp.SetVerbosityLevel(1);
        nlsp.SetInitialGuess(nsww.getCurrentX());
        nlsp.SetFuncNormTol(f_err_lim);
        nlsp.SetMaxNewtonStep(mpR / 10 * sqrt(nlp.m_dofs));
        nlsp.SetParameterInt("MaxSetupCalls", 1);
        nlsp.SetParameterInt("MaxSubSetupCalls", 1);
        nlsp.SetScaledStepTol(1e-7);
        nlsp.SetMaxBetaFails(40);
        nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
        World3d::Timer sym_time, com_time;
        com_time.reset();
        nlsp.SetInfoHandlerFn([&time = sym_time, &com_time, &lgfile](const char *module, const char *function, char *msg) {
            std::cout << "[" << module << "] " << function << "\n   " << msg << "\n" << "time = " << time.elapsed()
                      << " com_time = " << com_time.elapsed() << std::endl;
            lgfile << "[" << module << "] " << function << "\n   " << msg << "\n" << "time = " << time.elapsed()
                   << " com_time = " << com_time.elapsed() << std::endl;
        });
        nlsp.setInterIterationMonitorFunc(monitor_fcn);
        auto result_str_nlsp = [](NonLinearSolverKinsol &nlsp) {
            std::stringstream sstr;
            sstr << "\t#NonLinIts = " << nlsp.GetNumNolinSolvIters() << " Residual = " << nlsp.GetResidualNorm()
                 << "\n\t#linIts = " << nlsp.GetNumLinIters() << " #funcEvals = " << nlsp.GetNumFuncEvals()
                 << " #jacEvals = " << nlsp.GetNumJacEvals()
                 << "\n\t#convFails = " << nlsp.GetNumLinConvFails() << " #betaCondFails = "
                 << nlsp.GetNumBetaCondFails() << " #backtrackOps = " << nlsp.GetNumBacktrackOps() << "\n";
            sstr << "Reason: " << nlsp.GetReason();
            return sstr.str();
        };

        //initialize custom full_step NEWTON solver
        NonLinearSolverCustom nlsc;
        nlsc.setProblem(nlp);
        nlsc.setInterIterationMonitorFunc(monitor_fcn);
        nlsc.SetLinearSolver(ls);
        nlsc.SetInitialGuess(nsww.getCurrentX());
        nlsc.SetNSWorldWrapper(nsww);
        double nlsc_step = solverItData[0].delta_full_step;
        int nlsc_maxnits = solverItData[0].max_ns_its;
        nlsc.tau_algo = [&nlsc_step, &lgfile](int it, double rel_resid, double abs_resid, NonLinearSolverCustom& ns)->double{
            std::cout << "\tit = " << it << ": norm = " << abs_resid << " rel = " << rel_resid << std::endl;
            lgfile << "\tit = " << it << ": norm = " << abs_resid << " rel = " << rel_resid << std::endl;
            return nlsc_step;
        };
        nlsc.stop_cond = [&nlsc_maxnits, f_err_lim](int it, double rel, double abs, int& reason, NonLinearSolverCustom& ns){
            if (abs < f_err_lim) {
                return reason = 0, true;
            } else if (abs > 1e15 || std::isnan(abs)) {
                return reason = -1, true;
            }
            if (it >= nlsc_maxnits) {
                return reason = 2, true;
            }
            return reason = 1, false;
        };
        auto result_str_nlsc = [](NonLinearSolverCustom &nlsc) {
            std::stringstream sstr;
            sstr << "\t#NonLinIts = " << nlsc.GetNumNolinSolvIters() << " Residual = " << nlsc.GetResidualNorm()
                 << "\n\t#linIts = " << nlsc.GetNumLinIters() << " #funcEvals = " << nlsc.GetNumFuncEvals()
                 << " #jacEvals = " << nlsc.GetNumJacEvals() <<  "\n";
            sstr << "Reason: " << nlsc.GetReason();
            return sstr.str();
        };

        //initialize relaxation solver
        auto relaxate_solver = [Ht, &w, &sym_time, &lgfile, oid](double &delta, int nits,
                                                                                     int pfreq) {
            if (nits <= 0) return;
            std::cout << "delta = " << delta << "\n";
            lgfile << "delta = " << delta << "\n";
            static double diverge_lim = 1e15;
            static double err = 1e-4;
            auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &nits, &pfreq, &sym_time,
                                                         &lgfile,
                                                         true);
            StaticForceApplierP sfa(delta);
            sfa.addDamper(Damper(0.01 * Ht, 0.8, 15, 1000000).setInitDxScale(3.0)
                                  .setReInitDxParams(-1, -1).setMinDeltaFactor(1e-5));
            w.setForceApplier(sfa);
            sym_time.reset();
            w.Simulation(stopCond);
            double factor = w.obj(
                    oid)._next_x_updater.target<StaticForceApplierP>()->_dampers[0].target<Damper>()->getUsedDeltaFactor();
            delta *= factor;
            if (1 - factor > 1e-2) std::cout << "New relaxation delta = " << delta << std::endl;
        };

        int func_evals = 0, jac_evals = 0;
        auto print_iterate_info = [&fe = func_evals, &je = jac_evals, ait, &com_time, &lgfile](int it){
            std::cout << "ITERATE INFO: it: " << it << " com_time = " << com_time.elapsed() <<
                      " accum_func_evals = " << fe << " accum_jac_evals = " << je << std::endl;
            lgfile << "ITERATE INFO: it: " << it << " com_time = " << com_time.elapsed() <<
                   " accum_func_evals = " << fe << " accum_jac_evals = " << je << std::endl;
        };

        if (view){
            view_window = new std::thread(start_debug_gui, argc, argv);
            w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
        }
        int it = ait, s_id = 0;
        for (; it < maxnit && !slvFlag; ++it) {
            s_id = it - ait;
            if (it - ait >= solverItData.size()) s_id = solverItData.size()-1;
            auto csi = solverItData[s_id];
            relaxate_solver(delta_relax, csi.relax_its, csi.relax_pfreq);
            func_evals += csi.relax_its;

            sym_time.reset();
            if (csi.ns_type == 0){
//                ls.SetParameterReal("drop_tolerance", droptol);
//                ls.SetParameterReal("reuse_tolerance", 0.1 * droptol);
                nlsp.SetInitialGuess(nsww.getCurrentX());
                nlsp.SetNumMaxIters(csi.max_ns_its);
                slvFlag = nlsp.Solve();
                lgfile << result_str_nlsp(nlsp) << std::endl;
                std::cout << result_str_nlsp(nlsp) << std::endl;
                func_evals += nlsp.GetNumFuncEvals(), jac_evals += nlsp.GetNumJacEvals();
            } else if (csi.ns_type == 1){
//                ls.SetParameterReal("drop_tolerance", 1e-3);
//                ls.SetParameterReal("reuse_tolerance", 0.1 * 1e-3);
                ls.SetParameterReal("absolute_tolerance", 1e-9);
                nlsc_step = csi.delta_full_step;
                nlsc_maxnits = csi.max_ns_its;
                nlsc.SetInitialGuess(nsww.getCurrentX());
                slvFlag = nlsc.Solve();
                slvFlag  = (nlsc.getReasonFlag() == 0);
                lgfile << result_str_nlsc(nlsc) << std::endl;
                std::cout << result_str_nlsc(nlsc) << std::endl;
                func_evals += nlsc.GetNumFuncEvals(), jac_evals += nlsc.GetNumJacEvals();
            }

            w.obj(oids[0]).save(
                    to_save  + w.obj(oids[0]).name + "_it=" + std::to_string(it) + ".vtk");
            ObSaveTag(w.obj(oids[0]),
                      to_save + w.obj(oids[0]).name + "_it=" + std::to_string(it) + "_tag_x.txt",
                      "v:x");
            print_iterate_info(it);
        }
        copy_update();
        std::cout << (solverItData[s_id].ns_type == 0 ? result_str_nlsp(nlsp) : result_str_nlsc(nlsc)) << std::endl;
        lgfile << (solverItData[s_id].ns_type == 0 ? result_str_nlsp(nlsp) : result_str_nlsc(nlsc)) << std::endl;

        if (!slvFlag) {
            std::cout << "Convergence is not reached\n Break" << std::endl;
            lgfile << "Convergence is not reached\n Break" << std::endl;
            exit(-1);
        }
    }

    /*process and save result*/
    if (process_result){
        w.obj(oids[0]).save(extr_save + w.obj(oids[0]).name + ".txt");
        ObSaveTag(w.obj(oids[0]), extr_save  + w.obj(oids[0]).name + "_x0.txt", "v:point");
        for (int i = 0; i < oids.size(); ++i) {
            w.obj(oids[i]).save(to_save + w.obj(oids[i]).name + ".vtk");
        }
        double l_free = 0, l_free_init = 0, l_free_flat = 0;
        {
            auto &leaf = w.obj(oids[0]);
            auto init_x = leaf.m_mesh.property_map<V_ind, Point>(init_aniso_crd).first;
            for (auto e: leaf.m_mesh.edges()) {
                auto f = face_around(leaf.m_mesh, e);
                if (f[0].second && f[1].second)
                    continue;
                auto v = vert_around(leaf.m_mesh, e);
                auto vop = vert_opposite(leaf.m_mesh, e);
                if ((leaf.m_boundary[v[0]] & FREE) && (leaf.m_boundary[v[1]] & FREE)) {
                    double loc = sqrt((leaf.m_x[v[0]] - leaf.m_x[v[1]]).squared_length());
                    double loc_init = sqrt((leaf.m_x0[v[0]] - leaf.m_x0[v[1]]).squared_length());
                    double loc_flat = sqrt((init_x[v[0]] - init_x[v[1]]).squared_length());
                    l_free += loc, l_free_init += loc_init, l_free_flat += loc_flat;
                }
            }
        }
        Object3D dup;
        Messurer mes(dup, w.obj(oids[0]), w.obj(oids[1]), w.obj(oids[2]), {0, 0, 1});
        mes.setMargin(1*0.01*mpR);
        auto specialCase = [&w, oids](const Messurer::ObjectShape& from, const Messurer::ObjectShape& to, std::vector<Messurer::CollisionInfo::CPD>& info)->bool{
            int from_i = -1, to_i = -1;
            for (int i = 0; i < 3; ++i){
                if (&w.obj(oids[i]) == from.obj) from_i = i;
                if (&w.obj(oids[i]) == to.obj) to_i = i;
            }
            auto check_cond = [oi = from_i](Point p){
                double tau = atan2(p.y(), p.x());
                if (tau < 0) tau += 2*M_PI;
                if (oi == 0 && tau >= M_PI/3 && tau <= 5*M_PI/3) return true;
                if (oi != 0 && !(tau > M_PI/3*(2* oi - 1)+1e-5 && tau < M_PI/3*(1 + 2*oi)-1e-5)) return true;
                return false;
            };
            bool changed = false;
            for (int v = 0; v < from.mesh.points.size(); ++v){
                if (check_cond(from.mesh.points[v])){
                    Messurer::CollisionInfo::CPD c;
                    c.v = V_ind(v);
                    c.f = F_ind(info[0].f);
                    c.fp = from.mesh.points[v];
                    c.dist2 = 0;
                    info.push_back(c);

                    changed = true;
                }
            }

            return changed;
        };
        mes.computeCollidingBnd(4, specialCase);
        mes.computeMidPlanes();
        for (int i = 0; i < 3; ++i)
            for (int j = i+1; j < 3; ++j){
                std::array<Vector, 3> pl_nrm{Vector(-sqrt(3)/2, 0.5, 0), Vector(-sqrt(3)/2, -0.5, 0), Vector(0, -1, 0)};
                auto n = pl_nrm[i + j - 1];
                mes.m_planes(i, j) = {Messurer::Plane_3(CGAL::ORIGIN,  n), true};
                mes.m_planes(j, i) = {Messurer::Plane_3(CGAL::ORIGIN, -n), true};
            }

        for (int i = 0; i < 3; ++i) {
            auto colfilter = [&mes, &col = mes.m_colMap[i]](F_ind f) {
                auto lbl = col.face_lbl[col.remap[f]];
                bool res = true;
                for (int k = 0; k < 3; ++k)
                    res &= ((lbl & (7 << 3 * k)) > 0);
                return res;
            };
            auto nocolfilter = [&colfilter](F_ind f) { return !colfilter(f); };
            Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_init.stl", "v:point");
            Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + mes.getCusps()[i]->name + "_mes_init_col.stl",
                                      "v:point",colfilter);
            Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + mes.getCusps()[i]->name + "_mes_init_nocol.stl",
                                      "v:point",nocolfilter);
            Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + mes.getCusps()[i]->name + "_mes_flat_col.stl",
                                      init_aniso_crd,colfilter);
            Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + mes.getCusps()[i]->name + "_mes_flat_nocol.stl",
                                      init_aniso_crd,nocolfilter);
            Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_real_nocol.stl",
                                      "v:x", nocolfilter);
            Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_real_col.stl",
                                      "v:x", colfilter);
        }
        mes.computeHalfCoaptScans();
        for (int i = 0; i < 1; ++i)
            for (int j = 0; j < 3; ++j) {
                if (i == j) continue;
                Object3D(mes.m_hcScan(i, j).to_Mesh()).save(
                        extr_save + "half_scan_" + std::to_string(i) + "_" + std::to_string(j) + ".stl");
            }
        auto distr = mes.computeHalfCoaptDistrib(500);
        mes.uniteCollidingCuspHalves();
        for (int i = 0; i < 1; ++i) Object3D(mes.m_coaptScan[i].to_Mesh()).save(extr_save + "unite_scan_" + std::to_string(i) + ".stl");
        auto fdistr = mes.computeCoaptDistribution({1000, 1000, 1000});
        Messurer::saveCoaptDistribCSV(fdistr, extr_save + "full_distrib.csv");
        auto Hc = mes.computeHc();
        auto H = distr.getCoaptH();
        auto coaptStat = mes.computeCoaptStatus();
        bool isClosed = ((coaptStat & (1 << 6)) > 0);
        auto bill = mes.computeBillowing();
        auto areas = mes.computeColArea();
        auto print_results = [&](std::ostream& out){
            out << "Hcoapt = " << H << ", Hcentral = " << Hc << "\n"
                                                                "coaptStatus = " << coaptStat << " isClosed = " << isClosed << "\n";
            bill.print(out);
            out << "Coapt areas: \n\t0:" << areas[0] << "\n\t1:" << areas[1] << "\n\t2:" << areas[2] << "\n";
            out << "L_free = " << l_free << " L_free_init = " << l_free_init << " L_free_flat = " << l_free_flat << "\n";
        };
        print_results(std::cout);
        print_results(lgfile);

        Messurer::saveHalfCoaptDistribCSV(distr, extr_save + "_" + "half_distrib.csv");

        std::ofstream comtable;
        auto s = status(main_save + "comres.csv");
        if (!exists(s)) {
            comtable.open(main_save + "comres.csv");
            comtable << "res_dir; Hcpt; Hcent; Acpt; billowing; L_free;"
                        "scenario; Ht; P; mu; k1; k2; kappa; theta; mesh_h; f_abs_err; maxnits; delta_init; sim_status\n";
        }
        else
            comtable.open(main_save + "comres.csv", std::ios_base::app);

        comtable << "\"" << to_save << "\"" << ";" << H << ";" << Hc << ";" << (areas[0].colArea + areas[1].colArea + areas[2].colArea)/3 << ";"
                 << bill.getBillowing(0) << ";" << l_free << ";"
                 << scenario << ";" << Ht << ";" << P << ";" << mu << ";" << k1 << ";" << k2 << ";" << kappa << ";" << theta << ";"
                 << mesh_h[templ] << ";" << f_err_lim << ";" << maxnit << ";" << delta_relax << ";" << (do_solve ? slvFlag : 2) << std::endl;

        comtable.close();
    }

    if (view){
        view_window->join();
        delete view_window;
    }

    return 0;
}

int ModelLeafBench(int argc, char* argv[]){
    struct ModelLeafParams{
        double R = 15;
        double phi = 2*M_PI/3 - 0.4;
        double Hc = 1.5;
        double Hb = 15;
        double Omega = -40.0/180.0 * M_PI;
    };
    int scenario = 0; // 0 - membrane, 1 - shape, 2 - shape + BC
    double Ht = 0.4;
    double mu = 9e2; //kPa
    double theta = M_PI_4, k1 = 1.1e3 /*kPa*/, k2 = 1e3, kappa = 0.15;//.29;
    double P = 90.0_mmHg / 1e3; //kPa
    ModelLeafParams mp{15, 2*M_PI/3 - 0.4, 1.5, 15, -40.0/180.0 * M_PI};
    double mesh_h = 0.5;
    double err = 1e-10, diverge_lim = 1e10;
    int maxits = 500;
    bool regenerate = false;
    string init_config_prefix = "_4_2";
    bool only_process = false;
    int set_relax_its = -1, set_pfreq = -1;
    std::vector<double> pressures = {};

    string main_save = "../result/ModelLeafBench/", to_save, extr_dir = "extra/";
    {
        std::ostringstream oss;
        auto tt = std::time(nullptr);
        auto tm = *std::localtime(&tt);
        oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
        to_save = "default/";//oss.str() + "/";//"default/";//TODO: repair this
    }
    to_save = main_save + to_save;
    string to_gen = "../generated/";
    string prefix = "";

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << "\n";
            std::cout << " Command line options: " << "\n";
            std::cout << "  -sc, --scenario     <Choosing of the running scenario: 0 - membrane, 1 - shape, 2 - shape + BC>" << "\n";
            std::cout << "  -ht, --thickness    <Leafs thickness>" << "\n";
            std::cout << "  -pr, --pressure     <Diastolic pressure>" << "\n";
            std::cout << "  -mu, --mu           <Shear modulus>" << "\n";
            std::cout << "  -a , --aniso        <Parameters of anisotropy: k1, k2, kappa, theta>" << "\n";
            std::cout << "  -mh, --mesh_h       <Step of generated mesh>" << "\n";
            std::cout << "  -t , --target       <Directory to save results>" << "\n";
            std::cout << "  -tm, --maindir      <Common directory to save results from different simulations>" << "\n";
            std::cout << "  -p , --prefix       <Prefix to file name>" << "\n";
            std::cout << "  -gd, --gendir       <Path to save generated codes>" << "\n";
            std::cout << "  -sp, --slvstop      <Solver stop criterion parameters: err, maxits, diverge_lim>" << "\n";
            std::cout << "  -c , --case         <Set default parameters for the case>" << "\n";
            std::cout << "  -op, --only_process <Only process read data without calculations>" << "\n";
            std::cout << "  -r,  --regenerate   <Regenerate symbolic functions>" << "\n";
            std::cout << "  -its,--setrelaxits  <Set custom relax its and print frequency>" << "\n";
            std::cout << "  -ipr,--initpressure <Set initial relaxation iteration with lower pressure>" << "\n";
            std::cout << "  -vp, --valve        <Valve parameters: R, phi, Hc, Hb, Omega>" << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-ipr") == 0 || strcmp(argv[i], "--initpressure") == 0) {
            bool next = true;
            while (next && i + 1 < argc) {
                try {
                    if (i + 1 < argc) {
                        double p = stod(argv[++i]);
                        pressures.push_back(p * 1.0_mmHg / 1.0e3);
                    }
                } catch (std::invalid_argument &e) {
                    i--;
                    next = false;
                } catch (std::exception &e) {
                    std::cout << "Waited number of relax iterations value but error happens: '" << e.what() << "'\n";
                    --i;
                    next = false;
                }
            }
            continue;
        }
        if (strcmp(argv[i], "-its") == 0 || strcmp(argv[i], "--setrelaxits") == 0) {
            try {
                if (i + 1 < argc) set_relax_its = stoi(argv[++i]);
                if (i + 1 < argc) set_pfreq = stoi(argv[++i]);
            }  catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited number of relax iterations value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--case") == 0) {
            try {
                if (i + 1 < argc) {
                    int cs = stoi(argv[++i]);
                    if (cs < 0 || cs >= 96){
                        throw std::runtime_error("Wrong value of case, case should be integer in [0, 95]");
                    }
                    int cluster = (cs / 48); cs %= 48;
                    int tetacs = cs % 3; cs /= 3;
                    int kcs = cs % 4; cs /= 4;
                    int mucs = 2*cluster + cs % 2; cs /= 2;
                    int sc_cs = cs % 2;
                    scenario = (sc_cs == 1) ? 2 : 0;
                    switch (tetacs) {
                        case 0: theta = 0; break;
                        case 1: theta = M_PI_4; break;
                        case 2: theta = M_PI_2; break;
                    }
                    switch (kcs) {
                        case 0: kappa = 0; break;
                        case 1: kappa = 0.2; break;
                        case 2: kappa = 0.29; break;
                        case 3: kappa = 1.0/3; break;
                    }
                    switch (mucs) {
                        case 0: mu = 9e1; break;
                        case 1: mu = 9e2; break;
                        case 2: mu = 3e2; break;
                        case 3: mu = 3e3; break;
                    }
                    to_save = main_save + "case" + argv[i] + "/";
                }
            } catch (std::exception& e){
                std::cout << "Waited number of scenario but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-sc") == 0 || strcmp(argv[i], "--scenario") == 0) {
            try {
                if (i + 1 < argc) {
                    scenario = stoi(argv[++i]);
                    if (scenario < 0 || scenario > 3){
                        throw std::runtime_error("Wrong value of scenario, you should choose one of the following:\n\t0 - membrane\n\t1 - shape\n\t2 - shape + BC\n");
                    }
                }
            } catch (std::exception& e){
                std::cout << "Waited number of scenario but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-ht") == 0 || strcmp(argv[i], "--thickness") == 0) {
            try {
                if (i + 1 < argc)
                    Ht = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited thickness value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--regenerate") == 0) {
            regenerate = true;
            continue;
        }
        if (strcmp(argv[i], "-op") == 0 || strcmp(argv[i], "--only_process") == 0) {
            only_process = true;
            continue;
        }
        if (strcmp(argv[i], "-pr") == 0 || strcmp(argv[i], "--pressure") == 0) {
            try {
                if (i + 1 < argc)
                    P = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited thickness value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-mu") == 0 || strcmp(argv[i], "--mu") == 0) {
            try {
                if (i + 1 < argc)
                    mu = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited mu value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--aniso") == 0) {
            int j = i;
            const char* v[] = {"k1", "k2", "kappa", "theta"};
            try {
                if (i + 1 < argc) k1 = stof(argv[++i]);
                if (i + 1 < argc) k2 = stof(argv[++i]);
                if (i + 1 < argc) kappa = stof(argv[++i]);
                if (i + 1 < argc) theta = stof(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-vp") == 0 || strcmp(argv[i], "--valve") == 0) {
            int j = i;
            const char* v[] = {"R", "phi", "Hc", "Hb", "Omega"};
            try {
                if (i + 1 < argc) mp.R = stof(argv[++i]);
                if (i + 1 < argc) mp.phi = stof(argv[++i]);
                if (i + 1 < argc) mp.Hc = stof(argv[++i]);
                if (i + 1 < argc) mp.Hb = stof(argv[++i]);
                if (i + 1 < argc) mp.Omega = stof(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-mh") == 0 || strcmp(argv[i], "--mesh_h") == 0) {
            try {
                if (i + 1 < argc)
                    mesh_h = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited mesh size value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-sp") == 0 || strcmp(argv[i], "--slvstop") == 0) {
            int j = i;
            const char* v[] = {"err", "maxits", "diverge_lim"};
            try {
                if (i + 1 < argc) err = stof(argv[++i]);
                if (i + 1 < argc) maxits = stoi(argv[++i]);
                if (i + 1 < argc) diverge_lim = stoi(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }

        if (strcmp(argv[i], "-gd") == 0 || strcmp(argv[i], "--gendir") == 0) {
            if (i + 1 < argc) to_gen = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0) {
            if (i + 1 < argc) to_save = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-tm") == 0 || strcmp(argv[i], "--maindir") == 0) {
            if (i + 1 < argc) main_save = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) {
            if (i + 1 < argc) to_save = argv[++i];
            continue;
        }
    }

    auto print_input = [&](std::ostream& out) -> std::ostream& {
        return out << "Input parameters: " <<
                   "\n\tscenario = " << scenario <<
                   "\n\tHt = " << Ht <<
                   "\n\tP = "  << P <<
                   "\n\tmu = " << mu <<
                   "\n\tk1 = " << k1 << ", k2 = " << k2 << ", kappa = " << kappa << ", theta = " << theta <<
                   "\n\tmp = {"
                   " R = "    << mp.R <<
                   ", phi = "    << mp.phi <<
                   ", Hc = "     << mp.Hc <<
                   ", Hb = "    << mp.Hb <<
                   ", Omega = " << mp.Omega <<
                   "}" <<
                   "\n\tmesh_h = " << mesh_h <<
                   "\n\terr = " << err << ", maxits = " << maxits << ", diverge_lim = " << diverge_lim <<
                   "\n\tmain_dir = \"" << main_save << "\", save_dir = \"" << to_save << "\", gen_dir = \"" << to_gen << "\"" <<
                   "\n\tprefix = \"" << prefix << "\"" << std::endl;
    };

    print_input(std::cout);
    {
        auto make_dir = [](std::string save) {
            path path = save;
            auto s = status(path);
            if (!status_known(s) || s.type() != file_type::directory)
                create_directory(path);
        };
        make_dir(main_save);
        make_dir(to_save);
        make_dir(to_save + extr_dir);
        make_dir(to_gen);
    }

    std::string extr_save = to_save + extr_dir + prefix;
    std::ofstream lgfile(to_save + "log.txt", std::ios_base::app);
    print_input(lgfile);

    Object3D leaf;
    bool new_mesh = false;
    string init_aniso_crd = "v:point";
    if (!new_mesh) {
        bool regen_mesh = false;
        if (regen_mesh) {
            computeInitialMesh(ComputeInitMesh().set_args(argc, argv));
        }
        leaf.name = "ModelLeaf";
        leaf.read("../result/ModelLeafBench/template/leaf_init" + init_config_prefix + ".txt");
        init_aniso_crd = "v:initial";
        auto init_x = leaf.m_mesh.add_property_map<V_ind, Point>(init_aniso_crd).first;
        ObReadTag(leaf, "../result/ModelLeafBench/template/leaf_init_x0" + init_config_prefix + ".txt",
                  init_aniso_crd);//"v:point");//"v:initial");
        for (auto v: leaf.m_mesh.vertices()) leaf.m_next_x[v] = leaf.m_x[v], leaf.m_F[v] = CGAL::NULL_VECTOR;
    } else {
        leaf = computeInitialLeaf(mp.R, mp.phi, mp.Hc, mp.Hb, mp.Omega, mesh_h);
        leaf.name = "ModelLeaf";
    }
    auto aniso_f = leaf.m_mesh.add_property_map<F_ind, std::array<DReal , 3>>("f:aniso_f").first;
    auto aniso_s = leaf.m_mesh.add_property_map<F_ind, std::array<DReal , 3>>("f:aniso_s").first;
    auto vaniso_f = leaf.m_mesh.add_property_map<V_ind, Vector>("v:aniso_f").first;
    auto vaniso_s = leaf.m_mesh.add_property_map<V_ind, Vector>("v:aniso_s").first;
    if (init_aniso_crd == "v:point"){
        for (auto f: leaf.m_mesh.faces()){
            aniso_f[f] = std::array<DReal , 3>{cos(theta),  sin(theta), 0};
            aniso_s[f] = std::array<DReal , 3>{cos(theta), -sin(theta), 0};
        }
    } else {
        auto _init_x = leaf.m_mesh.property_map<V_ind, Point>(init_aniso_crd);
        if (!_init_x.second) throw std::runtime_error("Tag \"" + init_aniso_crd + "\" not found");
        auto& init_x = _init_x.first;
        for (auto f: leaf.m_mesh.faces()) {
            auto& o = leaf;
            auto v = vert_around(leaf.m_mesh, f);
            Eigen::Matrix3d Q, P;
            Q <<    o.m_x0[v[0]][0], o.m_x0[v[1]][0], o.m_x0[v[2]][0],
                    o.m_x0[v[0]][1], o.m_x0[v[1]][1], o.m_x0[v[2]][1],
                    o.m_x0[v[0]][2], o.m_x0[v[1]][2], o.m_x0[v[2]][2];
            P <<    init_x[v[0]][0], init_x[v[1]][0], init_x[v[2]][0],
                    init_x[v[0]][1], init_x[v[1]][1], init_x[v[2]][1],
                    init_x[v[0]][2], init_x[v[1]][2], init_x[v[2]][2];
            Eigen::Vector3d n0 = (P.col(1) - P.col(0)).cross(P.col(2) - P.col(0));
            double Ap = n0.stableNorm() / 2;
            n0 /= (2*Ap);
            Eigen::Matrix3d Dp;
            for (int i = 0; i < 3; ++i) {
                Dp.col(i) = n0.cross(P.col((i + 2) % 3) - P.col((i + 1) % 3));
            }
            Dp /= (2*Ap);
            Eigen::Vector3d Mfp(cos(theta),  sin(theta), 0), Msp(cos(theta),  -sin(theta), 0);
            auto F = Q * Dp.transpose();
            Eigen::Vector3d Mf = F * Mfp, Ms = F * Msp;
            Mf.normalize(), Ms.normalize();
            aniso_f[f] = std::array<DReal , 3>{Mf[0],  Mf[1], Mf[2]};
            aniso_s[f] = std::array<DReal , 3>{Ms[0],  Ms[1], Ms[2]};
        }
    }
    for (auto v: leaf.m_mesh.vertices()){
        vaniso_f[v] = vaniso_s[v] = CGAL::NULL_VECTOR;
        auto fs = face_around(leaf.m_mesh, v);
        for (auto f: fs){
            Vector ff(aniso_f[f][0], aniso_f[f][1], aniso_f[f][2]), ss(aniso_s[f][0], aniso_s[f][1], aniso_s[f][2]);
            vaniso_f[v] += ff / fs.size();
            vaniso_s[v] += ss / fs.size();
        }

    }

    static const int FIXED = 2, FREE = 1|4|8;
    DirichletCond dc = [](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & FIXED) return 0;
        //if (obj.m_boundary[v] & FREE) return 3;
        return 7;
    };
    leaf.setDirichletCond(dc);
    auto bdata = computeBndDataOnCilAligned(leaf,
                                            [](Object3D& leaf, E_ind e, V_ind* v, V_ind vop){
                                                return (leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED) && !(leaf.m_boundary[vop] & FREE);
                                            });
    if (!new_mesh) {
        leaf.apply_updaters();
        for (auto &d: bdata) {
            auto e = d.first;
            auto f = face_around(leaf.m_mesh, e)[0].first;
            auto b = Eigen::Map<Eigen::Vector3d>(d.second.normal);
            Eigen::Vector3d n(leaf.m_normal[f][0], leaf.m_normal[f][1], leaf.m_normal[f][2]);
            b = b - b.dot(n) * n;
            b.normalize();
        }
    }
    auto bdit = addBdataTag(leaf, bdata);

    int ait = -1;
    path path = to_save;
    auto _s = status(path);
    if (_s.type() == file_type::directory){
        for (auto& p: directory_iterator(to_save)){
            if (p.status().type() == file_type::regular){
                string name = p.path().filename();
                auto pos = name.find_last_of('=');
                if (pos != name.npos && name.size() - pos > 1) {
                    string nitstr = name.substr(name.find_last_of('=') + 1);
                    nitstr = nitstr.substr(0, nitstr.find_first_of('_'));
                    int nit = stoi(nitstr);
                    if (nit > ait) ait = nit;
                }
            }
        }
        if (ait >= 0){
            ObReadTag(leaf, to_save + "ModelLeaf_it="+ to_string(ait) +"_tag_x.txt","v:x");
        }
        ait++;
    }

    /*initialize forces*/
    SimplePressureLoad Pf(P);
    HGO_GSTModel Ef(Ht, mu, k1, k2, kappa, theta, to_gen, regenerate, 15);
    Ef.prepareJacobianFunction(regenerate); Ef.f.setVerbosity(2);
    BendingForce bf(Ef.f, regenerate);
    bf.prepareJacobianFunction(regenerate); bf.setVerbosity(2);
    SpurnPlane bottom_constraint(2*P, 0.01*mp.Hb,  SpurnPlane::Plane_3(Point(0, 0, -mp.Hb*1.4), Vector(0, 0, 1)));
    SpurnPlane left_constraint(1.1*P, 0.01*mp.R, SpurnPlane::Plane_3(CGAL::ORIGIN, Vector(-sin(-M_PI/3), cos(-M_PI/3), 0.0)));
    SpurnPlane right_constraint(1.1*P, 0.01*mp.R, SpurnPlane::Plane_3(CGAL::ORIGIN, -Vector(-sin(M_PI/3), cos(M_PI/3), 0.0)));
    std::map<E_ind, Vector> cdat; for(const auto& i: bdata)
        if (i.second.btype == BendingForce::BoundaryData::CLAMPED)
            cdat.insert({i.first, Vector{i.second.normal[0], i.second.normal[1], i.second.normal[2]}});
    ClampedBndPenalty clamped(0, cdat);
    /*set world*/
    World w;
    auto oid = w.addObject3D(std::move(leaf), 1);
    /*set symmetrical leafs*/
    Object3D leaf1 = w.obj(oid), leaf2 = w.obj(oid);
    leaf1.name = "ModelLeaf1"; leaf2.name = "ModelLeaf2";
    auto oid1 = w.addObject3D(std::move(leaf1), 1);
    auto oid2 = w.addObject3D(std::move(leaf2), 1);
    std::array<ObjectID, 3> oids {oid, oid1, oid2};
    auto copy_update = [&w, oids](){
        auto& mob = w.obj(oids[0]);
        for (int i = 1; i < 3; ++i) {
            Object3D& obj = w.obj(oids[i]);
            Eigen::AngleAxis<double> q(2.0 * M_PI / 3.0 * i, Eigen::Vector3d(0, 0, 1));
            for (auto v: obj.m_mesh.vertices()){
                Eigen::Vector3d p(mob.m_x[v][0], mob.m_x[v][1], mob.m_x[v][2]);
                Eigen::Vector3d pc = q * p;
                obj.m_next_x[v] = obj.m_x[v] = Point(pc[0], pc[1], pc[2]);
            }
        }
    };
    copy_update();
    w.objs()[oid1].second = 1;
    w.objs()[oid2].second = 1;

    bool slvFlag = false;
    if (!only_process) {
        /*set forces*/
        Force_ID fpi, fei, fbti, flti, frti, fbi, fcli;
        fpi = w.addForce(std::move(Pf), oid);
        fei = w.addForce(std::move(Ef), oid);
        //fbti = w.addForce(bottom_constraint, oid);
        flti = w.addForce(left_constraint, oid);
        frti = w.addForce(right_constraint, oid);
        if (scenario == 1 || scenario == 2) {
            fbi = w.addForce(std::move(bf), oid);
        }
        if (scenario == 2) {
            if (!new_mesh)
                w.obj(oid).m_forces[fbi].target<BendingForce>()->set_boundary_data(w.obj(oid), bdata);
            else
                fcli = w.addForce(Force(clamped), oid);
        }
        //initial saving
        if (ait == 0)
            w.obj(oids[0]).save(to_save + prefix + w.obj(oids[0]).name + "_init_state" + ".vtk", w.obj(oids[0]).m_x0);

        /*set solver*/
        NSWorldWrapper nsww(w);
        NLProblem nlp = NSWorldWrapper::makeProblem(nsww);
        NonLinearSolverKinsol nlsp(nlp);
        LinearSolver ls("inner_mptiluc");
        ls.setVerbosityLevel(1);
        ls.SetParameterReal("relative_tolerance", 1e-20);
        double droptol = 8e-3;//2e-5;
        ls.SetParameterReal("drop_tolerance", droptol);
        ls.SetParameterReal("reuse_tolerance", 0.1 * droptol);
        nlsp.SetLinearSolver(&ls);
        static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
        nlsp.SetVerbosityLevel(1);
        nlsp.SetInitialGuess(nsww.getCurrentX());
        nlsp.SetFuncNormTol(1e-5);
        nlsp.SetMaxNewtonStep(mp.R / 10 * sqrt(nlp.m_dofs));
        nlsp.SetParameterInt("MaxSetupCalls", 1);
        nlsp.SetParameterInt("MaxSubSetupCalls", 1);
        nlsp.SetScaledStepTol(1e-7);
        nlsp.SetMaxBetaFails(40);
        World3d::Timer sym_time, com_time;
        com_time.reset();
        nlsp.SetInfoHandlerFn([&time = sym_time, &com_time, &lgfile](const char *module, const char *function, char *msg) {
            std::cout << "[" << module << "] " << function << "\n   " << msg << "\n" << "time = " << time.elapsed()
                      << " com_time = " << com_time.elapsed() << std::endl;
            lgfile << "[" << module << "] " << function << "\n   " << msg << "\n" << "time = " << time.elapsed()
                      << " com_time = " << com_time.elapsed() << std::endl;
        });
        auto monitor_fcn = [&nsww, &copy_update](const double *x) {
            nsww.setX(x);
            copy_update();
            nsww.RenderIteration();
            return 0;
        };
        nlsp.setInterIterationMonitorFunc(monitor_fcn);
        auto result_str = [](NonLinearSolverKinsol &nlsp) {
            std::stringstream sstr;
            sstr << "\t#NonLinIts = " << nlsp.GetNumNolinSolvIters() << " Residual = " << nlsp.GetResidualNorm()
                 << "\n\t#linIts = " << nlsp.GetNumLinIters() << " #funcEvals = " << nlsp.GetNumFuncEvals()
                 << " #jacEvals = " << nlsp.GetNumJacEvals()
                 << "\n\t#convFails = " << nlsp.GetNumLinConvFails() << " #betaCondFails = "
                 << nlsp.GetNumBetaCondFails() << " #backtrackOps = " << nlsp.GetNumBacktrackOps() << "\n";
            sstr << "Reason: " << nlsp.GetReason();
            return sstr.str();
        };
        /*solve problem*/
//    thread t(start_debug_gui, argc, argv);
//    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());

        if (scenario == 2 && new_mesh) {
            w.obj(oid).m_forces[fcli].target<ClampedBndPenalty>()->setStiffnes(0);
            std::cout << "SET STIFFNESS = " << 0 << std::endl;
            lgfile << "SET STIFFNESS = " << 0 << std::endl;
        }
        auto relaxate_solver = [Ht, &w, &err, &diverge_lim, &sym_time, &lgfile, oid](double &delta, int nits,
                                                                                     int pfreq) {
            auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &nits, &pfreq, &sym_time,
                                                         &lgfile,
                                                         true);
            StaticForceApplierP sfa(delta);
            sfa.addDamper(Damper(0.01 * Ht, 0.8, 15, 1000000).setInitDxScale(3.0)
                                  .setReInitDxParams(-1, -1).setMinDeltaFactor(1e-5));
            w.setForceApplier(sfa);
            sym_time.reset();
            w.Simulation(stopCond);
            double factor = w.obj(
                    oid)._next_x_updater.target<StaticForceApplierP>()->_dampers[0].target<Damper>()->getUsedDeltaFactor();
            delta *= factor;
            if (1 - factor > 1e-2) std::cout << "New relaxation delta = " << delta << std::endl;
        };
        int func_evals = 0, jac_evals = 0;
        auto print_iterate_info = [&fe = func_evals, &je = jac_evals, ait, &com_time, &lgfile](int it){
            std::cout << "ITERATE INFO: it: " << it + ait << " com_time = " << com_time.elapsed() <<
                          " accum_func_evals = " << fe << " accum_jac_evals = " << je << std::endl;
            lgfile << "ITERATE INFO: it: " << it + ait << " com_time = " << com_time.elapsed() <<
                      " accum_func_evals = " << fe << " accum_jac_evals = " << je << std::endl;
        };
        double delta = 1e-4;
        if (scenario == 1 || scenario == 2) {
            if (mu > 3e2) delta = 1e-5;
            slvFlag = false;

            if (!pressures.empty()) {
                for (int i = 0, cnt = pressures.size() - 1; i <= cnt; ++i) {
                    w.obj(oid).m_forces[fpi].target<SimplePressureLoad>()->setPressure(
                            [pf = pressures[i]](Object3D &, F_ind) { return pf; });
                    std::cout << "MSG: Set pressure to " << pressures[i] << " mmHg" << std::endl;
                    if (i == 0) {
                        nlsp.SetNumMaxIters(maxits);
                        nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
                        sym_time.reset();
                        slvFlag = nlsp.Solve();
                        std::cout << result_str(nlsp) << std::endl;
                    } else {
                        nlsp.SetNumMaxIters(maxits);
                        nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
                        sym_time.reset();
                        slvFlag = nlsp.Solve();
                        std::cout << result_str(nlsp) << std::endl;
                        if (!slvFlag) {
                            relaxate_solver(delta, 10000, 200);
                            if (nsww.compCorrentForceNorm() < nlsp.GetResidualNorm()) {
                                nlsp.SetInitialGuess(nsww.getCurrentX());
                                slvFlag = nlsp.Solve();
                            }
                        }
                    }
                }
                w.obj(oid).m_forces[fpi].target<SimplePressureLoad>()->setPressure(
                        [pf = P](Object3D &, F_ind) { return pf; });
            }

            int main_it = 1000;
            int relax_maxits = 5000, newt_maxits = 40;
            int pfreq = 200;
            bool force_relax_params = false;
            auto setRelaxItsPrm = [&force_relax_params, &relax_maxits, &pfreq, &cits = set_relax_its, &cpfreq = set_pfreq](int maxits, int freq){
                if (!force_relax_params) {
                    relax_maxits = (cits < 0) ? maxits : cits;
                    pfreq = (cpfreq < 0) ? freq : cpfreq;
                } else {
                    relax_maxits = maxits;
                    pfreq = freq;
                }
            };
            setRelaxItsPrm(5000, 200);
            for (int it = 0; it < main_it && !slvFlag; ++it) {
                //nlsp.SetScaledStepTol(1e-5);
                if (it < 2) nlsp.SetNumMaxIters(maxits);
                else nlsp.SetNumMaxIters(newt_maxits);
                nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
                sym_time.reset();
                if (it + ait != 0) {
                    slvFlag = nlsp.Solve();
                    lgfile << result_str(nlsp) << std::endl;
                    func_evals += nlsp.GetNumFuncEvals(), jac_evals += nlsp.GetNumJacEvals();
                }


                double newt_resid = nlsp.GetResidualNorm();
                bool make_relaxation = !slvFlag;
                force_relax_params = false;
                while (make_relaxation) {
                    make_relaxation = false;
                    if (it != 0 && !force_relax_params) {
                        if (nlsp.GetReasonFlag() != KIN_MAXITER_REACHED) {
                            setRelaxItsPrm(1000, 50);
                        } else {
                            setRelaxItsPrm(1000, 50);
                        }
                    }
                    relaxate_solver(delta, relax_maxits, pfreq);
                    func_evals += relax_maxits;
                    double relax_resid = nsww.compCorrentForceNorm();
                    if (relax_resid < newt_resid || nlsp.GetReasonFlag() != KIN_MAXITER_REACHED) {
                        if (it+ait == 0 || relax_resid / newt_resid < 2e1)
                            nlsp.SetInitialGuess(nsww.getCurrentX());
                        else {
                            nsww.setX(N_VGetArrayPointer(nlsp.GetVectorX()));
                            force_relax_params = true;
                            make_relaxation = true;
                            setRelaxItsPrm(relax_maxits/2, pfreq);
                            if (relax_maxits <= 15) {
                                setRelaxItsPrm(15, pfreq);
                                delta *= 0.8;
                            }
                        }
                    }
                }
                w.obj(oids[0]).save(
                        to_save + prefix + w.obj(oids[0]).name + "_it=" + std::to_string(it + ait) + ".vtk");
                ObSaveTag(w.obj(oids[0]),
                          to_save + prefix + w.obj(oids[0]).name + "_it=" + std::to_string(it + ait) + "_tag_x.txt",
                          "v:x");
                print_iterate_info(it);
            }
            copy_update();
            std::cout << result_str(nlsp) << std::endl;
            lgfile << result_str(nlsp) << std::endl;

            if (new_mesh && scenario == 2) {
                double stiff = 0.5;
                for (int i = 1; i < 1000; ++i) {
                    if (i <= 8) stiff *= 2;
                    else stiff += 100;
                    w.obj(oid).m_forces[fcli].target<ClampedBndPenalty>()->setStiffnes(stiff);
                    std::cout << "SET STIFFNESS = " << stiff << std::endl;
                    lgfile << "SET STIFFNESS = " << stiff << std::endl;
                    nlsp.SetNumMaxIters(maxits);
                    sym_time.reset();
                    slvFlag = nlsp.Solve();
                    func_evals += nlsp.GetNumLinFuncEvals(), jac_evals += nlsp.GetNumJacEvals();
                    std::cout << result_str(nlsp) << std::endl;
                    lgfile << result_str(nlsp) << std::endl;
                    print_iterate_info(i);
                }
            }
        }
        if (scenario == 0) {
            slvFlag = false;
            delta = 1e-4;
            int main_it = 100;
            int relax_maxits = 10000, newt_maxits = 400;
            int pfreq = 200;
            auto setRelaxItsPrm = [&relax_maxits, &pfreq, &cits = set_relax_its, &cpfreq = set_pfreq](int maxits, int freq){
                relax_maxits = (cits < 0) ? maxits : cits;
                pfreq = (cpfreq < 0) ? freq : cpfreq;
            };
            setRelaxItsPrm(10000, 200);
            for (int it = 0; it < main_it && !slvFlag; ++it) {
                nlsp.SetScaledStepTol(1e-5);
                if (it == 0) nlsp.SetNumMaxIters(maxits);
                else nlsp.SetNumMaxIters(newt_maxits);
                nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
                sym_time.reset();
                if (it + ait != 0) {
                    slvFlag = nlsp.Solve();
                    lgfile << result_str(nlsp) << std::endl;
                    func_evals += nlsp.GetNumFuncEvals(), jac_evals += nlsp.GetNumJacEvals();
                }

                double newt_resid = nlsp.GetResidualNorm();
                if (!slvFlag) {
                    auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &relax_maxits, &pfreq, &sym_time,
                                                                 &lgfile,
                                                                 true);
                    StaticForceApplierP sfa(delta);
                    sfa.addDamper(Damper(0.01 * Ht, 0.8, 30, 1000000).setInitDxScale(3.0)
                                          .setReInitDxParams(-1, -1));
                    w.setForceApplier(sfa);
                    sym_time.reset();
                    w.Simulation(stopCond);
                    double factor = w.obj(
                            oid)._next_x_updater.target<StaticForceApplierP>()->_dampers[0].target<Damper>()->getUsedDeltaFactor();
                    delta *= factor;
                    func_evals += relax_maxits;
                    if (1 - factor > 1e-2) std::cout << "New relaxation delta = " << delta << std::endl;

                    double relax_resid = nsww.compCorrentForceNorm();
                    nlsp.SetInitialGuess(nsww.getCurrentX());
                }
                w.obj(oids[0]).save(
                        to_save + prefix + w.obj(oids[0]).name + "_it=" + std::to_string(it + ait) + ".vtk");
                ObSaveTag(w.obj(oids[0]),
                          to_save + prefix + w.obj(oids[0]).name + "_it=" + std::to_string(it + ait) + "_tag_x.txt",
                          "v:x");
                print_iterate_info(it);
            }
            copy_update();
            std::cout << result_str(nlsp) << std::endl;
            lgfile << result_str(nlsp) << std::endl;
        }

        if (!slvFlag) {
            std::cout << "Convergence is not reached\n Break" << std::endl;
            lgfile << "Convergence is not reached\n Break" << std::endl;
            exit(-1);
        }

//        t.join();
    }

    /*process and save result*/
    w.obj(oids[0]).save(extr_save + w.obj(oids[0]).name + ".txt");
    ObSaveTag(w.obj(oids[0]), extr_save  + w.obj(oids[0]).name + "_x0.txt", "v:point");
    for (int i = 0; i < oids.size(); ++i) {
        w.obj(oids[i]).save(to_save + prefix + w.obj(oids[i]).name + ".vtk");
    }
    Object3D dup;
    Messurer mes(dup, w.obj(oids[0]), w.obj(oids[1]), w.obj(oids[2]), {0, 0, 1});
    mes.setMargin(1*0.01*mp.R);
    //std::function<bool(const ObjectShape& from, const ObjectShape& to, std::vector<CollisionInfo::CPD>& info)>
    auto specialCase = [&w, oids](const Messurer::ObjectShape& from, const Messurer::ObjectShape& to, std::vector<Messurer::CollisionInfo::CPD>& info)->bool{
        int from_i = -1, to_i = -1;
        for (int i = 0; i < 3; ++i){
            if (&w.obj(oids[i]) == from.obj) from_i = i;
            if (&w.obj(oids[i]) == to.obj) to_i = i;
        }
        auto check_cond = [oi = from_i](Point p){
            double tau = atan2(p.y(), p.x());
            if (tau < 0) tau += 2*M_PI;
            if (oi == 0 && tau >= M_PI/3 && tau <= 5*M_PI/3) return true;
            if (oi != 0 && !(tau > M_PI/3*(2* oi - 1)+1e-5 && tau < M_PI/3*(1 + 2*oi)-1e-5)) return true;
            return false;
        };
        bool changed = false;
        for (int v = 0; v < from.mesh.points.size(); ++v){
            if (check_cond(from.mesh.points[v])){
                Messurer::CollisionInfo::CPD c;
                c.v = V_ind(v);
                c.f = F_ind(info[0].f);
                c.fp = from.mesh.points[v];
                c.dist2 = 0;
                info.push_back(c);

                changed = true;
            }
        }

        return changed;
    };
    mes.computeCollidingBnd(4, specialCase);//specialCase
    mes.computeMidPlanes();
    for (int i = 0; i < 3; ++i)
    for (int j = i+1; j < 3; ++j){
        std::array<Vector, 3> pl_nrm{Vector(-sqrt(3)/2, 0.5, 0), Vector(-sqrt(3)/2, -0.5, 0), Vector(0, -1, 0)};
        auto n = pl_nrm[i + j - 1];
        mes.m_planes(i, j) = {Messurer::Plane_3(CGAL::ORIGIN,  n), true};
        mes.m_planes(j, i) = {Messurer::Plane_3(CGAL::ORIGIN, -n), true};
    }

    for (int i = 0; i < 3; ++i) {
        auto colfilter = [&mes, &col = mes.m_colMap[i]](F_ind f) {
            auto lbl = col.face_lbl[col.remap[f]];
            bool res = true;
            for (int k = 0; k < 3; ++k)
                res &= ((lbl & (7 << 3 * k)) > 0);
            return res;
        };
        auto nocolfilter = [&colfilter](F_ind f) { return !colfilter(f); };
        Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_init.stl", "v:point");
        Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + prefix + mes.getCusps()[i]->name + "_mes_init_col.stl", "v:point",
                                  colfilter);
        Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + prefix + mes.getCusps()[i]->name + "_mes_init_nocol.stl", "v:point",
                                  nocolfilter);
        Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_real_nocol.stl", "v:x", nocolfilter);
        Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_real_col.stl", "v:x", colfilter);
    }
    mes.computeHalfCoaptScans();
    for (int i = 0; i < 1; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            Object3D(mes.m_hcScan(i, j).to_Mesh()).save(
                    extr_save + "half_scan_" + std::to_string(i) + "_" + std::to_string(j) + ".stl");
        }
    auto distr = mes.computeHalfCoaptDistrib(500);
    mes.uniteCollidingCuspHalves();
    for (int i = 0; i < 1; ++i) Object3D(mes.m_coaptScan[i].to_Mesh()).save(extr_save + "unite_scan_" + std::to_string(i) + ".stl");
    auto fdistr = mes.computeCoaptDistribution({1000, 1000, 1000});
    Messurer::saveCoaptDistribCSV(fdistr, extr_save + "full_distrib.csv");
    auto Hc = mes.computeHc();
    auto H = distr.getCoaptH();
    auto coaptStat = mes.computeCoaptStatus();
    bool isClosed = ((coaptStat & (1 << 6)) > 0);
    auto bill = mes.computeBillowing();
    auto areas = mes.computeColArea();
    auto print_results = [&](std::ostream& out){
        out << "Hcoapt = " << H << ", Hcentral = " << Hc << "\n"
                                                            "coaptStatus = " << coaptStat << " isClosed = " << isClosed << "\n";
        bill.print(out);
        out << "Coapt areas: \n\t0:" << areas[0] << "\n\t1:" << areas[1] << "\n\t2:" << areas[2] << "\n";
    };
    print_results(std::cout);
    print_results(lgfile);

    Messurer::saveHalfCoaptDistribCSV(distr, extr_save + "_" + "half_distrib.csv");

    std::ofstream comtable;
    auto s = status(main_save + "comres.csv");
    if (!exists(s)) {
        comtable.open(main_save + "comres.csv");
        comtable << "res_dir; prefix; Hcpt; Hcent; Acpt;"
                    "scenario; Ht; P; mu; k1; k2; kappa; theta; mp_R; mp_phi; mp_Hc; mp_Hb; mp_Omega; mesh_h; err; maxits; sim_status\n";
    }
    else
        comtable.open(main_save + "comres.csv", std::ios_base::app);

    comtable << "\"" << to_save << "\"" << ";" << "\"" << prefix << "\"" << ";" << H << ";" << Hc << ";" << (areas[0].colArea + areas[1].colArea + areas[2].colArea)/3 << ";"
             << scenario << ";" << Ht << ";" << P << ";" << mu << ";" << k1 << ";" << k2 << ";" << kappa << ";" << theta << ";" << mp.R << ";" << mp.phi
             << ";" << mp.Hc << ";" << mp.Hb << ";" << mp.Omega << ";" << mesh_h << ";" << err << ";" << maxits << ";" << (only_process ? 2 : slvFlag) << std::endl;

    comtable.close();

    return 0;
}

int VanLoonBench(int argc, char* argv[]){
    return computeInitialMesh(ComputeInitMesh().set_args(argc, argv));
    return SaveDisplaceMap();
    return ModelLeafBench1(argc, argv);
    return ModelLeafBench(argc, argv);
//    auto am = generate_HalfEllipseWithRectangle(10, 4, 10, 0.2);
//    Object3D(convert_to_Mesh(am)).save("../result/VanLoonBench/temp/res.stl");
//    return 0;
    thread t(start_debug_gui, argc, argv);

    int scenario = 0; // 0 - membrane, 1 - shape, 2 - shape + BC
    double Ht = 0.4;
    double mu = 9e2; //kPa
    double theta = M_PI_4, k1 = 1.1e3 /*kPa*/, k2 = 1e3, kappa = 0.15;//.29;
    double P = 90.0_mmHg / 1e3; //kPa
//    ThubricarValvePrms tvp{13.5, 13.5, 15.5, 3.65, 4 / 180.0 * M_PI, 3, Ht};
//    ThubricarValvePrms tvp{12, 12, 17.8, 1, 0 / 180.0 * M_PI, 3, Ht};
    ThubricarValvePrms tvp{12, 12, 17.8, 1, 10 / 180.0 * M_PI, 3, 0.01*Ht};
    double mesh_h = 0.4;//0.8
    double err = 1e-5, diverge_lim = 1e10;
    int maxits = 2e5, pfreq = 5e2;
    double delta = (0.03 + 0.87 * (3 * kappa)) / (mu + 2 * k1);
    bool use_newton = true;

    string main_save = "../result/VanLoonBench/", to_save, extr_dir = "extra/";
    {
        std::ostringstream oss;
        auto tt = std::time(nullptr);
        auto tm = *std::localtime(&tt);
        oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
        to_save = "default/";//oss.str() + "/";//"default/";//TODO: repair this
    }
    to_save = main_save + to_save;
    string to_gen = "../generated/";
    string prefix = "";

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << "\n";
            std::cout << " Command line options: " << "\n";
            std::cout << "  -sc, --scenario   <Choosing of the running scenario: 0 - membrane, 1 - shape, 2 - shape + BC>" << "\n";
            std::cout << "  -ht, --thickness  <Leafs thickness>" << "\n";
            std::cout << "  -p , --pressure   <Diastolic pressure>" << "\n";
            std::cout << "  -mu, --mu         <Shear modulus>" << "\n";
            std::cout << "  -a , --aniso      <Parameters of anisotropy: k1, k2, kappa, theta>" << "\n";
            std::cout << "  -mh, --mesh_h     <Step of generated mesh>" << "\n";
            std::cout << "  -t , --target     <Directory to save results>" << "\n";
            std::cout << "  -tm, --maindir    <Common directory to save results from different simulations>" << "\n";
            std::cout << "  -p , --prefix     <Prefix to file name>" << "\n";
            std::cout << "  -gd, --gendir     <Path to save generated codes>" << "\n";
            std::cout << "  -sp, --slvstop    <Solver stop criterion parameters: err, maxits, pfreq>" << "\n";
            std::cout << "  -vp, --valve      <Valve parameters: R1, R2, H, dh, alpha, nleafs>" << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-sc") == 0 || strcmp(argv[i], "--scenario") == 0) {
            try {
                if (i + 1 < argc)
                    scenario = stoi(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited number of scenario but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-ht") == 0 || strcmp(argv[i], "--thickness") == 0) {
            try {
                if (i + 1 < argc)
                    Ht = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited thickness value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--pressure") == 0) {
            try {
                if (i + 1 < argc)
                    P = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited thickness value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-mu") == 0 || strcmp(argv[i], "--mu") == 0) {
            try {
                if (i + 1 < argc)
                    mu = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited mu value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--aniso") == 0) {
            int j = i;
            const char* v[] = {"k1", "k2", "kappa", "theta"};
            try {
                if (i + 1 < argc) k1 = stof(argv[++i]);
                if (i + 1 < argc) k2 = stof(argv[++i]);
                if (i + 1 < argc) kappa = stof(argv[++i]);
                if (i + 1 < argc) theta = stof(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-vp") == 0 || strcmp(argv[i], "--valve") == 0) {
            int j = i;
            const char* v[] = {"R1", "R2", "H", "dh", "alpha", "nleafs"};
            try {
                if (i + 1 < argc) tvp.R1 = stof(argv[++i]);
                if (i + 1 < argc) tvp.R2 = stof(argv[++i]);
                if (i + 1 < argc) tvp.H = stof(argv[++i]);
                if (i + 1 < argc) tvp.dh = stof(argv[++i]);
                if (i + 1 < argc) tvp.alpha = stof(argv[++i]);
                if (i + 1 < argc) tvp.nleafs = stoi(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-mh") == 0 || strcmp(argv[i], "--mesh_h") == 0) {
            try {
                if (i + 1 < argc)
                    mesh_h = stof(argv[++i]);
            } catch (std::exception& e){
                std::cout << "Waited mesh size value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-sp") == 0 || strcmp(argv[i], "--slvstop") == 0) {
            int j = i;
            const char* v[] = {"err", "maxits", "pfreq"};
            try {
                if (i + 1 < argc) err = stof(argv[++i]);
                if (i + 1 < argc) maxits = stoi(argv[++i]);
                if (i + 1 < argc) pfreq = stoi(argv[++i]);
            } catch (std::invalid_argument& e){
                i--;
            } catch (std::exception& e){
                std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }

        if (strcmp(argv[i], "-gd") == 0 || strcmp(argv[i], "--gendir") == 0) {
            if (i + 1 < argc) to_gen = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0) {
            if (i + 1 < argc) to_save = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-tm") == 0 || strcmp(argv[i], "--maindir") == 0) {
            if (i + 1 < argc) main_save = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) {
            if (i + 1 < argc) to_save = argv[++i];
            continue;
        }
    }

    auto print_input = [&](std::ostream& out) -> std::ostream& {
        return out << "Input parameters: " <<
                         "\n\tscenario = " << scenario <<
                         "\n\tHt = " << Ht <<
                         "\n\tP = "  << P <<
                         "\n\tmu = " << mu <<
                         "\n\tk1 = " << k1 << ", k2 = " << k2 << ", kappa = " << kappa << ", theta = " << theta <<
                         "\n\ttvp = {"
                          " R1 = "    << tvp.R1 <<
                         ", R2 = "    << tvp.R2 <<
                         ", H = "     << tvp.H <<
                         ", dh = "    << tvp.dh <<
                         ", alpha = " << tvp.alpha <<
                         ", nleafs = "  << tvp.nleafs <<
                         "}" <<
                         "\n\tmesh_h = " << mesh_h <<
                         "\n\terr = " << err << ", maxits = " << maxits << ", pfreq = " << pfreq <<
                         "\n\tmain_dir = \"" << main_save << "\", save_dir = \"" << to_save << "\", gen_dir = \"" << to_gen << "\"" <<
                         "\n\tprefix = \"" << prefix << "\"" << std::endl;
    };

    print_input(std::cout);
    {
        auto make_dir = [](std::string save) {
            path path = save;
            auto s = status(path);
            if (!status_known(s) || s.type() != file_type::directory)
                create_directory(path);
        };
        make_dir(main_save);
        make_dir(to_save);
        make_dir(to_save + extr_dir);
        make_dir(to_gen);
    }

    std::ofstream lgfile(to_save + "log.txt");
    print_input(lgfile);

    static const int FIXED = 2, FREE = 1|4|8;
    std::vector<Object3D> _valve = generateThubricarCloseValve(tvp, mesh_h);
    DirichletCond dc = [](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & FIXED) return 0;
        return 7;
    };
    for (auto& l: _valve) l.setDirichletCond(dc);
    for (auto lid = 0; lid < _valve.size(); ++lid) _valve[lid].name = "leaf" + to_string(lid);
    SimplePressureLoad Pf(P);
    HGO_GSTModel Ef(Ht, mu, k1, k2, kappa, theta, to_gen, false);
    Ef.prepareJacobianFunction(false);

    World w;
    vector<ObjectID> vid; vid.reserve(_valve.size());
    std::array<vector<Force_ID>, 3> frid;
    for (auto& leaf: _valve) vid.push_back(w.addObject3D(std::move(leaf), 1));
    double colMargin = Ht / 2;
//    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BulletCollisionManager>();
//    w.setCollider(move(colMan));
//    for (auto id: vid) reinterpret_cast<BulletCollisionManager*>(w.getCollider())->set_margin(id, colMargin);
    for (int i = 0; i < 3; ++i) {
        double phi = M_PI / tvp.nleafs * (2*i - 1);
        SpurnPlane spurnPlane1(P, colMargin, SpurnPlane::Plane_3(CGAL::ORIGIN, Vector{-sin(phi), cos(phi), 0}));
        phi = M_PI / tvp.nleafs * (2*i + 1);
        SpurnPlane spurnPlane2(P, colMargin, SpurnPlane::Plane_3(CGAL::ORIGIN, Vector{sin(phi), -cos(phi), 0}));
        w.addForce(spurnPlane1, vid[i]);
        w.addForce(spurnPlane2, vid[i]);
    }

    if (scenario > 0) {
        for (int i = 0; i < 3; ++i) {
            BendingForce bf(Ef.f, (i==0) && false);
            bf.prepareJacobianFunction((i==0) && false);
            if (scenario == 2) {
                map<E_ind, BendingForce::BoundaryData::Entry> boundary;
                map<E_ind, BendingForce::BoundaryData::Entry> init_boundary;
                Vector shft{cos(2*M_PI*i / tvp.nleafs), sin(2*M_PI*i / tvp.nleafs), 0};
                shft *= tvp.ht / (2 * sin(M_PI / tvp.nleafs));
                {
                    auto &obj = _valve[i];
                    for (auto e: obj.m_mesh.edges()) {
                        auto f = face_around(obj.m_mesh, e);
                        if (f[0].second && f[1].second)
                            continue;
                        auto v = vert_around(obj.m_mesh, e);
                        auto vop = vert_opposite(obj.m_mesh, e);
                        if ((obj.m_boundary[v[0]] & FIXED) && (obj.m_boundary[v[1]] & FIXED) && !(obj.m_boundary[vop.first[0]] & FREE)) {
                            Vector side[2] = {obj.m_x[v[0]] - obj.m_x[v[1]], obj.m_x[v[0]] - obj.m_x[vop.first[0]]};
                            side[0] /= sqrt(side[0].squared_length()); side[1] /= sqrt(side[1].squared_length());
                            Vector l_init = CGAL::cross_product(CGAL::cross_product(side[0], side[1]), side[0]);
                            l_init /= sqrt(l_init.squared_length());
                            if (l_init * side[1] < 0 ) l_init *= -1;
                            l_init *= 0.2; //TODO: remove this
                            Point p[2] = {obj.m_x[v[0]] - shft, obj.m_x[v[1]] - shft};
                            double phi = (atan2(p[0].y(), p[0].x()) + atan2(p[1].y(), p[1].x())) / 2;
                            double z = (p[0].z() + p[1].z()) / 2;
                            double R = tvp.R1 + (tvp.R2 - tvp.R1) * z / tvp.H;
                            Vector ephi(-sin(phi), cos(phi), 0);
                            Vector el = Point(tvp.R1 * cos(phi), tvp.R1 * sin(phi), 0) -
                                        Point(tvp.R2 * cos(phi), tvp.R2 * sin(phi), tvp.H);
                            el /= sqrt(el.squared_length());
                            Vector er = CGAL::cross_product(el, ephi);
                            Vector lt = obj.m_x[v[1]] - obj.m_x[v[0]];
                            lt /= sqrt(lt.squared_length());
                            Vector eb = CGAL::cross_product(er, lt);
                            if (eb.z() > 0) eb *= -1;
                            eb *= 0.2; //TODO: remove this
                            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED,
                                                                            {eb[0], eb[1], eb[2]});
                            init_boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED,
                                                                                 {l_init[0], l_init[1], l_init[2]});
                        } else {
                            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
                            init_boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
                        }
                    }
                }
                bf.set_boundary_data(w.obj(vid[i]), boundary);
            }
            auto bfid = w.addForce(std::move(bf), vid[i]);
            frid[i].push_back(bfid);
        }
    }
    for (int i = 0; i < 3; ++i) {
        auto pfid = w.addForce(Pf, vid[i]);
        auto efid = w.addForce(Ef, vid[i]);
        frid[i].push_back(pfid);
        frid[i].push_back(efid);
    }

    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto time = World3d::Timer();
    int status = 0;
    auto stopCond = World3d::StaticStopCondition(err, diverge_lim, maxits, &pfreq, &time, &lgfile, true);

    if (use_newton) {
        NewtonSolverStepAlgo nssa;
#ifdef USE_INMOST_SOLVERS
        nssa._solver.solver = std::move(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#else
        nssa._solver.solver = std::move(LinearSolver("eigen"));
#endif
        nssa._solver.solver.setVerbosityLevel(0);
        nssa.newton_algo = ChangableNewtonStepAlgo().setLogOut(lgfile).setTauAlgo(
                [](int it, double rel, double abs, NewtonSolverStepAlgo& ns) -> double{
                    //initial relaxation
                    if (it < 25) {
                        if (it == 0) return 1e-4;
                        if (it <  5) return 5e-4;
                        if (it < 10) return 1e-3;
                        if (it < 15) return 1e-2;
                        return 6e-2;
                    }
                    //stability borders
                    if (rel > 0.03) {
                        if (rel > 0.1) return 1e-1;
                        return 0.3;
                    }
                    //newton step for small enough residual
                    return 0.5;
                });

        if (scenario == 0) {
            StaticForceApplierP sfa(delta);
            sfa.addDamper(Damper(0.01 * Ht, 0.8, 100, 1000000));
            w.setForceApplier(sfa);
            maxits = 1000;
            status = w.Simulation(stopCond);
            if (status >= 0) status = 0;
            else std::cout << "Error during relaxation", exit(-1);
        }
        maxits = 2000, pfreq = 1;
        w.setStepSimulationAlgo(std::move(nssa));
    } else {
        StaticForceApplierP sfa(delta);
        sfa.addDamper(Damper(0.01*Ht, 0.8, 100, 1000000));
        w.setForceApplier(sfa);
    }

    time.reset();
    status = w.Simulation(stopCond);

    std::string extr_save = to_save + extr_dir + prefix;
    for (int i = 0; i < vid.size(); ++i) {
        w.obj(vid[i]).save(extr_save + w.obj(vid[i]).name + ".txt");
        ObSaveTag(w.obj(vid[i]), extr_save  + w.obj(vid[i]).name + "_x0.txt", "v:point");
        w.obj(vid[i]).save(to_save + prefix + w.obj(vid[i]).name + ".stl");
    }
    Object3D dup;
    Messurer mes(dup, w.obj(vid[0]), w.obj(vid[1]), w.obj(vid[2]), {0, 0, -1});
    mes.setMargin(2*1.2*colMargin);
    mes.computeCollidingBnd(4);
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
        Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_init.stl", "v:point");
        Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + prefix + mes.getCusps()[i]->name + "_mes_init_col.stl", "v:point",
                            colfilter);
        Messurer::saveObjectShape(mes.m_colission_shapes[i], to_save + prefix + mes.getCusps()[i]->name + "_mes_init_nocol.stl", "v:point",
                            nocolfilter);
        Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_real_nocol.stl", "v:x", nocolfilter);
        Messurer::saveObjectShape(mes.m_colission_shapes[i], extr_save + mes.getCusps()[i]->name + "_mes_real_col.stl", "v:x", colfilter);
    }
    mes.computeHalfCoaptScans();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            Object3D(mes.m_hcScan(i, j).to_Mesh()).save(
                    extr_save + "half_scan_" + std::to_string(i) + "_" + std::to_string(j) + ".stl");
        }
    auto distr = mes.computeHalfCoaptDistrib(500);
    mes.uniteCollidingCuspHalves();
    for (int i = 0; i < 3; ++i) Object3D(mes.m_coaptScan[i].to_Mesh()).save(extr_save + "unite_scan_" + std::to_string(i) + ".stl");
    auto fdistr = mes.computeCoaptDistribution({1000, 1000, 1000});
    Messurer::saveCoaptDistribCSV(fdistr, extr_save + "full_distrib.csv");
    auto Hc = mes.computeHc();
    auto H = distr.getCoaptH();
    auto coaptStat = mes.computeCoaptStatus();
    bool isClosed = ((coaptStat & (1 << 6)) > 0);
    auto bill = mes.computeBillowing();
    auto areas = mes.computeColArea();
    auto print_results = [&](std::ostream& out){
        out << "Hcoapt = " << H << ", Hcentral = " << Hc << "\n"
                                                               "coaptStatus = " << coaptStat << " isClosed = " << isClosed << "\n";
        bill.print(out);
        out << "Coapt areas: \n\t0:" << areas[0] << "\n\t1:" << areas[1] << "\n\t2:" << areas[2] << "\n";
    };
    print_results(std::cout);
    print_results(lgfile);

    Messurer::saveHalfCoaptDistribCSV(distr, extr_save + "_" + "half_distrib.csv");


    std::ofstream comtable;
    auto s = ::status(main_save + "comres.csv");
    if (!exists(s)) {
        comtable.open(main_save + "comres.csv");
        comtable << "res_dir; prefix; Hcpt; Hcent; Acpt;"
                    "scenario; Ht; P; mu; k1; k2; kappa; theta; tvp_R1; tvp_R2; tvp_H; tvp_dh; tvp_alpha; tvp_nleafs; mesh_h; err; maxits; delta; sim_status\n";
    }
    else
        comtable.open(main_save + "comres.csv", std::ios_base::app);

    comtable << "\"" << to_save << "\"" << ";" << "\"" << prefix << "\"" << ";" << H << ";" << Hc << ";" << (areas[0].colArea + areas[1].colArea + areas[2].colArea)/3 << ";"
            << scenario << ";" << Ht << ";" << P << ";" << mu << ";" << k1 << ";" << k2 << ";" << kappa << ";" << theta << ";" << tvp.R1 << ";" << tvp.R2
            << ";" << tvp.H << ";" << tvp.dh << ";" << tvp.alpha << ";" << tvp.nleafs << ";" << mesh_h << ";" << err << ";" << maxits <<";"<< delta
            << ";" << status << std::endl;

    comtable.close();

    t.join();
    return 0;
}
