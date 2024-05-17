//
// Created by alex on 25.10.2021.
//



#include "MeshGenerators/Helper.h"
#include "Object3D.h"
#include "World.h"
#include "Forces/Forces.h"

#include "NSWorldWrapper.h"
#include "Solvers/NonLinearSolverKinsol.h"
#include "Solvers/NonLinearSolverCustom.h"
#include "Benchmarks.h"

using namespace World3d;

//theta - parameter of selection
AniMesh generate_selected_center_half_sphere(double R, double theta, double size){
    AniMesh am;

    const int nVVert = 4,  nLine = 5,  nSurface = 2;
    double VVert[nVVert*3] = {-2*R, -2*R, -2*R, 2*R, 2*R, 2*R, 0.0};
    int LineD[nLine*3] = {2, 1, 1,  1, 4, 1,  2, 3, 1,  3, 2, 1,  4, 1, 1};
    int LineP[nLine*2] = {1, 1,     1, 2,     1, 3,     1, 3,     1, 2}; /*param functions*/
    int exportCurves[] = {0,        2,        4,        4,        2};
    double LineT[] = {theta, M_PI_2,  0, M_PI,  0, M_PI,  M_PI, 2*M_PI,  M_PI, 2*M_PI};
    int SurfL[] = {
            6, 1, 1, 0, 0,
            2, 1, 1, 0, 0
    };
    int SurfI[] = {
            4, 1,  3, 1,  1, 0,  2, 0,  5, 0, 1, 1,
            3, 0,  4, 0
    };
    double SurfT[4*2] = {-1.1, 1.1, -1.1, 1.1,  -1.1, 1.1, -1.1, 1.1}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [rt = sin(theta)](int i, double t, double *pu, double *pv){
        switch(i){
            case 1: { *pu = sin(t)   , *pv = 0.0      ; break; }
            case 2: { *pu = cos(t)   , *pv = sin(t)   ; break; }
            case 3: { *pu = rt*cos(t), *pv = rt*sin(t); break; }
        }
        return 1;
    };
    asd.bounsurf = [R](int i, double u, double v, double *px, double *py, double *pz) {
        double z = 1.0 - u*u - v*v;
        if (fabs(u) < DBL_EPSILON) u = 0;
        if (fabs(v) < DBL_EPSILON) v = 0;
        z = (z < 0.0 || fabs(z) < 10*DBL_EPSILON) ?  0.0 : R * sqrt(z);

        *px = R * u,  *py = R * v, *pz = ((i % 2) ? 1 : -1) * z;
        return 1;
    };
    makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
    am.face_label.resize(am.faces.size() / 3, 0);
    double min_z = R * cos(theta);
    for (int i = 0; i < am.faces.size() / 3; ++i){
        double cz = 0;
        for (int n = 0; n < 3; ++n) cz += am.vertices[3*(am.faces[3*i+n]-1) + 2];
        cz /= 3;
        if (cz > min_z) am.face_label[i] = 4;
        else am.face_label[i] = 2;
    }

    mark_vertices(am);

    return am;
}

int testHalfSphereStatDef(int argc, char* argv[]){
    string res_path = "../result/StatDefSph/HalfSph/";
    string gendir = "../generated";
    bool view = false;
    double P = 2, mu = 1e7;
    double delta = 6e-7, err = 1e-4, abs_err = 1e-10;
    double delta_relax = delta, relax_its = 1000, relax_pfreq = 100;
    int newton_max_its = 100;
    double H0 = 0.001, R = 1, theta = M_PI/3, mesh_h = 0.1;
    int cit = 0;
    std::ofstream logf(res_path + "log.txt");
    std::vector<double> mesh_h_s = {0.1, 0.1/2, 0.1/3.9, 0.1/7.9, 0.1/16};
    std::vector<double> mu_s = {1e4, 1e5, 1e6, 1e7, 1e8};
    for (auto _mesh_h: mesh_h_s){
        mesh_h = _mesh_h;
        auto am = generate_selected_center_half_sphere(R, theta, mesh_h);
        for (auto _mu: mu_s){
            mu = _mu;
            std::cout << "\n\n" << "mu = " << mu << " mesh_h = " << mesh_h << std::endl;
            logf << "\nmu = " << mu << " mesh_h = " << mesh_h << std::endl;
            Object3D _obj{convert_to_Mesh(am, "v:boundary_lbl", "f:face_type")};
            _obj.name = "sph";
            World w;
            auto id = w.addObject3D(std::move(_obj), 1);
            Object3D &obj = w.obj(id);

            auto f_type = obj.m_mesh.property_map<F_ind, int>("f:face_type").first;
            auto Texct = obj.m_mesh.add_property_map<F_ind, std::array<double, 6>>("f:Texct").first;
            auto Teval = obj.m_mesh.add_property_map<F_ind, std::array<double, 6>>("f:Teval").first;
            auto compute_dif = [&Texct, &Teval, &obj, &logf]() {
                auto Ap_ = set_Ap(obj);
                std::array<double, 6> dif = {0}, difC = {0};
                for (auto f: obj.m_mesh.faces())
                    for (int k = 0; k < 6; ++k) {
                        double ldif = abs(Texct[f][k] - Teval[f][k]);
                        dif[k] +=  ldif * Ap_[f];
                        if (ldif > difC[k]) difC[k] = ldif;
                    }
                std::cout << "∫|Texct-Teval|dΩ [6] = { ", logf << "\t∫|Texct-Teval|dΩ [6] = { ";
                for (int i = 0; i < 6; ++i) std::cout << dif[i] << ", ", logf << dif[i] << ", ";
                std::cout << "}" << std::endl, logf << "}" << std::endl;
                std::cout << "||Texct-Teval||_C [6] = { ", logf << "\t||Texct-Teval||_C [6] = { ";
                for (int i = 0; i < 6; ++i) std::cout << difC[i] << ", ", logf << difC[i] << ", ";
                std::cout << "}" << std::endl, logf << "}" << std::endl;
            };
            auto compute_dx = [&Texct, &Teval, &obj, &logf]() {
                std::array<double, 3> dx = {0};
                for (auto v: obj.m_mesh.vertices()) {
                    auto ldx = (obj.m_x[v] - obj.m_x0[v]);
                    for (int k = 0; k < 3; ++k) if (abs(ldx[k]) > dx[k]) dx[k] = abs(ldx[k]);
                }
                std::cout << "||X_inv - X||_inf = {" << dx[0] << ", " << dx[1] << ", " << dx[2] << "}" << std::endl;
                logf << "\t||X_inv - X||_inf = {" << dx[0] << ", " << dx[1] << ", " << dx[2] << "}" << std::endl;
            };
            auto filter_results = [&Texct, &Teval, &f_type, &obj]() {
                for (auto f: obj.m_mesh.faces())
                    if (f_type[f] != 4) std::copy(Texct[f].begin(), Texct[f].end(), Teval[f].begin());
            };
            auto aniso_s = obj.m_mesh.add_property_map<F_ind, std::array<DReal , 3>>("f:aniso_s").first;
            for (auto f: obj.m_mesh.faces()) {
                auto v = vert_around(obj.m_mesh, f);
                auto t = ((obj.m_x0[v[0]] - CGAL::ORIGIN) + (obj.m_x0[v[1]] - CGAL::ORIGIN) +
                          (obj.m_x0[v[2]] - CGAL::ORIGIN)) / 3;
                double phi = atan2(t[1], t[0]);
                double teta = atan2(hypot(t[0], t[1]), t[2]);
                double tt = P * R / 2;
                auto sq = [](auto x) { return x * x; };
                aniso_s[f] = std::array<DReal , 3>{cos(teta)*cos(phi), cos(teta)*sin(phi), -sin(teta)};
                Texct[f][0] = tt * (sq(cos(teta)) * sq(cos(phi)) + sq(sin(phi)));
                Texct[f][1] = tt * (sq(cos(teta)) * sq(sin(phi)) + sq(cos(phi)));
                Texct[f][2] = tt * (sq(sin(teta)));
                //sqrt(2) to achieve ||Teval||_2 = ||T3d||_F
                Texct[f][3] = tt * (-sin(2 * phi) * sq(sin(teta)) / 2) * sqrt(2);
                Texct[f][4] = tt * (-sin(2 * teta) * cos(phi) / 2) * sqrt(2);
                Texct[f][5] = tt * (-sin(2 * teta) * sin(phi) / 2) * sqrt(2);
            }
            obj.set_is_movable_checker([](Object3D &obj, const V_ind &v) -> bool { return !(obj.m_boundary[v] & 2); });
            SimplePressureLoad Pr(P);
            //Force elastic = NeoGookModel(mu, H0, gendir, cit == 0);
            double alpha = 0.05, beta = 1.0;
            Force elastic = MooneyRivlinModel2(H0, mu, alpha, beta, {NAN, NAN, NAN}, gendir, cit == 0, "f:aniso_s");
            elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(cit == 0);
            HEDataGatherer gather(elastic.target<HyperElasticForceBase>()->f, cit == 0);
            w.addForce(elastic, id);
            w.addForce(Pr, id);

            /*init solvers*/
            NSWorldWrapper nsww(w);
            NLProblem nlp = NSWorldWrapper::makeProblem(nsww);
            LinearSolver ls("inner_mptiluc");
            ls.setVerbosityLevel(1);
            ls.SetParameterReal("relative_tolerance", 1e-20);
            double droptol = 15e-4;
            ls.SetParameterReal("drop_tolerance", droptol);
            ls.SetParameterReal("reuse_tolerance", 0.1 * droptol);
            auto monitor_fcn = [&nsww](const double *x) {
                nsww.setX(x);
                nsww.RenderIteration();
                return 0;
            };

            //initialize KINSOL solver
            NonLinearSolverKinsol nlsp(nlp);
            nlsp.SetLinearSolver(&ls);
            static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
            nlsp.SetVerbosityLevel(1);
            nlsp.SetInitialGuess(nsww.getCurrentX());
            nlsp.SetFuncNormTol(abs_err);
            nlsp.SetMaxNewtonStep(10 * H0 * sqrt(obj.m_mesh.num_vertices()) * 1e10);
            nlsp.SetParameterInt("MaxSetupCalls", 1);
            nlsp.SetParameterInt("MaxSubSetupCalls", 1);
            nlsp.SetScaledStepTol(1e-10);
            nlsp.SetMaxBetaFails(40);
            nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
            World3d::Timer sym_time, com_time;
            com_time.reset();
            nlsp.SetInfoHandlerFn([&time = sym_time, &com_time](const char *module, const char *function, char *msg) {
                std::cout << "[" << module << "] " << function << "\n   " << msg << "\n" << "time = " << time.elapsed()
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

            //initialize relaxation solver
            auto relaxate_solver = [Ht = H0, &w, &sym_time, oid = id](double &delta, int nits,
                                                                      int pfreq) {
                if (nits <= 0) return;
                std::cout << "delta = " << delta << "\n";
                static double diverge_lim = 1e15;
                static double err = 1e-4;
                auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &nits, &pfreq, &sym_time,
                                                             nullptr,
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
            auto print_iterate_info = [&fe = func_evals, &je = jac_evals, &com_time](int it) {
                std::cout << "ITERATE INFO: it: " << it << " com_time = " << com_time.elapsed() <<
                          " accum_func_evals = " << fe << " accum_jac_evals = " << je << std::endl;
            };

            std::thread *view_window;
            if (view) {
                view_window = new std::thread(start_debug_gui, argc, argv);
                w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
            }

            bool compute_custom = true;
            if (compute_custom) {
                relaxate_solver(delta_relax, relax_its, relax_pfreq);
                nlsp.SetInitialGuess(nsww.getCurrentX());
                nlsp.SetNumMaxIters(newton_max_its);
                bool slvFlag = nlsp.Solve();
                std::cout << result_str_nlsp(nlsp) << std::endl;
                func_evals += nlsp.GetNumFuncEvals(), jac_evals += nlsp.GetNumJacEvals();
                print_iterate_info(0);

                gather.gatherData(obj, [&Teval](const HEDataGatherer::CompDat &dat) {
                    std::copy(dat.T, dat.T + 6, Teval[dat.f].begin());
                    for (int k = 3; k < 6; ++k) Teval[dat.f][k] *= sqrt(2); //to achieve ||Teval||_2 = ||T3d||_F
                });
            } else {
                StaticDefiniteGather::Traits traits;
                StaticDefiniteGather statgath;
                traits.setH(H0).setMu(mu).setDelta(delta).setErr(abs_err);
                traits.setSolveMethod(StaticDefiniteGather::Traits::SIMPLE_NEWTON).setMaxIters(15e5).setFreq(1);
#ifdef USE_INMOST_SOLVERS
                traits.setLinearSolver(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#else
                traits.setLinearSolver(LinearSolver("eigen"));
#endif
                statgath.setTraits(traits).setObject(obj);
                statgath.m_t.setP(P);
                statgath.computeTensionedState();
                statgath.gatherMyData(
                        [&Teval](const StaticDefiniteGather::ObjsData &dat) -> void {
                            std::pair<int, int> q3[] = {{0, 0},
                                                        {1, 1},
                                                        {2, 2},
                                                        {0, 1},
                                                        {0, 2},
                                                        {1, 2}},
                                    q2[] = {{0, 0},
                                            {1, 1},
                                            {1, 0}};
                            for (int i = 0; i < 6; ++i) Teval[dat.f][i] = dat.T(q3[i].first, q3[i].second);
                            for (int k = 3; k < 6; ++k) Teval[dat.f][k] *= sqrt(2); //to achieve ||Teval||_2 = ||T3d||_F
                        });
            }

            compute_dx();
            compute_dif();
            filter_results();
            compute_dif();
            auto Tdif = obj.m_mesh.add_property_map<F_ind, std::array<double, 6>>("f:Texct-Teval").first;
            for (auto f: obj.m_mesh.faces()) for (int k = 0; k < 6; ++k) Tdif[f][k] = Texct[f][k] - Teval[f][k];
            stringstream ss;
            ss << "h=" << std::setprecision(5) << mesh_h << std::setprecision(1) << "_mu=" << std::scientific << mu;
            string postfix = "_" + ss.str();
            logf << "Save result in " << res_path + "tension" + postfix + ".vtk" << std::endl;
            obj.save(res_path + "tension" + postfix + ".vtk");

            if (view) {
                view_window->join();
                delete view_window;
            }

            cit++;
        }
    }

    return 0;
}
