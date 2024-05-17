#include "CircleTest.h"
#include <boost/program_options.hpp>

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

namespace po = boost::program_options;

int main(int argc, char* argv[]){
    double R = 1, H = 0.07, p_star = 0.0, mu = 0.35;
    double dp_star = 0.1, ps_max = 1.2;
    double h_step = 0.05;
    int eforce_type = 0;
    static const int NHK_TP = 0, ZHOU_TP = 1, OGDEN_BR_TP = 2;
    std::string gendir = "../../../generated"; bool regenerate = true;
    std::string result_dir = "./";
    std::string mesh_fname = "res.vtk", data_fname = "distribution.csv";
    bool with_gui_window = true, with_prerelaxation = true;
    bool extend_data = false;

    try {
        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
                ("version,v", "print version string")
                ("help,h", "produce help message")
                ;
        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
                ("material_model,m", po::value<int>(&eforce_type)->default_value(0), "Choosing of the material model: 0 - NHK, 1 - ZHOU, 2 - OGDEN_BRAIN")
                ("dps", po::value<double>(&dp_star)->default_value(0.1), "Step of dimensionless p changing")
                ("max_ps", po::value<double>(&ps_max)->default_value(1.2), "Max dimensionless p value")
                ("mesh_step", po::value<double>(&h_step)->default_value(0.05), "Dimensionless step of generated mesh, h = R*h_dimensionless")
                ("view", "Create ImGui debug view window")
                ("with_prerelax", "Make initial prerelaxation")
                ("regen_expressions", "Regenerate expressions")
                ("extend_data", "Extend database by new data")
                ("root_directory", po::value<std::string>(&result_dir)->default_value("./"), "Root directory for saving and reading of results")
                ("mesh_name", po::value<std::string>(&mesh_fname)->default_value("res.vtk"), "Name of file to save mesh in root directory")
                ("date_name", po::value<std::string>(&data_fname)->default_value("distribution.csv"), "Name of file to save gathered data in root directory")
                ("radius,R", po::value<double>(&R)->default_value(1), "Radius R of circle, mm")
                ("thickness,t", po::value<double>(&H)->default_value(0.07), "Leafs thickness, mm")
                ("shear_modulus", po::value<double>(&mu)->default_value(0.35), "Shear modulus mu, kPa")
                ("gen_directory", po::value<std::string>(&gendir)->default_value("../../../generated"), "Path to work with generated expressions")
                ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config);

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << cmdline_options << "\n";
            return 0;
        }
        if (vm.count("version")){
            cout << "ConjugatePairs version 0.1.0\n";
            return 0;
        }
        with_gui_window = vm.count("view");
        with_prerelaxation = vm.count("with_prerelax");
        regenerate = vm.count("regen_expressions");
        extend_data = vm.count("extend_data");
        {
            path p(result_dir);
            if (!exists(status(p)) || !is_directory(status(p)))
                throw std::runtime_error("root_directory = \"" + result_dir + "\" is not exists");
        }
        if (result_dir.back() != '/') result_dir += "/";
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return -1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        return -2;
    }

    p_star = dp_star < ps_max ? dp_star : ps_max;
    double p_val = mu * H / R * p_star;
    auto _obj = Object3D(convert_to_Mesh(generate_circle(R, R*h_step), "v:boundary_lbl")).setName("circle");
    _obj.set_is_movable_checker([](Object3D& obj, const V_ind& v) { return !(obj.m_boundary[v] % 2); });
    Force Fe;
    switch (eforce_type){
        case NHK_TP: {
            Fe = NeoGookModel(mu, H, gendir, regenerate);//NeoGookModelOpt(mu, H);
            Fe.target<HyperElasticForceBase>()->prepareJacobianFunction(regenerate);
            break;
        }
        case ZHOU_TP:{
            Fe = ZhouModel(mu, H, gendir, regenerate);
            Fe.target<HyperElasticForceBase>()->prepareJacobianFunction(regenerate);
            break;
        }
        case OGDEN_BR_TP:{
            double alpha = -25.3;
            Fe = OgdenBPModel(mu, H, alpha, gendir, regenerate);
            Fe.target<HyperElasticForceBase>()->prepareJacobianFunction(regenerate);
            break;
        }
        default:
            throw std::runtime_error("Faced unknown index of force type = " + std::to_string(eforce_type));
    }
    Force Pr = SimplePressureLoad([&p_val](Object3D&, F_ind)->DReal{ return p_val; });

    World w;
    std::shared_ptr<std::thread> front_end;
    if (with_gui_window) {
        front_end = make_shared<std::thread>(start_debug_gui, argc, argv);
        auto renderer = std::make_unique<World3d::DefaultRenderer>();
        renderer->set_interconnector(&g_interconnection);
        w.setRenderer(std::move(renderer));
    }
    auto id = w.addObject3D(move(_obj), 1);
    auto pf_id = w.addForce(Pr, id);
    auto ef_id = w.addForce(Fe, id);
    // Force bending = BendingForce(NeoGookModel(mu, H, gendir, false).f, false);
    // bending.target<BendingForcetrue>()->prepareJacobianFunction(false);
    // w.addForce(bending, id);
    World3d::Timer sym_time;
    w.UpdateForces();
    double abs_err_init = w.getWorldNorm() / p_star;

    if (with_prerelaxation) {
        double delta = 1.0;
        int maxits = std::min(int(1e3*0.05/(delta*h_step)), 10000);
        double err = 1e-4;
        double diverge_lim = 1e15;
        int pfreq = 100;
        auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &maxits, &pfreq, &sym_time, nullptr, true);
        StaticForceApplierP sfa(delta);
        w.setForceApplier(sfa);

        w.Simulation(stopCond);
    }

    NSWorldWrapper nsww(w);
    NonLinearSolverKinsol nlsp(NSWorldWrapper::makeProblem(nsww));
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
    nlsp.SetNumMaxIters(500);
    nlsp.SetMaxNewtonStep(R / 10 * sqrt(nlsp.m_prob.m_dofs));
    bool slvFlag = false;
    nlsp.SetParameterInt("MaxSetupCalls", 1);
    nlsp.SetParameterInt("MaxSubSetupCalls", 1);
    nlsp.SetScaledStepTol(1e-7);
    nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
    //nlsp.SetMaxBetaFails(80);
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
    auto watch_num = w.obj(id).m_forces[ef_id].target<HyperElasticForceBase>()->f.add_watch(DefaultHyperElasticForce::compute_HS_tensor_expr);
    std::ofstream out;
    if (exists(status(path(result_dir + data_fname))) && extend_data){
        out.open(result_dir + data_fname, std::ios::app);
    } else {
        out.open(result_dir + data_fname, std::ios::trunc);
        out << "ps, R, delta, epsilon, gamma, pi, sigma, tau_s\n";
    }
    auto conjugate_pairs_collect = [&](){
        auto& o = w.obj(id);
        auto* force = w.obj(id).m_forces[ef_id].target<HyperElasticForceBase>();
        auto L_ = set_L(o);
        auto D_ = set_D_vecs(o);
        auto x_ = o.m_x;
        auto p_ = o.m_x0;
        auto Ap_ = set_Ap(o);
        auto S_ = set_S(o);
        auto S2d_ = set_S2d(o, check_flat_initial_template(o));
            
        for (auto f: o.m_mesh.faces()){
            auto v = vert_around(o.m_mesh, f);
            auto Stensor = force->f.eval_watch(o, f, watch_num);
            auto Ap = Ap_[f], Aq = sqrt(S_[f].squared_length());
            Eigen::Matrix<DReal, 3, 3> Q, P, D;
       
            Q << x_[v[0]][0], x_[v[1]][0], x_[v[2]][0],
                x_[v[0]][1], x_[v[1]][1], x_[v[2]][1],
                x_[v[0]][2], x_[v[1]][2], x_[v[2]][2];
            P << p_[v[0]][0], p_[v[1]][0], p_[v[2]][0],
                p_[v[0]][1], p_[v[1]][1], p_[v[2]][1],
                p_[v[0]][2], p_[v[1]][2], p_[v[2]][2];
            D << D_[f][0][0], D_[f][1][0], D_[f][2][0],
                D_[f][0][1], D_[f][1][1], D_[f][2][1],
                D_[f][0][2], D_[f][1][2], D_[f][2][2];
            double radius = (P * Eigen::Matrix<DReal, 3, 1>{1.0/3, 1.0/3, 1.0/3}).norm();
            Eigen::Matrix<DReal, 2, 3> L = Eigen::Map<Eigen::Matrix<DReal, 2, 3>>(L_[f].data()); //L = S2d * D
            Eigen::Matrix<DReal, 3, 2> F2d = Q * L.transpose();
            Eigen::Matrix<DReal, 2, 2> C2d = F2d.transpose() * F2d;
            Eigen::Matrix<DReal, 3, 2> R;
            R.col(0) = F2d.col(0) / F2d.col(0).norm();
            R.col(1) = F2d.col(1) - F2d.col(1).dot(R.col(0))*R.col(0); R.col(1) /= R.col(1).norm();
            Eigen::Matrix<DReal, 2, 2> Ftld, Ftld1 = R.transpose() * F2d; 
            Ftld << sqrt(C2d(0, 0)), C2d(0, 1) / sqrt(C2d(0, 0)),
                    0,  sqrt(C2d(1, 1) - C2d(0, 1) * C2d(0, 1) / C2d(0, 0));
            if ((Ftld - Ftld1).norm() > Ftld.norm() * 1e-7){
                std::cout << f << "\n";
                std::cout << "Ftld = \n" << Ftld << "Ftld1 = \n" << Ftld1 << std::endl;
            }       
            std::array<DReal, 3> eta = 
                {0.5 * log(Ftld(0, 0) * Ftld(1,1)),
                 0.5 * log(Ftld(0, 0) / Ftld(1,1)),
                 Ftld(0, 1) / Ftld(0,0) };    
            Eigen::Matrix<DReal, 2, 2> Ttilde = Ftld1 * Stensor * Ftld1.transpose(); 
            std::array<DReal, 3> tau = 
                {Ttilde(0, 0) + Ttilde(1, 1),
                 Ttilde(0, 0) - Ttilde(1, 1),
                 Ftld(0, 0) / Ftld(1,1) * Ttilde(0, 1) };
            out << std::scientific << std::setprecision(15) 
                << p_star << ", " << radius 
                << ", " << eta[0] << ", " << eta[1] << ", " << eta[2] 
                << ", " << tau[0] << ", " << tau[1] << ", " << tau[2] << "\n";    
        }
    };
    
    p_star = 0;
    while (abs(p_star - ps_max) > std::min(dp_star, ps_max)*1e-10){
        p_star += std::min(ps_max - p_star, dp_star);
        p_val = mu * H / R * p_star;
        slvFlag = nlsp.Solve(); std::cout << "ps = " << p_star << "\n" << result_str_nlsp(nlsp) << std::endl;
        if (!slvFlag) break;
        conjugate_pairs_collect();
    }
    out.close();

    w.obj(id).save(result_dir + mesh_fname);

    if (with_gui_window){
        front_end->join();
        front_end.reset();
    }

    return 0;
}