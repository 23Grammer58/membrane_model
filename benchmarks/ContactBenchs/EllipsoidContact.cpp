//
// Created by alex on 19.11.2021.
//
#include "../BenchCommon.h"

#include "AVSim/Core/MeshGen/Helper.h"
#include "AVSim/Core/Object3D.h"
#include "AVSim/Core/World.h"
#include "AVSim/Core/Forces/Forces.h"
#include "AVSim/Core/Collision/CGALCollissionManager.h"
#include "AVSim/Core/MeshGen/StretchFigure.h"
#include "AVSim/Core/Collision/BulletCollisionManager.h"

using namespace std;
using namespace World3d;

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

struct ProblemDiscr{
    enum TestType{
        Ellipsoid = 0,
        ThreeLeaf,
        Overlap,
        END
    };
    string main_path = "../result/ContactTests/";
    string gendir = "../generated";
    TestType tt = Ellipsoid;
    double H0 = 0.02, mesh_h = 0.1;
    std::vector<double> geometry = {1, 0.3, 0.5};
    double P = 1, mu = 50;
    double delta = 0.01, err = 1e-7, diverge_lim = 1e15;
    int pfreq = 100, maxits = 1e5;
    double margin = 0.01;
    bool is_penalty = false;
    bool with_view = false;

    bool only_save_input = false;

    std::vector<std::pair<string, string>> serialize(){
        std::vector<std::pair<string, string>> res;
        switch (tt) {
            case Ellipsoid: res.emplace_back("TestType", "Ellipsoid"); break;
            case ThreeLeaf: res.emplace_back("TestType","ThreeLeaf"); break;
            case Overlap: res.emplace_back("TestType", "Overlap") ; break;
            default: res.emplace_back("TestType", "NoType");
        }
        std::ostringstream oss;
        oss << "{ ";
        if (!geometry.empty()) oss << geometry[0];
        for (int i = 1; i < geometry.size(); ++i) oss << ", " << geometry[i];
        oss << " }";
        res.emplace_back("geometry", oss.str());
#define ADD(X) res.emplace_back(#X,to_string(X))
        ADD(H0);
        ADD(mesh_h);
        ADD(P);
        ADD(mu);
        ADD(margin);
        ADD(is_penalty);
        ADD(delta);
        ADD(maxits);
        ADD(pfreq);
        ADD(err);
        ADD(diverge_lim);
        ADD(with_view);
        res.emplace_back("main_path", main_path);
        res.emplace_back("gendir", gendir);
#undef ADD
        return res;
    }

    static ProblemDiscr Init(int argc, char* argv[]){
        ProblemDiscr p;
        bool m_set_delta = false;
        bool m_set_margin = false;

        auto print_help_message = [](std::string prefix = "", bool is_exit = false){
            std::cout <<prefix<< "Help message: " << "\n";
            std::cout <<prefix<< " Command line options: -cs #cs [other params]" << "\n";
            std::cout <<prefix<< "  -cs, --case         <Choose case and set default parameters for it: 0 - Ellipsoid, 1 - ThreeLeaf, 2 - Overlap>" << "\n";
            std::cout <<prefix<< "  -ht, --thickness    <Thickness>" << "\n";
            std::cout <<prefix<< "  -pr, --pressure     <Pressure>" << "\n";
            std::cout <<prefix<< "  -mu, --elastic_mod  <Shear modulus>" << "\n";
            std::cout <<prefix<< "  -mh, --mesh_size    <Step of generated mesh>" << "\n";
            std::cout <<prefix<< "  -dt, --delta        <Relaxation solver parameter>" << "\n";
            std::cout <<prefix<< "  -gm, --geometry     <Geometrical parameters for simulations:\n"
                      <<prefix<< "                           Ellipsoid: R, b, d\n"
                      <<prefix<< "                           ThreeLeaf: R, d\n"
                      <<prefix<< "                           Overlap: R, d_right, d_left, dist>\n";
            std::cout <<prefix<< "  -sp, --slvstop      <Solver stop criterion parameters: maxnits, f_abs_err, diverge_lim>" << "\n";
            std::cout <<prefix<< "  -mg, --margin       <Collision shape thickness>" << "\n";
            std::cout <<prefix<< "  -vw, --view         <Create ImGui debug view window>" << "\n";
            std::cout <<prefix<< "  -fq, --print_freq   <Frequency of printing info during simulation (one print per print_freq iterations)>" << "\n";
            std::cout <<prefix<< "  -cl, --col_type     <Collider type: 0 - constraint, 1 - penalty>" << "\n";
            std::cout <<prefix<< "  -md, --maindir      <Common directory to save results>" << "\n";
            std::cout <<prefix<< "  -gd, --gendir       <Path to save generated codes>" << "\n";
            std::cout <<prefix<< "  -si, --save_input   <Only save input and cancel>" << "\n";
            std::cout <<prefix<< "  -h,  --help         <Print this message and cancel>" << "\n";
            std::cout << std::endl;
            if (is_exit) exit(0);
        };
        if (argc == 1) print_help_message("", true);
        for (int i = 1; i < argc; i++) {
            //Print help message and exit
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                print_help_message("", true);
            }
            if (strcmp(argv[i], "-cs") == 0 || strcmp(argv[i], "--case") == 0){
                try {
                    int var = -1;
                    if (i + 1 < argc) {
                        var = stoi(argv[++i]);
                    } else throw std::runtime_error("Not found case number");
                    switch (var) {
                        case Ellipsoid: {
                            p.tt = Ellipsoid;
                            p.H0 = 0.02, p.mesh_h = 0.1;
                            p.geometry = {1, 0.3, 0.5};
                            p.P = 1, p.mu = 50;
                            p.delta = 0.01, p.err = 1e-7, p.diverge_lim = 1e15;
                            p.pfreq = 100, p.maxits = 1e5;
                            p.margin = p.H0 / 2;
                            break;
                        }
                        case ThreeLeaf: {
                            p.tt = ThreeLeaf;
                            p.H0 = 0.05, p.mesh_h = 0.6;
                            p.geometry = {10, 3};
                            p.P = 80 * 133.322, p.mu = 5e6;
                            p.delta = 0.025/(p.mu*p.H0), p.err = 1e-7, p.diverge_lim = 1e15;
                            p.pfreq = 100, p.maxits = 1e5;
                            p.margin = p.H0/2;
                            break;
                        }
                        case Overlap: {
                            p.tt = Overlap;
                            p.H0 = 0.2, p.mesh_h = 0.6;
                            p.geometry = {10, 5, 0, 1};
                            p.P = 80 * 133.322, p.mu = 5e6 * 0.05 / p.H0;
                            p.delta = 0.025/(p.mu*p.H0), p.err = 1e-7, p.diverge_lim = 1e15;
                            p.pfreq = 100, p.maxits = 2e5;
                            p.margin = p.H0/2;

                            break;
                        }
                        default: throw std::runtime_error("Wrong case: Should be 0 - Ellipsoid, 1 - ThreeLeaf, 2 - Overlap, but get " + to_string(var));
                    }
                } catch (std::exception& e){
                    std::cerr << "Waited " << "case" << " but error happens: '" << e.what() << "'\n"; --i;
                    print_help_message("", true);
                }
                continue;
            } else if (i == 1){
                std::cerr << "Wrong usage: first argument should be \"-cs\" #cs, but found \"" << argv[i] << "\"\n";
                print_help_message("", true);
            }
#define READVAR(SHRT, LNG, var, procvar, varname) \
            if (strcmp(argv[i], SHRT) == 0 || strcmp(argv[i], LNG) == 0) {\
                try {\
                    if (i + 1 < argc) {\
                        var = procvar(argv[++i]);\
                    }\
                } catch (std::exception& e){\
                    std::cout << "Waited " << varname << " but error happens: '" << e.what() << "'\n"; --i;\
                }\
                continue;\
            }
#define READVARANDSETFLAG(SHRT, LNG, var, procvar, varname, FLAG) \
            if (strcmp(argv[i], SHRT) == 0 || strcmp(argv[i], LNG) == 0) {\
                try {\
                    if (i + 1 < argc) {\
                        var = procvar(argv[++i]);\
                        FLAG = true;\
                    }\
                } catch (std::exception& e){\
                    std::cout << "Waited " << varname << " but error happens: '" << e.what() << "'\n"; --i;\
                }\
                continue;\
            }
#define READSTR(SHRT, LNG, var) \
            if (strcmp(argv[i], SHRT) == 0 || strcmp(argv[i], LNG) == 0) {\
                if (i + 1 < argc) var = argv[++i];\
                continue;\
            }
#define READVEC(SHRT, LNG, VEC, PROCELEM, VECNAME) \
            if (strcmp(argv[i], SHRT) == 0 || strcmp(argv[i], LNG) == 0) { \
                VEC.clear();                       \
                for (int j = 0; j < 10 && (i+1 < argc && (argv[i+1][0] != '-' || std::isdigit(argv[i+1][1]))); ++j) {\
                    try {\
                        VEC.push_back(PROCELEM(argv[++i]));\
                    } catch (std::exception& e) {\
                        std::cout << "Waited " << VECNAME << "[" << j << "] but error happens: '" << e.what() << "'\n"; --i;\
                    }\
                }\
                continue;\
            }

            READVAR("-ht", "--thickness", p.H0, stod, "H0");
            READVAR("-pr", "--pressure", p.P, stod, "P");
            READVAR("-mu", "--elastic_mod", p.mu, stod, "mu");
            READVAR("-mh", "--mesh_size", p.mesh_h, stod, "mesh_h");
            READVARANDSETFLAG("-dt", "--delta", p.delta, stod, "delta", m_set_delta);
            READVEC("-gm", "--geometry", p.geometry, stod, "geometry");
            if (strcmp(argv[i], "-sp") == 0 || strcmp(argv[i], "--slvstop") == 0) {
                int j = i;
                const char* v[] = {"maxints", "f_ab_err", "diverge_lim"};
                try {
                    if (i + 1 < argc) p.maxits = stoi(argv[++i]);
                    if (i + 1 < argc) p.err = stod(argv[++i]);
                    if (i + 1 < argc) p.diverge_lim = stod(argv[++i]);
                } catch (std::invalid_argument& e){
                    i--;
                } catch (std::exception& e){
                    std::cout << "Waited " << v[i - j] << " value but error happens: '" << e.what() << "'\n"; --i;
                }
                continue;
            }
            READVARANDSETFLAG("-mg", "--margin", p.margin, stod, "margin", m_set_margin);
            auto get_bool = [](std::string val){ int b = stoi(val); return b != 0; };
            READVAR("-vw", "--view", p.with_view, get_bool, "view");
            READVAR("-fq", "--print_freq", p.pfreq, stoi, "pfreq");
            if (strcmp(argv[i], "-cl") == 0 || strcmp(argv[i], "--col_type") == 0) {
                try {
                    int var = -1;
                    if (i + 1 < argc) {
                        var = stoi(argv[++i]);
                        switch (var) {
                            case 0: p.is_penalty = false; break;
                            case 1: p.is_penalty = true; break;
                            default:
                                std::cerr << "Unknown collision type: \"" << var << "\". Should be 0 - constraint, 1 - penalty.\n";
                                print_help_message("", true);
                        }
                    }
                } catch (std::exception& e){
                    std::cout << "Waited " << "collision_type" << " but error happens: '" << e.what() << "'\n"; --i;
                }
                continue;
            }
            READSTR("-md", "--maindir", p.main_path);
            READSTR("-gd", "--gendir", p.gendir);
#undef READVEC
#undef READSTR
#undef READVARANDSETFLAG
#undef READVAR
            if (strcmp(argv[i], "-si") == 0 || strcmp(argv[i], "--save_input") == 0){
                p.only_save_input = true;
                continue;
            }
            if (true){
                std::cerr << "Faced unparsed value: \"" + std::string(argv[i]) + "\"" << std::endl;
                print_help_message("", true);
            }
        }
        if (argc > 1) {
            if (!m_set_margin) p.margin = p.H0 / 2;
            if (!m_set_delta) {
                switch (p.tt) {
                    case Ellipsoid: p.delta = 0.01 / (p.mu * p.H0); break;
                    case ThreeLeaf: p.delta = 0.025 / (p.mu * p.H0); break;
                    case Overlap: p.delta = 0.025 / (p.mu * p.H0);  break;
                }
            }
        }
        auto ps = p.serialize();
        std::cout << "ProblemDiscr{\n";
        for (auto& i: ps) std::cout << "\t" << i.first << " = " << i.second << "\n";
        std::cout << "}\n" << std::endl;

        return p;
    }
};

int ellipsoidContact(int argc, char* argv[], ProblemDiscr p){
    string res_path = p.main_path + "EllipsoidTest/";
    string gendir = p.gendir;
    bool view = p.with_view;

    if (p.geometry.size() < 3) throw std::runtime_error("Geometry should consist 3 vals");
    double H0 = p.H0, R = p.geometry[0], b = p.geometry[1], d = p.geometry[2], mesh_h = p.mesh_h;
    double P = p.P, mu = p.mu;//50+
    double delta = p.delta, err = p.err, diverge_lim = p.diverge_lim;
    int pfreq = p.pfreq, maxits = p.maxits;
    double bullet_margin = p.margin;
    bool is_penalty = p.is_penalty;
    if (is_penalty) bullet_margin *= 2;

    if (!exists(res_path)) create_directory(res_path);
    auto am = generate_EllipsoidPart(R, b, d, mesh_h);
    Object3D _obj0{convert_to_Mesh(am, "v:boundary_lbl")};
    _obj0.name = "obj0";
    for (int v = 0; v < am.vertices.size()/3; ++v) am.vertices[3*v + 2] = -am.vertices[3*v + 2] + 3*b;
    invertFaceOrientation(am);
    Object3D _obj1{convert_to_Mesh(am, "v:boundary_lbl")};
    _obj1.name = "obj1";
//    _obj0.save(res_path + "obj0.vtk");
//    _obj1.save(res_path + "obj1.vtk");
    auto is_movable = [](Object3D &obj, const V_ind &v) -> bool { return !(obj.m_boundary[v] & 2); };
    World w;
    std::array<ObjectID, 2> id = {w.addObject3D(std::move(_obj0), 1), w.addObject3D(std::move(_obj1), 1)};
    for (auto i: id) w.obj(i).set_is_movable_checker(is_movable);
    SimplePressureLoad Pr(P);
    Force elastic = NeoGookModelOpt(mu, H0);
    w.addForce(elastic);
    w.addForce(Pr);
    struct PerformTimers{
        double find_col_tm = 0, sol_col_tm = 0, updt_f_tm = 0, appl_f_tm = 0, appl_nxt_tm = 0;
        double ffct = 0, fsct = 0, fuft = 0, faft = 0, fant = 0;
        World3d::Timer m_timer;
    };
    PerformTimers perform;

    auto stepAlgo = [&p = perform](World* w)->int{
        int status = 0;
        p.m_timer.reset();
        if ((status = w->UpdateForces()) < 0) return status;
        p.updt_f_tm += p.m_timer.elapsed(); p.m_timer.reset();
        if ((status = w->ApplyForces()) < 0) return status;
        p.appl_f_tm += p.m_timer.elapsed(); p.m_timer.reset();

        w->getCollider()->findCollisions();
        p.find_col_tm += p.m_timer.elapsed(); p.m_timer.reset();
        w->getCollider()->solveCollisions();
        p.sol_col_tm += p.m_timer.elapsed(); p.m_timer.reset();

        if ((status = w->ApplyNext()) < 0) return status;
        p.appl_nxt_tm += p.m_timer.elapsed(); p.m_timer.reset();
        return status;
    };

    string prefix = "bul_constr_", postfix;
    if (is_penalty) prefix = "bul_penalty_";
    {
        std::ostringstream oss;
        auto tt = std::time(nullptr);
        auto tm = *std::localtime(&tt);
        oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
        postfix = "_" + oss.str();
    }
    if (p.only_save_input){
        for (auto i: id) w.obj(i).save(res_path + w.obj(i).name + "_R=" + to_string(R) + "_b=" +to_string(b) + "_d=" + to_string(d) + ".vtk");
        return 0;
    }
    std::ofstream shifts(res_path + prefix + "shift" + postfix + ".csv");
    shifts << "it, delta_x, wall_sol - wall, find_col_time, sol_col_time, update_force_time, apply_force_time, apply_next_time\n";
    string stat_name = res_path + "stat" + ".cvs";
    bool is_new = exists(stat_name);
    {
        std::ofstream stat(stat_name, std::ios::app);
        if (!is_new)
            stat << "prefix, postfix, H0, R, b, d, mh, P, mu, delta, maxit, param"
                 << ", find_col_time, sol_col_time, update_force_time, apply_force_time, apply_next_time, full_time\n";
    }
    std::ostringstream stat;
    stat << prefix << ", " << postfix << ", "
         << H0 << "," << R << ", " << b << ", " << d << ", " << mesh_h << ", "
         << P << ", " << mu << ", " << delta << ", " << maxits << ", {margin=" << bullet_margin << "}, ";
//    unique_ptr<CollisionManagerBase> colMan = std::make_unique<CGALColManVFfull>();
//    w.setCollider(move(colMan));
//    for (auto i: id) reinterpret_cast<CGALColManVFfull*>(w.getCollider())->set_margin(i, bullet_margin);

    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BulletCollisionManager>();
    if (is_penalty)
        reinterpret_cast<BulletCollisionManager*>(colMan.get())->setCustomSolveColAlgo(PenaltyCollisionAlgo());
    w.setCollider(move(colMan));
    for (auto i: id) reinterpret_cast<BulletCollisionManager*>(w.getCollider())->set_margin(i, bullet_margin);

    std::thread *view_window;
    if (view) {
        view_window = new std::thread(start_debug_gui, argc, argv);
        w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    }

    StaticForceApplierP sfa(delta);
    //sfa.addDamper(Damper(0.01*H, 0.8, 100, 100000).setReInitDxParams(1e7, 0.4));
    w.setForceApplier(sfa);

    World3d::Timer sym_time;
    auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &maxits, &pfreq, &sym_time, nullptr, true);
    stopCond.addExternalCall([freq = pfreq, &shifts, &obj0 = w.obj(id[0]), wal_pos = 1.5*b-0.5*H0, &p = perform](StepSimInfo& info, World *w) -> bool {
        p.ffct += p.find_col_tm; p.fsct += p.sol_col_tm; p.fuft += p.updt_f_tm; p.faft += p.appl_f_tm; p.fant += p.appl_nxt_tm;
        static double resid_init = 1.0;
        double resid = w->getWorldShift();
        V_ind max_v = *std::max_element(obj0.m_mesh.vertices().begin(), obj0.m_mesh.vertices().end(), [&obj0](V_ind v0, V_ind v1) { return obj0.m_x[v0][2] < obj0.m_x[v1][2];});
        double wall_place = wal_pos  - obj0.m_x[max_v][2];

        if (info.it > 0) shifts << info.it << ", " << resid << ", " << wall_place
            << ", " << p.find_col_tm << ", " << p.sol_col_tm << ", " << p.updt_f_tm << ", " << p.appl_f_tm << ", " << p.appl_nxt_tm << "\n";
        p.find_col_tm = p.sol_col_tm = p.updt_f_tm = p.appl_f_tm = p.appl_nxt_tm = 0;
        if (info.it == 1) resid_init = w->getWorldShift();
        if (info.it > 0 && info.it % freq  == 0) {
            double eps = resid / resid_init;
            std::cout << "Shift[" << info.it << "] = " << resid << " (" << eps << ")" << "\n";
        }
        return true;
    });
    w.setStepSimulationAlgo(stepAlgo);
    w.Simulation(stopCond);
    stat << perform.ffct << ", " << perform.fsct << ", " << perform.fuft << ", " << perform.faft << ", " << perform.fant << ", " << sym_time.elapsed() << std::endl;
    {
        std::ofstream stat_cp(stat_name, std::ios::app);
        stat_cp << stat.str();
        stat_cp.close();
    }
    shifts.close();
    sym_time.reset();

    for (auto i: id) w.obj(i).save(res_path + prefix + w.obj(i).name + postfix + ".vtk");

    if (view) {
        view_window->join();
        delete view_window;
    }

    return 0;
}

static AniMesh generateLeaf(double R, double d, double mesh_h){
    auto len_func = [R, d](double t)->double{
        auto bell = [m = 3*M_PI_2, w = M_PI_4, r = M_PI_2](double x){
            return (exp(- (x-m)*(x-m) / (w * w)) -  exp(- (r * r) / (w * w))) / (1 - exp(- (r * r) / (w * w) ));
        };
        return R*abs(sin(t)) + d*bell(t);
    };
    auto curve = [R](double t) -> Point{ return Point{R*cos(t), R*sin(t), 0}; };
    auto asd = StretchFigure()
            .setStretchDirection({0, 1, 0})
            .setStretchLengthFunc(len_func)
            .setCurveFunc(curve)
            .setPntParams({M_PI, 2*M_PI}, {4, 4, 4})
            .prepareFigure();
    AniMesh am;
    Ani3dMeshOut amo(am);
    makeAft3dSurface(asd.getAniSurfDiscr(), Ani3dFSize(mesh_h).fsize, amo);
    mark_vertices(am);
    return am;
}

int threeLeafContact(int argc, char* argv[], ProblemDiscr p){
    string res_path = p.main_path + "ThreeLeafTest/";
    string gendir = p.gendir;
    bool view = p.with_view;

    if (p.geometry.size() < 2) throw std::runtime_error("Geometry should consist 2 vals");
    double H0 = p.H0, R = p.geometry[0], d = p.geometry[1], mesh_h = p.mesh_h;
    double P = p.P, mu = p.mu;//50+
    double delta = p.delta, err = p.err, diverge_lim = p.diverge_lim;
    int pfreq = p.pfreq, maxits = p.maxits;
    double bullet_margin = p.margin;
    bool is_penalty = p.is_penalty;
    if (is_penalty) bullet_margin *= 2;

    if (!exists(res_path)) create_directory(res_path);
    std::array<Object3D, 3> _obj;
    {
        std::array<AniMesh, 3> am;
        am[0] = generateLeaf(R, d, mesh_h);
        am[1] = am[2] = am[0];
        double Rc = 3 * (2*R+H0) / (2 * M_PI) * (1 + 1e-5);
        for (int i = 0; i < 3; ++i) {
            for (int n = 0; n < am[i].vertices.size() / 3; ++n) {
                double x0 = am[i].vertices[3 * n + 0], y0 = am[i].vertices[3 * n + 1];
                double phi = (x0 / Rc + 2 * M_PI / 3 * i);
                am[i].vertices[3 * n + 0] = Rc * cos(phi);
                am[i].vertices[3 * n + 1] = Rc * sin(phi);
                am[i].vertices[3 * n + 2] = y0;
            }
        }

        for (int i = 0; i < 3; ++i) {
            _obj[i].m_mesh = get_invert_Mesh(convert_to_Mesh(am[i], "v:boundary_lbl"));
            _obj[i].reset_mesh();
            _obj[i].name = "leaf" + to_string(i);
//            _obj[i].save(res_path + _obj[i].name + "_init.vtk");
//            exit(-1);
        }
    }
    static const int CLAMPED = 4, FREE = 1;
    auto dir_cond = [](Object3D& obj, const V_ind& v) -> unsigned char{ return (obj.m_boundary[v] & CLAMPED) ? 0 : 7; };
    World w;
    std::array<ObjectID, 3> id = {-1};
    for (int i = 0; i < 3; ++i) id[i] = w.addObject3D(move(_obj[i]), 1);
    for (auto i: id) w.obj(i).setDirichletCond(dir_cond);
    SimplePressureLoad Pr(P);
    Force elastic = NeoGookModelOpt(mu, H0);
    w.addForce(elastic);
    w.addForce(Pr);
    struct PerformTimers{
        double find_col_tm = 0, sol_col_tm = 0, updt_f_tm = 0, appl_f_tm = 0, appl_nxt_tm = 0;
        double ffct = 0, fsct = 0, fuft = 0, faft = 0, fant = 0;
        World3d::Timer m_timer;
    };
    PerformTimers perform;

    auto stepAlgo = [&p = perform](World* w)->int{
        int status = 0;
        p.m_timer.reset();
        if ((status = w->UpdateForces()) < 0) return status;
        p.updt_f_tm += p.m_timer.elapsed(); p.m_timer.reset();
        if ((status = w->ApplyForces()) < 0) return status;
        p.appl_f_tm += p.m_timer.elapsed(); p.m_timer.reset();

        w->getCollider()->findCollisions();
        p.find_col_tm += p.m_timer.elapsed(); p.m_timer.reset();
        w->getCollider()->solveCollisions();
        p.sol_col_tm += p.m_timer.elapsed(); p.m_timer.reset();

        if ((status = w->ApplyNext()) < 0) return status;
        p.appl_nxt_tm += p.m_timer.elapsed(); p.m_timer.reset();
        return status;
    };

    string prefix = "bul_constr_", postfix;
    if (is_penalty) prefix = "bul_penalty_";
    {
        std::ostringstream oss;
        auto tt = std::time(nullptr);
        auto tm = *std::localtime(&tt);
        oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
        postfix = "_" + oss.str();
    }
    if (p.only_save_input){
        for (auto i: id) w.obj(i).save(res_path + w.obj(i).name + "_R=" + to_string(R) + "_d=" + to_string(d) + ".vtk");
        return 0;
    }
    std::ofstream shifts(res_path + prefix + "shift" + postfix + ".csv");
    shifts << "it, delta_x, find_col_time, sol_col_time, update_force_time, apply_force_time, apply_next_time\n";
    string stat_name = res_path + "stat" + ".cvs";
    bool is_new = exists(stat_name);
    {
        std::ofstream stat(stat_name, std::ios::app);
        if (!is_new)
            stat << "prefix, postfix, H0, R, d, mh, P, mu, delta, maxit, param"
                 << ", find_col_time, sol_col_time, update_force_time, apply_force_time, apply_next_time, full_time"
                 << std::endl;
        stat.close();
    }
    std::ostringstream stat;
    stat << prefix << ", " << postfix << ", "
         << H0 << "," << R  << ", " << d << ", "<< mesh_h << ", "
         << P << ", " << mu << ", " << delta << ", " << maxits << ", {margin=" << bullet_margin << "}, ";
//    unique_ptr<CollisionManagerBase> colMan = std::make_unique<CGALColManVFfull>();
//    w.setCollider(move(colMan));
//    for (auto i: id) reinterpret_cast<CGALColManVFfull*>(w.getCollider())->set_margin(i, bullet_margin);

    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BulletCollisionManager>();
    if (is_penalty)
        reinterpret_cast<BulletCollisionManager*>(colMan.get())->setCustomSolveColAlgo(PenaltyCollisionAlgo());
    w.setCollider(move(colMan));
    for (auto i: id) reinterpret_cast<BulletCollisionManager*>(w.getCollider())->set_margin(i, bullet_margin);

    std::thread *view_window;
    if (view) {
        view_window = new std::thread(start_debug_gui, argc, argv);
        w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    }

    StaticForceApplierP sfa(delta);
    //sfa.addDamper(Damper(0.01*H, 0.8, 100, 100000).setReInitDxParams(1e7, 0.4));
    w.setForceApplier(sfa);

    World3d::Timer sym_time;
    auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &maxits, &pfreq, &sym_time, nullptr, true);
    stopCond.addExternalCall([freq = pfreq, &shifts, &obj0 = w.obj(id[0]), &p = perform](StepSimInfo& info, World *w) -> bool {
        p.ffct += p.find_col_tm; p.fsct += p.sol_col_tm; p.fuft += p.updt_f_tm; p.faft += p.appl_f_tm; p.fant += p.appl_nxt_tm;
        static double resid_init = 1.0;
        double resid = w->getWorldShift();

        if (info.it > 0) shifts << info.it << ", " << resid
                                << ", " << p.find_col_tm << ", " << p.sol_col_tm << ", " << p.updt_f_tm << ", " << p.appl_f_tm << ", " << p.appl_nxt_tm << "\n";
        p.find_col_tm = p.sol_col_tm = p.updt_f_tm = p.appl_f_tm = p.appl_nxt_tm = 0;
        if (info.it == 1) resid_init = w->getWorldShift();
        if (info.it > 0 && info.it % freq  == 0) {
            double eps = resid / resid_init;
            std::cout << "Shift[" << info.it << "] = " << resid << " (" << eps << ")" << "\n";
        }
        return true;
    });
    w.setStepSimulationAlgo(stepAlgo);
    w.Simulation(stopCond);
    stat << perform.ffct << ", " << perform.fsct << ", " << perform.fuft << ", " << perform.faft << ", " << perform.fant << ", " << sym_time.elapsed() << std::endl;
    {
        std::ofstream stat_cp(stat_name, std::ios::app);
        stat_cp << stat.str();
        stat_cp.close();
    }
    shifts.close();
    sym_time.reset();

    for (auto i: id) w.obj(i).save(res_path + prefix + w.obj(i).name + postfix + ".vtk");

    if (view) {
        view_window->join();
        delete view_window;
    }

    return 0;
}


int OverlapContact(int argc, char* argv[], ProblemDiscr p){
    string res_path = p.main_path + "OverlapTest/";
    string gendir = p.gendir;
    bool view = p.with_view;

    if (p.geometry.size() < 4) throw std::runtime_error("Geometry should consist 4 vals");
    double H0 = p.H0, R = p.geometry[0], d[2] = {p.geometry[1], p.geometry[2]}, dist =  p.geometry[3], mesh_h = p.mesh_h;
    double P = p.P, mu = p.mu;//50+
    double delta = p.delta, err = p.err, diverge_lim = p.diverge_lim;
    int pfreq = p.pfreq, maxits = p.maxits;
    double bullet_margin = p.margin;
    bool is_penalty = p.is_penalty;
    if (is_penalty) bullet_margin *= 2;

    if (!exists(res_path)) create_directory(res_path);
    std::array<Object3D, 2> _obj;
    {
        std::array<AniMesh, 2> am;
        am[0] = generateLeaf(R, d[0], mesh_h);
        am[1] = generateLeaf(R, d[1], mesh_h);
        double Rc = 3 * (2*R+H0) / (2 * M_PI) * (1 + 1e-5);
        for (int i = 0; i < 2; ++i) {
            for (int n = 0; n < am[i].vertices.size() / 3; ++n) {
                double x0 = am[i].vertices[3 * n + 0], y0 = am[i].vertices[3 * n + 1];
                double phi = (x0 / Rc + M_PI * i);
                am[i].vertices[3 * n + 0] = Rc * cos(phi) + (i ? 1 : -1) *((Rc - H0)/2 - dist/2);
                am[i].vertices[3 * n + 1] = Rc * sin(phi);
                am[i].vertices[3 * n + 2] = y0;
            }
        }

        for (int i = 0; i < 2; ++i) {
            _obj[i].m_mesh = get_invert_Mesh(convert_to_Mesh(am[i], "v:boundary_lbl"));
            _obj[i].reset_mesh();
            _obj[i].name = "leaf" + to_string(i);
//            _obj[i].save(res_path + _obj[i].name + "_init.vtk");
//            exit(-1);
        }
    }
    static const int CLAMPED = 4, FREE = 1;
    auto dir_cond = [](Object3D& obj, const V_ind& v) -> unsigned char{ return (obj.m_boundary[v] & CLAMPED) ? 0 : 7; };
    World w;
    std::array<ObjectID, 2> id = {-1};
    for (int i = 0; i < id.size(); ++i) id[i] = w.addObject3D(move(_obj[i]), 1);
    for (auto i: id) w.obj(i).setDirichletCond(dir_cond);
    SimplePressureLoad Pr(P);
    Force elastic = NeoGookModelOpt(mu, H0);
    w.addForce(elastic);
    w.addForce(Pr);
    struct PerformTimers{
        double find_col_tm = 0, sol_col_tm = 0, updt_f_tm = 0, appl_f_tm = 0, appl_nxt_tm = 0;
        double ffct = 0, fsct = 0, fuft = 0, faft = 0, fant = 0;
        World3d::Timer m_timer;
    };
    PerformTimers perform;

    auto stepAlgo = [&p = perform](World* w)->int{
        int status = 0;
        p.m_timer.reset();
        if ((status = w->UpdateForces()) < 0) return status;
        p.updt_f_tm += p.m_timer.elapsed(); p.m_timer.reset();
        if ((status = w->ApplyForces()) < 0) return status;
        p.appl_f_tm += p.m_timer.elapsed(); p.m_timer.reset();

        w->getCollider()->findCollisions();
        p.find_col_tm += p.m_timer.elapsed(); p.m_timer.reset();
        w->getCollider()->solveCollisions();
        p.sol_col_tm += p.m_timer.elapsed(); p.m_timer.reset();

        if ((status = w->ApplyNext()) < 0) return status;
        p.appl_nxt_tm += p.m_timer.elapsed(); p.m_timer.reset();
        return status;
    };

    string prefix = "bul_constr_", postfix;
    if (is_penalty) prefix = "bul_penalty_";
    {
        std::ostringstream oss;
        auto tt = std::time(nullptr);
        auto tm = *std::localtime(&tt);
        oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
        postfix = "_" + oss.str();
    }
    if (p.only_save_input){
        for (auto i: id) w.obj(i).save(res_path + w.obj(i).name
            + "_R=" + to_string(R) + "_d0=" +to_string(d[0]) + "_d1=" + to_string(d[1]) + "_dist=" + to_string(dist) + ".vtk");
        return 0;
    }
    std::ofstream shifts(res_path + prefix + "shift" + postfix + ".csv");
    shifts << "it, delta_x, find_col_time, sol_col_time, update_force_time, apply_force_time, apply_next_time\n";
    string stat_name = res_path + "stat" + ".cvs";
    bool is_new = exists(stat_name);
    {
        std::ofstream stat(stat_name, std::ios::app);
        if (!is_new)
            stat << "prefix, postfix, H0, R, d0, d1, dist, mh, P, mu, delta, maxit, param"
                 << ", find_col_time, sol_col_time, update_force_time, apply_force_time, apply_next_time, full_time"
                 << std::endl;
        stat.close();
    }
    std::ostringstream stat;
    stat << prefix << ", " << postfix << ", "
         << H0 << "," << R  << ", " << d[0] << ", "<< d[1] << ", " << dist << ", " << mesh_h << ", "
         << P << ", " << mu << ", " << delta << ", " << maxits << ", {margin=" << bullet_margin << "}, ";
//    unique_ptr<CollisionManagerBase> colMan = std::make_unique<CGALColManVFfull>();
//    w.setCollider(move(colMan));
//    for (auto i: id) reinterpret_cast<CGALColManVFfull*>(w.getCollider())->set_margin(i, bullet_margin);

    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BulletCollisionManager>();
    if (is_penalty)
        reinterpret_cast<BulletCollisionManager*>(colMan.get())->setCustomSolveColAlgo(PenaltyCollisionAlgo());
    w.setCollider(move(colMan));
    for (auto i: id) reinterpret_cast<BulletCollisionManager*>(w.getCollider())->set_margin(i, bullet_margin);

    std::thread *view_window;
    if (view) {
        view_window = new std::thread(start_debug_gui, argc, argv);
        w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    }

    StaticForceApplierP sfa(delta);
    //sfa.addDamper(Damper(0.01*H, 0.8, 100, 100000).setReInitDxParams(1e7, 0.4));
    w.setForceApplier(sfa);

    World3d::Timer sym_time;
    auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &maxits, &pfreq, &sym_time, nullptr, true);
    stopCond.addExternalCall([freq = pfreq, &shifts, &obj0 = w.obj(id[0]), &p = perform](StepSimInfo& info, World *w) -> bool {
        p.ffct += p.find_col_tm; p.fsct += p.sol_col_tm; p.fuft += p.updt_f_tm; p.faft += p.appl_f_tm; p.fant += p.appl_nxt_tm;
        static double resid_init = 1.0;
        double resid = w->getWorldShift();

        if (info.it > 0) shifts << info.it << ", " << resid
                                << ", " << p.find_col_tm << ", " << p.sol_col_tm << ", " << p.updt_f_tm << ", " << p.appl_f_tm << ", " << p.appl_nxt_tm << "\n";
        p.find_col_tm = p.sol_col_tm = p.updt_f_tm = p.appl_f_tm = p.appl_nxt_tm = 0;
        if (info.it == 1) resid_init = w->getWorldShift();
        if (info.it > 0 && info.it % freq  == 0) {
            double eps = resid / resid_init;
            std::cout << "Shift[" << info.it << "] = " << resid << " (" << eps << ")" << "\n";
        }
        return true;
    });
    w.setStepSimulationAlgo(stepAlgo);
    w.Simulation(stopCond);
    stat << perform.ffct << ", " << perform.fsct << ", " << perform.fuft << ", " << perform.faft << ", " << perform.fant << ", " << sym_time.elapsed() << std::endl;
    {
        std::ofstream stat_cp(stat_name, std::ios::app);
        stat_cp << stat.str();
        stat_cp.close();
    }
    shifts.close();
    sym_time.reset();

    for (auto i: id) w.obj(i).save(res_path + prefix + w.obj(i).name + postfix + ".vtk");

    if (view) {
        view_window->join();
        delete view_window;
    }

    return 0;
}

int ContactTest(int argc, char* argv[]){
    auto p = ProblemDiscr::Init(argc, argv);
    switch (p.tt) {
        case ProblemDiscr::Ellipsoid: return ellipsoidContact(argc, argv, p);
        case ProblemDiscr::ThreeLeaf: return threeLeafContact(argc, argv, p);
        case ProblemDiscr::Overlap: return OverlapContact(argc, argv, p);
    }
    return -1;
}