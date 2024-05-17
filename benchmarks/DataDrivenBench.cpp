//
// Created by alex on 02.02.2021.
//

#include <iostream>
#include <string>
#include <stdexcept>
#include <Eigen/Dense>
#include "Benchmarks.h"
#include "TriangularMeshHelpers.h"

#ifdef USE_MATIO
#include <matio.h>
#endif

using namespace World3d;

static int start_gui(int argc, char** argv){
#ifdef USE_MAGNUM_GUI
    GuiApplication app({argc, argv});
    return app.exec();
#else
    return 0;
#endif
}

static bool ObSaveTag(const Object3D& obj, string filename, string tag){
    ofstream ob(filename, ios::binary | ios::trunc);
    if (!ob) return false;
    auto m_x = obj.m_mesh.property_map<V_ind, Point>(tag).first;
    for (auto& v: obj.m_mesh.vertices()){
        std::array<double, 3> P = {m_x[v][0], m_x[v][1], m_x[v][2]};
        ob.write(reinterpret_cast<const char *> ((P.data())), 3 * sizeof(double));
    }
    ob.close();
    return true;
};

static bool ObReadTag(Object3D& obj, string filename, string tag){
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
};

int testDataDriven(int argc, char* argv[]){
    thread t(start_gui, argc, argv);
    AniMesh am = generate_circle(1, 0.05);
    for (auto& i: am.vertex_label) i &= 1;
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = "circle";
    translate_Mesh(m, {0, 0, 2.0});
    m = get_invert_Mesh(m);
    Object3D obj1{m};
    obj1.name = "inv_circle";
    World world;
    world.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BulletCollisionManager>();
    world.setCollider(move(colMan));
    auto id = world.addObject3D(move(obj), 1);
    auto id1 = world.addObject3D(move(obj1), 1);
    reinterpret_cast<BulletCollisionManager*>(world.getCollider())->set_margin(0.01);
    SimpleMassSpringModel MSM(3.0e5, 5e-3);
    SimplePressureLoad Pr(20.0_mmHg);
    std::string dir = "../generated";
    NeoGookModel ngk(1.0e6/3, 5e-3, dir, true);
    SVKirchhoffModel svk(5e-3*4, 0, 1e6/2, dir, true);
//    BendingForce bsvk(svk.f, true);
//    DataDrivenAnalytic svk_twink(svk.f, true);
    DataDrivenAnalytic ngk_twink(ngk.f, true);
    world.addForce(std::move(ngk), id);
    world.addForce(std::move(ngk_twink), id1);

    world.addForce(std::move(Pr));
    world.setForceApplier(StaticForceApplier(2.9e-5/4));

    world.Simulation([](StepSimInfo& info, World* w)->bool{
        int it = info.it;
        if (it % 1000 == 0) std::cout << "it = " << it << std::endl;
        if (it >= 20000) return true;
        return false;
    });
    t.join();
    return 0;
}

template<typename TMat>
void gatherDefGrad(Object3D& obj, TMat& F_d){
    auto Lp = set_L(obj);
    auto Nel = obj.m_mesh.num_faces();
    for (auto f: obj.m_mesh.faces()) {
        auto vv = vert_around(obj.m_mesh, f);
        auto L = Eigen::Map<Eigen::Matrix<World3d::DReal, 2, 3>>(Lp[f].data());
        Eigen::Matrix<World3d::DReal, 3, 3> Q;
        Q <<    obj.m_x[vv[0]][0], obj.m_x[vv[1]][0], obj.m_x[vv[2]][0],
                obj.m_x[vv[0]][1], obj.m_x[vv[1]][1], obj.m_x[vv[2]][1],
                obj.m_x[vv[0]][2], obj.m_x[vv[1]][2], obj.m_x[vv[2]][2];
        Eigen::Matrix<World3d::DReal, 3, 2> F = Q * L.transpose();
        F_d(f.idx() + 0 * Nel, 0) = F(0, 0), F_d(f.idx() + 2 * Nel, 0) = F(0, 1);
        F_d(f.idx() + 3 * Nel, 0) = F(1, 0), F_d(f.idx() + 1 * Nel, 0) = F(1, 1);
    }
}

struct SigmaStabilityOptions{
    std::string res_path = "../result/DataDrivenBench/StabilityTest/";
    //mesh options
    double R = 10, mesh_h = 0.5;
    //direct modeling options
    int get_scenario = 0;
    double P = 5.0e3;
    double mu = 1.0e6, H = 5e-2;
    double delta = 1.5e-6;
    int freq = 1000, maxits = 220000;
    double err = 1e-7;
    std::string dir = "../generated";
    //inverse modeling options
    int inv_scenario = 0;
    double inv_mu = 1.0e12;
    double inv_delta = 1.5e-12;
    int inv_freq = freq, inv_maxits = maxits;
    double inv_err = 1e-6;

    Force choose_elastic(int scenario, double mu) {
        Force elastic;
        switch (scenario) {
            case 0:
                elastic = NeoGookModel(mu, H, dir, false);
                break;
            case 1:
                elastic = SVKirchhoffModel(H, 0, mu, dir, false);
                break;
            case 2:
                elastic = MooneyRivlinModel2(H, mu / 4, 0.05, 0.3, {1, 0, 0}, dir, false);
                break;
            default:
                throw std::runtime_error("Faced undefined scenario");
        }
        return elastic;
    };
    static std::string get_short_pref(int scenario){
        switch (scenario) {
            case 0: return "ngk";
            case 1: return "svk";
            case 2: return "mrm";
            default: throw std::runtime_error("Faced undefined scenario");
        };
        return "";
    };
    static std::function<bool(StepSimInfo&, World* w)> stopCondition(int& freq, double& err, World3d::Timer& time, int& maxits){
        return [&freq, &err, &time, &maxits](StepSimInfo& info, World* w)->bool{
            static double resid_init;
            int it = info.it;
            if (it == 1) resid_init = w->getWorldNorm("v:force");
            if (it > 0 && it % freq == 0) {
                double resid = w->getWorldNorm("v:force");
                double eps = resid / resid_init;
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                if (eps < err){
                    std::cout << "Algorithm is converged: \n";
                    return true;
                }
                if (it > maxits || std::isnan(eps)){
                    std::cout << "Algorithm is diverged or reached maximum of time - iters: \n";
                    return true;
                }
            }
            return false;
        };
    };
};

int processGeneratedState(int argc, char* argv[], SigmaStabilityOptions so = SigmaStabilityOptions()){
    static auto readTag = [](Object3D& obj, string filename, string tag){
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
    };

    Object3D robj;
    robj.read(so.res_path + "neogook_input"  + ".txt");
    readTag(robj, so.res_path + "neogook_input"  + "_m_x0.tag", "v:point");
    readTag(robj, so.res_path + "neogook_input"  + "_m_x.tag", "v:x");
    ResponseStatistics rs;
    Force elastic = std::move(so.choose_elastic(so.get_scenario, so.mu));
    HEDataGatherer gather(elastic.target<HyperElasticForceBase>()->f, true);
    gather.gatherData(robj, rs);

    StaticDefiniteGather::Traits traits;
    traits.setP(so.P).setH(so.H).setMu(so.inv_mu).setDelta(so.inv_delta).setErr(so.inv_err);
    traits.setSolveMethod(StaticDefiniteGather::Traits::SIMPLE_NEWTON).setMaxIters(100).setFreq(1);
    StaticDefiniteGather statgath(robj, traits);
    statgath.computeTensionedState().gatherResponseData(rs);
    std::string postfix = "";//"_JT";
    rs.save_distribution(so.res_path + "gather_" + "ngk" + "-" + so.get_short_pref(so.inv_scenario) + postfix + ".dist");

#ifdef USE_INMOST_SOLVERS
    {
        auto& ob = robj;
        using namespace INMOST;
        INMOST::Mesh* mesh = new INMOST::Mesh();
        vector<HandleType> new_nodes(ob.m_mesh.num_vertices());
        vector<HandleType> new_faces(ob.m_mesh.num_faces());
        for (auto v: ob.m_mesh.vertices()){
            std::array<Storage::real,3> xyz{(Storage::real)ob.m_x0[v][0], (Storage::real)ob.m_x0[v][1], (Storage::real)ob.m_x0[v][2]};
            new_nodes[v.idx()] = mesh->CreateNode(xyz.data())->GetHandle();
        }
        for (auto f: ob.m_mesh.faces()){
            auto v = vert_around(ob.m_mesh, f);
            ElementArray<Node> subarr(mesh);
            subarr.reserve(3);
            for (int j = 0; j < 3; ++j)
                subarr.push_back(new_nodes[v[j].idx()]);
            const int nodesnum[3] = {0,1,2};
            const int sizes[1] = {3};

            new_faces[f.idx()] = mesh->CreateCell(subarr, nodesnum, sizes, 1).first.GetHandle();
        }
        const int ntag = 4;
        std::array<Tag, ntag> kk = {
                mesh->CreateTag("|(R-Rcomp)/R|", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("R", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("R-Rcomp", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("Rcomp", DATA_REAL, NODE, FACE, 3)
        };

        for (auto n = mesh->BeginNode(); n != mesh->EndNode(); ++n){
            auto ff = face_around(ob.m_mesh, V_ind(n->LocalID()));
            std::array<std::array<double, 3>, ntag> res{0};
            for (auto f: ff)
                for (int k = 0; k < 3; ++k) {
                    res[0][k] += fabs((rs.stat[f.idx()][0].response[k] - rs.stat[f.idx()][1].response[k]) /
                                   rs.stat[f.idx()][0].response[k]);
                    res[1][k] += rs.stat[f.idx()][0].response[k];
                    res[2][k] += rs.stat[f.idx()][0].response[k] - rs.stat[f.idx()][1].response[k];
                    res[3][k] += rs.stat[f.idx()][1].response[k];
                }
            for (int i = 0; i < 3 && !ff.empty(); ++i) {
                for (int j = 0; j < ntag; ++j)
                    res[j][i] /= ff.size();
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < ntag; ++j)
                    n->RealArray(kk[j])[i] = res[j][i];
            }
        }

        mesh->SaveVTK(so.res_path + "neogook_res_dR" + postfix + "_flat.vtk");
        delete mesh;
    }
#endif

    return 0;
}

int testSigmaStability(int argc, char* argv[], SigmaStabilityOptions so = SigmaStabilityOptions()){
//    thread t(start_gui, argc, argv);

    AniMesh am = generate_circle(so.R, so.mesh_h);
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = "circle";
    const int CENTER = 2|4;
    auto is_movable = [] (Object3D& obj, const V_ind& v) -> bool { return !(obj.m_boundary[v] & 1); };
    obj.set_is_movable_checker(is_movable);

    World world;
//    world.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    ObjectID id = world.addObject3D(move(obj), 1);
//    V_ind vc;
//    for (auto v: world.obj(id).m_mesh.vertices())
//        if (world.obj(id).m_boundary[v] == CENTER) vc = v;
    SimplePressureLoad Pr(so.P);

    Force elastic = std::move(so.choose_elastic(so.get_scenario, so.mu));
    world.addForce(elastic, id);
    world.addForce(Pr, id);
    world.setForceApplier(StaticForceApplier(so.delta));

    auto time = World3d::Timer();
    auto stopCondition = [](auto& freq, auto& err, auto& time, auto& vc, auto& id, auto& maxits){
        return [&freq, &err, &time, &vc, &id, &maxits](StepSimInfo& info, World* w)->bool{
            static double resid_init;
            int it = info.it;
            if (it == 1) resid_init = w->getWorldNorm("v:force");
            if (it > 0 && it % freq == 0) {
                double resid = w->getWorldNorm("v:force");
                double eps = resid / resid_init;
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                std::cout << "watch { " << w->obj(id).m_x[vc].z() << " }" << std::endl;
                if (eps < err){
                    std::cout << "Algorithm is converged: \n";
                    return true;
                }
                if (it > maxits || std::isnan(eps)){
                    std::cout << "Algorithm is diverged or reached maximum of time - iters: \n";
                    return true;
                }
            }
            return false;
        };
    };
    world.Simulation(so.stopCondition(so.freq, so.err, time, so.maxits));
    static auto saveTag = [](const Object3D& obj, string filename, string tag){
        ofstream ob(filename, ios::binary | ios::trunc);
        if (!ob) return false;
        auto m_x = obj.m_mesh.property_map<V_ind, Point>(tag).first;
        for (auto& v: obj.m_mesh.vertices()){
            std::array<double, 3> P = {m_x[v][0], m_x[v][1], m_x[v][2]};
            ob.write(reinterpret_cast<const char *> ((P.data())), 3 * sizeof(double));
        }
        ob.close();
        return true;
    };
    static auto readTag = [](Object3D& obj, string filename, string tag){
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
    };
    world.obj(id).save(so.res_path + "neogook_input"  + ".txt");
    saveTag(world.obj(id), so.res_path + "neogook_input"  + "_m_x0.tag", "v:point");
    saveTag(world.obj(id), so.res_path + "neogook_input"  + "_m_x.tag", "v:x");

//    t.join();

    return 0;
}

Eigen::MatrixXd read_mat(std::string name, bool print = false){
    Eigen::MatrixXd mat;
#ifdef USE_MATIO
    mat_t *matfp;
    matvar_t *matvar;
    matfp = Mat_Open(name.c_str(),MAT_ACC_RDONLY);
    if (!matfp ) throw std::runtime_error("Error while opening MAT file \"" + name + "\"");
    //read only first variable
    matvar = Mat_VarReadNextInfo(matfp);
    if (!matvar) throw std::runtime_error("Error while reading MAT file \"" + name + "\"");
    if (matvar->data_type != MAT_T_DOUBLE && matvar->class_type != MAT_C_DOUBLE)
        throw std::runtime_error("Variable in MAT file has unsupported type");

    mat.resize(matvar->dims[0], matvar->dims[1]);
    int    start[2]={0,0};
    int    stride[2]={1,1};
    int    edge[2];
    edge[0] = matvar->dims[0], edge[1] = matvar->dims[1];

    Mat_VarReadData(matfp, matvar, mat.data(), start, stride, edge);

    if (print) Mat_VarPrint(matvar, 0), std::cout << "\n";
//        Mat_VarFree(matvar);
    Mat_Close(matfp);
#endif
    return mat;
}

class HGOModel1: public HyperElasticForceBase {
public:
    DReal _h, _c, _k1, _k2, _kappa;
    std::array<DReal, 3> _Av;
    HGOModel1(double h, double c, double k1, double k2, double kappa, std::array<double, 3> _Av, std::string gen_dir, bool regenerate = true):
            _h{h}, _c{c}, _k1{k1}, _k2{k2}, _kappa{kappa}, _Av{_Av}
    {
        type = "HGOModel1";
        SX H  = SX::sym("H" );
        SX C = SX::sym("C");
        SX K1  = SX::sym("k1" );
        SX K2  = SX::sym("k2" );
        SX kap = SX::sym("kap");
        auto lv = DefaultHyperElasticForce::getDefaultLocalVars();
        SX Ap = lv["Ap"];
        auto I = DefaultHyperElasticForce::getDefaultInvariants();
        std::array<DReal, 3> Av{_Av[0], _Av[1], _Av[2]};
        auto orth = DefaultHyperElasticForce::makeOrthotropicInvariant("I4s", Av);
        SX I1 = I["I1"];
        SX I3 = I["I3"];
        SX I4s = get<2>(orth)["I4s"];
        SX E = _kappa*(I1 + 1.0/I3) + (1 - 3*_kappa)*I4s;
        SX U = Ap * H * (C/2 * (I1 + 1.0/I3 - 3) + K1/(4*K2) * SX::exp(K2*SX::sq(E - 1)) - 1);
        auto iv = DefaultHyperElasticForce::getDefaultAvailableInputVars();
        iv.emplace_back("H", H, std::make_unique<ConstDataRequier>(_h));
        iv.emplace_back("C", C, std::make_unique<ConstDataRequier>(_c));
        iv.emplace_back("k1", K1, std::make_unique<ConstDataRequier>(_k1));
        iv.emplace_back("k2", K2, std::make_unique<ConstDataRequier>(_k2));
        iv.emplace_back("kap", kap, std::make_unique<ConstDataRequier>(_kappa));

        HH::ArgsDef args(std::move(iv), std::move(lv), std::move(I));
        args = HH::concat(args, orth);
        f = DefaultHyperElasticForce("HGOModel1", gen_dir, std::move(args), U, regenerate);
    }
    std::unique_ptr<HyperElasticForceBase> copy() override{
        return std::make_unique<HGOModel1>(_h, _c, _k1, _k2, _kappa, _Av, f.gen_dir, false);
    }
};

int testCircleStatDef(int argc, char* argv[]){
    thread t(start_gui, argc, argv);

    double R0 = 4, H0 = 0.004, lambda = 1.25, mesh_h = 0.1/2;
    double R = lambda * R0;
    double P = 2 / R;
    double mu = P * R / (2 * H0);
    double alpha = 0.0245, beta = 0.0297, teta = 0.347;
    int freq = 1, maxits = 100, status = 0;
    double err = 1e-4;
    std::string postfix = "neo_fine_" + std::to_string(static_cast<int>(100*lambda));
    string res_path = "../result/StatDefSph/Circle";
    bool resimulate = true;
    bool save_meta_info = false;
    Object3D obj;
    if (resimulate) {
        AniMesh am = generate_circle(R0, mesh_h);
        obj = Object3D{ convert_to_Mesh(am, "v:boundary_lbl") };
    } else {
        obj.read(res_path + "/circle_" + postfix + ".txt");
        if (!ObReadTag(obj, res_path + "/circle_x0_" + postfix + ".txt", "v:point"))
            throw std::runtime_error("Can't open " + res_path + "/circle_x0_" + postfix + ".txt");
//        ObReadTag(obj, res_path + "/circle_x_" + postfix + ".txt", "v:x");
    }
    obj.name = "circle";
    auto is_movable = [] (Object3D& obj, const V_ind& v) -> bool { return !(obj.m_boundary[v] & 1); };
    obj.set_is_movable_checker(is_movable);

//    Force elastic = MooneyRivlinModel(H0, mu, alpha, /*beta, {cos(teta), sin(teta), 0},*/ "../generated", true);
//    Force elastic = MooneyRivlinModel(H0, mu, alpha, /*beta, {cos(teta), sin(teta), 0},*/ "../generated", true);
    Force elastic = NeoGookModel(mu, H0, "../generated", true);
    elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(true);
    SimplePressureLoad Pr(P);

    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto oid = w.addObject3D(std::move(obj), 1);
    auto efid = w.addForce(std::move(elastic), oid);
    auto pid = w.addForce(std::move(Pr), oid);
    NewtonSolverStepAlgo nssa;
#ifdef USE_INMOST_SOLVERS
    nssa._solver.solver = std::move(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#else
    nssa._solver.solver = std::move(LinearSolver("eigen"));
#endif
    nssa._solver.solver.setVerbosityLevel(0);
    nssa.newton_algo = [](NewtonSolverStepAlgo& ns, World* w){
        ns.updateMeshData(w);
        int N = ns.getNodfs(w);
        ns.assembleSystem(w, N);
        double tau = 0.1;//*1e-4;//0.3;
        static double resid_init = -1;
        if (ns.m_it == 0) {
            resid_init = ns.computeRhsNorm(1);
            tau = 1e-4;
        } else {
            double eps = ns.computeRhsNorm(1) / resid_init;
            if (ns.m_it < 23){   //initial relaxation
                if (ns.m_it < 5) tau = 5e-4;
                else if (ns.m_it < 10) tau = 1e-2;
                else if (ns.m_it < 15) tau = 3e-2;
                else tau = 6e-2;
            }
//            if (ns.m_it < 5){   //initial relaxation
//                if (ns.m_it == 2) tau = 4e-3;
//                else tau = 1e-2;
//            }
            else if (eps > 100){  //stability borders
                tau = 1e-1;
            } else if (eps < 3){
                tau = 1;
            } else {
                tau = 0.3;
            }
        }
        ns.Solve();
        ns.applySolution(w, tau);
        ns.m_it++;
        return 0;
    };
//    w.setForceApplier(StaticForceApplier(5e-2));
//    freq = 500; maxits = 1e5;

    auto time = World3d::Timer();
    if (resimulate) {
        int temp_freq = 500, temp_maxit = 1000;
        swap(freq, temp_freq), swap(maxits, temp_maxit);
        w.setForceApplier(StaticForceApplier(5e-2));
        auto stopCond = [&freq, &err, &time, &maxits, &status](StepSimInfo& info, World *w) -> bool {
            static double resid_init;
            int it = info.it;
            if (it == 1) resid_init = w->getWorldNorm("v:force");
//            if (it > 0) sleep(10);
            if (it > 0 && it % freq == 0) {
                double resid = w->getWorldNorm("v:force");
                double eps = resid / resid_init;
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed()
                          << "\n";
                if (eps < err) {
                    std::cout << "Algorithm is converged: \n";
                    status = 0;
                    return true;
                }
                if (it >= maxits || std::isnan(eps)) {
                    std::cout << "Algorithm is diverged or reached maximum of time - iters: \n";
                    status = 1;
                    return true;
                }
            }
            return false;
        };
        w.Simulation(stopCond); //initial relaxation

        swap(freq, temp_freq), swap(maxits, temp_maxit);
        w.setStepSimulationAlgo(std::move(nssa));
        w.Simulation(stopCond);
        if (save_meta_info) {
            w.obj(oid).save(res_path + "/circle_" + postfix + ".txt");
            ObSaveTag(w.obj(oid), res_path + "/circle_x0_" + postfix + ".txt", "v:point");
            //        ObSaveTag(w.obj(oid), res_path + "/circle_x_" + postfix + ".txt", "v:x");
        }

    }


    StaticDefiniteGather::Traits traits;
    StaticDefiniteGather statgath;
    traits.setH(H0).setMu(1e8).setDelta(6e-12).setErr(1e-6);
    traits.setSolveMethod(StaticDefiniteGather::Traits::SIMPLE_NEWTON).setMaxIters(5).setFreq(1);
    traits.setP(P);
#ifdef USE_INMOST_SOLVERS
    traits.setLinearSolver(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#endif
    traits.m_nssa._solver.solver.setVerbosityLevel(1);
    statgath.setTraits(traits).setObject(w.obj(oid));
    statgath.m_t.m_nssa._solver.solver.SetParameterReal("relative_tolerance", 1e-20);
    statgath.m_t.m_nssa._solver.solver.SetParameterReal("absolute_tolerance", 1e-15);

    HEDataGatherer gather(w.obj(oid).m_forces[efid].target<HyperElasticForceBase>()->f, true);
    int nfaces = w.obj(oid).m_mesh.num_faces();
    Eigen::MatrixXd T(6*nfaces, 1), Td(6*nfaces, 1), Te(3*nfaces, 1), Ted(3*nfaces, 1);
    gather.gatherData(w.obj(oid), [&T, &Te, nf = nfaces](const HEDataGatherer::CompDat &dat) {
        for (int i = 0; i < 6; ++i) T(dat.f.idx() + i * nf, 0) = dat.T[i];
        for (int i = 0; i < 3; ++i) Te(dat.f.idx() + i * nf, 0) = dat.Tee[i];
    });
    statgath.computeTensionedState().gatherMyData([&T = Td, &Te = Ted, nf = nfaces](const StaticDefiniteGather::ObjsData &dat) -> void {
                std::pair<int, int> q3[] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}},
                        q2[] = {{0, 0}, {1, 1}, {1, 0}};
                if (dat.f.idx() < nf) {
                    for (int i = 0; i < 6; ++i) T(dat.f.idx() + i * nf, 0) = dat.T(q3[i].first, q3[i].second);
                    for (int i = 0; i < 3; ++i) Te(dat.f.idx() + i * nf, 0) = dat.Tee(q2[i].first, q2[i].second);
                }
            });

    std::cout   << "||T||_F = "  << T.norm() <<  " ||Td||_F = "  << Td.norm() <<   " ||T - Td||_F = "   << (T - Td).norm()   << "\n"
                << "||Te||_F = " << Te.norm() << " ||Ted||_F = " << Ted.norm() <<  " ||Te - Ted||_F = " << (Te - Ted).norm() << std::endl;
//    std::ofstream fout(res_path + "/see_Te_" + (resimulate ? "sim" : "save") + ".txt");
//    fout << Te << std::endl;
//    fout.close();

#ifdef USE_INMOST_SOLVERS
    {
        auto& ob = w.obj(oid);
        using namespace INMOST;
        INMOST::Mesh* mesh = new INMOST::Mesh();
        vector<HandleType> new_nodes(ob.m_mesh.num_vertices());
        vector<HandleType> new_faces(ob.m_mesh.num_faces());
        for (auto v: ob.m_mesh.vertices()){
            std::array<Storage::real,3> xyz{ob.m_x0[v][0], ob.m_x0[v][1], ob.m_x0[v][2]};
            new_nodes[v.idx()] = mesh->CreateNode(xyz.data())->GetHandle();
        }
        for (auto f: ob.m_mesh.faces()){
            auto v = vert_around(ob.m_mesh, f);
            ElementArray<Node> subarr(mesh);
            subarr.reserve(3);
            for (int j = 0; j < 3; ++j)
                subarr.push_back(new_nodes[v[j].idx()]);
            const int nodesnum[3] = {0,1,2};
            const int sizes[1] = {3};

            new_faces[f.idx()] = mesh->CreateCell(subarr, nodesnum, sizes, 1).first.GetHandle();
        }
        const int ntag = 4;
        const int addtag = 4;
        std::array<Tag, ntag+addtag> kk = {
                mesh->CreateTag("(Eexct-E)/Eexct", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("Eexct", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("Eexct-E", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("E", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("T", DATA_REAL, NODE, FACE, 6),
                mesh->CreateTag("Texct", DATA_REAL, NODE, FACE, 6),
                mesh->CreateTag("Texct-T", DATA_REAL, NODE, FACE, 6),
                mesh->CreateTag("(Texct-T)/Texct", DATA_REAL, NODE, FACE, 6)
        };

        double mx = 0; int ix = -1;
        int n_epoch = 0;
        auto EE = [&Te, nfaces, n = n_epoch](F_ind f, int k) { return Te(f.idx() + k*nfaces, n); };
        auto EEd = [&Ted, nfaces, n = n_epoch](F_ind f, int k) { return Ted(f.idx() + k*nfaces, n); };
        auto TT = [&T, nfaces, n = n_epoch](F_ind f, int k) { return T(f.idx() + (k + 0)*nfaces, n); };
        auto TTd = [&Td, nfaces, n = n_epoch](F_ind f, int k) { return Td(f.idx() + (k + 0)*nfaces, n); };
        for (auto n = mesh->BeginNode(); n != mesh->EndNode(); ++n){
            auto ff = face_around(ob.m_mesh, V_ind(n->LocalID()));
            std::array<std::array<double, 3>, ntag> res{0};
            std::array<std::array<double, 6>, addtag> resT{0};
            for (auto f: ff) {
                for (int k = 0; k < 3; ++k) {
                    res[0][k] += (EE(f, k) - EEd(f, k)) / EE(f, k);
                    res[1][k] += EE(f, k);
                    res[2][k] += EE(f, k) - EEd(f, k);
                    res[3][k] += EEd(f, k);
                }
                for (int k = 0; k < 6; ++k) {
                    resT[0][k] += TTd(f, k);
                    resT[1][k] += TT(f, k);
                    resT[2][k] += TT(f, k) - TTd(f, k);
                    resT[3][k] += (TT(f, k) - TTd(f, k)) / TT(f, k);
                }
            }
            for (int i = 0; i < 3 && !ff.empty(); ++i) {
                for (int j = 0; j < ntag; ++j)
                    res[j][i] /= ff.size();
            }
            for (int i = 0; i < 6 && !ff.empty(); ++i) {
                for (int j = 0; j < addtag; ++j)
                    resT[j][i] /= ff.size();
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < ntag; ++j)
                    n->RealArray(kk[j])[i] = res[j][i];
            }
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < addtag; ++j) {
                    int r = i;
                    if (i == 4) r = 5;
                    if (i == 5) r = 4;
                    n->RealArray(kk[j + ntag])[i] = resT[j][r];
                    if (i > 2) n->RealArray(kk[j + ntag])[i] *= sqrt(2);
                }
            }
        }
        mesh->SaveVTK(res_path + "/tension_" + postfix + ".vtk");
        delete mesh;
    }
#endif

    t.join();
    return 0;
}

int testSphereStatDef(int argc, char* argv[]){
    double lambda = 1.25; //R / R0
    double R0 = 4, H0 = 0.004;
    double R = lambda * R0;
    double P = 2 / R;
    double mu = P * R0 / (H0 * 2 * (pow(lambda, 6) - 1) / pow(lambda, 7));
    double mesh_h = 0.2;
    auto am = generate_eigth_sphere_part(R, mesh_h);
    Object3D obj{convert_to_Mesh(am, "v:boundary_lbl")};
    obj.name = "sph";
    string gendir = "../generated";
    string res_path = "../result/StatDefSph";
    obj.setDirichletCond([](Object3D& obj, const V_ind& v) -> unsigned char{ return obj.m_boundary[v]; });
    for (auto v: obj.m_mesh.vertices())
        obj.m_x0[v] = CGAL::ORIGIN + (obj.m_x[v] - CGAL::ORIGIN) * R0 / R;

    StaticDefiniteGather::Traits traits;
    StaticDefiniteGather statgath;
    traits.setH(H0).setMu(1e7).setDelta(6e-7).setErr(1e-4);
    traits.setSolveMethod(StaticDefiniteGather::Traits::SIMPLE_NEWTON).setMaxIters(15e5).setFreq(1);
#ifdef USE_INMOST_SOLVERS
    traits.setLinearSolver(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#else
    traits.setLinearSolver(LinearSolver("eigen"));
#endif
    statgath.setTraits(traits).setObject(obj);

    auto aniso_s = obj.m_mesh.add_property_map<F_ind, std::array<DReal , 3>>("f:aniso_s").first;
    for (auto f: obj.m_mesh.faces()){
        auto v = vert_around(obj.m_mesh, f);
        auto t = ((obj.m_x0[v[0]] - CGAL::ORIGIN) + (obj.m_x0[v[1]] - CGAL::ORIGIN) + (obj.m_x0[v[2]] - CGAL::ORIGIN)) / 3;
        double phi = atan2(t[1], t[0]);
        double teta = atan2(hypot(t[0], t[1]), t[2]);
        aniso_s[f] = std::array<DReal , 3>{cos(teta)*cos(phi), cos(teta)*sin(phi), -sin(teta)};
    }
//    double alpha = 0.0245, beta = 0.0297;
//    Force elastic = MooneyRivlinModel(H0, mu, alpha, beta, {NAN, NAN, NAN}, gendir, true,
//                                      "f:aniso_s");
    Force elastic = NeoGookModel(mu, H0, gendir);
    HEDataGatherer gather(elastic.target<HyperElasticForceBase>()->f, true);
//    statgath.m_t.m_nssa._solver.solver.SetParameterReal("relative_tolerance", 1e-14);
//    ResponseStatistics rs, rsd;
    int nfaces = obj.m_mesh.num_faces();
    int niter = 1;
    Eigen::MatrixXd T(9*nfaces, 1), Td(9*nfaces, 1);
    for (int n = 1; n <= niter; ++n){
//        if (false) {
            gather.gatherData(obj, [&T = T.col(n - 1).matrix(), nf = nfaces](const HEDataGatherer::CompDat &dat) {
                for (int i = 0; i < 6; ++i) T[dat.f.idx() + i * nf] = dat.T[i];
                for (int i = 0; i < 3; ++i) T[dat.f.idx() + (i + 6) * nf] = dat.Tee[i];
            });
//        gather.gatherData(w.obj(oid), rs);
//        } else {
            statgath.m_t.setP(P);
            statgath.computeTensionedState();
//        statgath.gatherResponseData(rsd);
        gather.gatherData(statgath.m_w.objs().begin()->second.first, [&T = Td.col(n - 1).matrix(), nf = nfaces](const HEDataGatherer::CompDat &dat) {
            for (int i = 0; i < 6; ++i) T[dat.f.idx() + i * nf] = dat.T[i];
            for (int i = 0; i < 3; ++i) T[dat.f.idx() + (i + 6) * nf] = dat.Tee[i];
        });
            statgath.gatherMyData(
                    [&T = Td, nf = nfaces](const StaticDefiniteGather::ObjsData &dat) -> void {
                        std::pair<int, int> q3[] = {{0, 0},
                                                    {1, 1},
                                                    {2, 2},
                                                    {0, 1},
                                                    {0, 2},
                                                    {1, 2}},
                                q2[] = {{0, 0},
                                        {1, 1},
                                        {1, 0}};
                        if (dat.f.idx() < nf) {
//                        std::cout << "||Td||^2_F = " << dat.T.eval().norm() << std::endl;
                            for (int i = 0; i < 6; ++i) T(dat.f.idx() + i * nf, 0) = dat.T(q3[i].first, q3[i].second);
                            for (int i = 0; i < 3; ++i)
                                T(dat.f.idx() + (i + 6) * nf, 0) = dat.Tee(q2[i].first, q2[i].second);
                            for (int i = 0; i < 9; ++i)
                                if (isnan(T(dat.f.idx() + i * nf, 0)))
                                    std::cout << "Here nan" << std::endl;
                        auto Td = [&T, l = dat.f.idx(), nf](int i) { return T(l + i * nf, 0); };
                        auto sq = [](auto x) { return x * x; };
                        }
                    });
        std::cout   << "||T||_F = "  << T.block(0,0,6*nfaces,1).norm() <<  " ||Td||_F = "  << Td.block(0,0,6*nfaces,1).norm() <<   " ||T - Td||_F = "   << (T - Td).block(0,0,6*nfaces,1).norm()   << std::endl;
        std::cout   << "||Te||_F = "  << T.block(6*nfaces,0,3*nfaces,1).norm() <<  " ||Ted||_F = "  << Td.block(6*nfaces,0,3*nfaces,1).norm() <<   " ||Te - Ted||_F = "   << (T - Td).block(6*nfaces,0,3*nfaces,1).norm()   << std::endl;
//        }
        {
            int nf = nfaces;
            for (auto f: obj.m_mesh.faces()) {
                auto v = vert_around(obj.m_mesh, f);
                auto t = ((obj.m_x[v[0]] - CGAL::ORIGIN) + (obj.m_x[v[1]] - CGAL::ORIGIN) + (obj.m_x[v[2]] - CGAL::ORIGIN)) / 3;
                double phi = atan2(t[1], t[0]);
                double teta = atan2(hypot(t[0], t[1]), t[2]);
                //{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}
                double tt = P * R / 2;
                auto sq = [](auto x) { return x * x; };
                T(f.idx() + 0*nf, n-1) = tt * ( sq(cos(teta)) * sq(cos(phi)) + sq(sin(phi)) );
                T(f.idx() + 1*nf, n-1) = tt * ( sq(cos(teta)) * sq(sin(phi)) + sq(cos(phi)) );
                T(f.idx() + 2*nf, n-1) = tt * ( sq(sin(teta)) );
                T(f.idx() + 3*nf, n-1) = tt * ( -sin(2*phi) * sq(sin(teta)) / 2 );
                T(f.idx() + 4*nf, n-1) = tt * ( -sin(2*teta) * cos(phi) / 2 );
                T(f.idx() + 5*nf, n-1) = tt * ( -sin(2*teta) * sin(phi) / 2 );
//                std::cout << "||T||_F = " << sq(T(f.idx() + 0*nf, n-1)) + sq(T(f.idx() + 1*nf, n-1)) + sq(T(f.idx() + 2*nf, n-1)) +
//                    2 * (sq(T(f.idx() + 3*nf, n-1)) + sq(T(f.idx() + 4*nf, n-1)) + sq(T(f.idx() + 5*nf, n-1))) << std::endl;
            }
        }
        std::cout   << "||T||_F = "  << T.norm() <<  " ||Td||_F = "  << Td.norm() <<   " ||T - Td||_F = "   << (T - Td).norm()   << std::endl;
    }
//    rs.stat.reserve(rsd.stat.size());
//    for (int fi = 0; fi < rsd.stat.size(); ++fi) {
//        auto e = rsd.stat[fi][0];
//        ResponseStatistics::Elem newe = e;
//        auto v = vert_around(obj.m_mesh, F_ind(fi));
//        auto t = ((obj.m_x[v[0]] - CGAL::ORIGIN) + (obj.m_x[v[1]] - CGAL::ORIGIN) + (obj.m_x[v[2]] - CGAL::ORIGIN)) / 3;
//        t /= sqrt(t.squared_length());
//
//        rs.stat.push_back({});
//    }
#ifdef USE_INMOST_SOLVERS
    {
        auto& ob = obj;//w.obj(oid);
        using namespace INMOST;
        INMOST::Mesh* mesh = new INMOST::Mesh();
        vector<HandleType> new_nodes(ob.m_mesh.num_vertices());
        vector<HandleType> new_faces(ob.m_mesh.num_faces());
        for (auto v: ob.m_mesh.vertices()){
            std::array<Storage::real,3> xyz{ob.m_x0[v][0], ob.m_x0[v][1], ob.m_x0[v][2]};
            new_nodes[v.idx()] = mesh->CreateNode(xyz.data())->GetHandle();
        }
        for (auto f: ob.m_mesh.faces()){
            auto v = vert_around(ob.m_mesh, f);
            ElementArray<Node> subarr(mesh);
            subarr.reserve(3);
            for (int j = 0; j < 3; ++j)
                subarr.push_back(new_nodes[v[j].idx()]);
            const int nodesnum[3] = {0,1,2};
            const int sizes[1] = {3};

            new_faces[f.idx()] = mesh->CreateCell(subarr, nodesnum, sizes, 1).first.GetHandle();
        }
        const int ntag = 4;
        const int addtag = 4;
        std::array<Tag, ntag+addtag> kk = {
                mesh->CreateTag("(Eexct-E)/Eexct", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("Eexct", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("Eexct-E", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("E", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("T", DATA_REAL, NODE, FACE, 6),
                mesh->CreateTag("Texct", DATA_REAL, NODE, FACE, 6),
                mesh->CreateTag("Texct-T", DATA_REAL, NODE, FACE, 6),
                mesh->CreateTag("(Texct-T)/Texct", DATA_REAL, NODE, FACE, 6)
        };

        double mx = 0; int ix = -1;
        int n_epoch = 0;
        auto EE = [&T, nfaces, n = n_epoch](F_ind f, int k) { return T(f.idx() + (k + 6)*nfaces, n); };
        auto EEd = [&Td, nfaces, n = n_epoch](F_ind f, int k) { return Td(f.idx() + (k + 6)*nfaces, n); };
        auto TT = [&T, nfaces, n = n_epoch](F_ind f, int k) { return T(f.idx() + (k + 0)*nfaces, n); };
        auto TTd = [&Td, nfaces, n = n_epoch](F_ind f, int k) { return Td(f.idx() + (k + 0)*nfaces, n); };
        for (auto n = mesh->BeginNode(); n != mesh->EndNode(); ++n){
            auto ff = face_around(ob.m_mesh, V_ind(n->LocalID()));
            std::array<std::array<double, 3>, ntag> res{0};
            std::array<std::array<double, 6>, addtag> resT{0};
            for (auto f: ff) {
//                if ((EE(f, 2) - EEd(f, 2)) / EE(f, 2) > mx){
//                    mx = (EE(f, 2) - EEd(f, 2)) / EE(f, 2);
//                    ix = f.idx();
//                }
                for (int k = 0; k < 3; ++k) {
                    res[0][k] += (EE(f, k) - EEd(f, k)) / EE(f, k);
                    res[1][k] += EE(f, k);
                    res[2][k] += EE(f, k) - EEd(f, k);
                    res[3][k] += EEd(f, k);
                }
                for (int k = 0; k < 6; ++k) {
                    resT[0][k] += TTd(f, k);
                    resT[1][k] += TT(f, k);
                    resT[2][k] += TT(f, k) - TTd(f, k);
                    resT[3][k] += (TT(f, k) - TTd(f, k)) / TT(f, k);
                }
            }
            for (int i = 0; i < 3 && !ff.empty(); ++i) {
                for (int j = 0; j < ntag; ++j)
                    res[j][i] /= ff.size();
            }
            for (int i = 0; i < 6 && !ff.empty(); ++i) {
                for (int j = 0; j < addtag; ++j)
                    resT[j][i] /= ff.size();
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < ntag; ++j)
                    n->RealArray(kk[j])[i] = res[j][i];
            }
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < addtag; ++j) {
                    int r = i;
                    if (i == 4) r = 5;
                    if (i == 5) r = 4;
                    n->RealArray(kk[j + ntag])[i] = resT[j][r];
                    if (i > 2) n->RealArray(kk[j + ntag])[i] *= sqrt(2);
                }
            }
        }
//        std::cout << "mx = "<< mx << " ix = " << ix
//                  << " TT(f, 2) =  " << EE(F_ind(ix), 2)  << " TTd(f, 2) = " << EEd(F_ind(ix), 2) << std::endl;
        mesh->SaveVTK(res_path + "/tension.vtk");
        delete mesh;
    }
#endif

    return 0;
}

int generateSyntheticData(int argc, char* argv[]){
#ifndef USE_MATIO
    return 1;
#endif
    std::string path_pref = "../data-bench/membrane_kinematics/";
    std::string res_path = "../result/DataDrivenBench/";

    Eigen::MatrixXd
            X0 = read_mat(path_pref + "x_0.mat", false),
            Y0 = read_mat(path_pref + "y_0.mat", false),
            Z0 = read_mat(path_pref + "z_0.mat", false),
            TRI = read_mat(path_pref + "TRI.mat", false);

    AniMesh am;
    am.vertices.resize(3*X0.rows());
    am.faces.resize(3*TRI.rows());
    for (int i = 0; i < X0.rows(); ++i)
        am.vertices[3*i + 0] = X0(i, 0), am.vertices[3*i + 1] = Y0(i, 0), am.vertices[3*i + 2] = Z0(i, 0);
    for (int i = 0; i < TRI.rows(); ++i)
        for (int k = 0; k < 3; ++k)
            am.faces[3*i + k] = static_cast<int>(TRI(i, k));
    Mesh m = convert_to_Mesh(am);
    remark_boundary(m);
    Object3D obj(m);

    Object3D obj1(res_path + "input_outer.txt");
    vector<int> remap(obj1.m_mesh.num_vertices());
    int new_vert = 0;
    am.vertices.reserve(am.vertices.size() + 3*obj1.m_mesh.num_vertices());
    for (int i = 0; i < obj1.m_mesh.num_vertices(); ++i) {
        int ind = obj1.m_boundary[V_ind(i)] >> 2;
        if (obj1.m_boundary[V_ind(i)] != 1 &&
            !(ind > obj.m_mesh.num_vertices() || (obj1.m_x[V_ind(i)] - obj.m_x[V_ind(ind)]).squared_length() > 1e-5)){
            remap[i] = ind;
        } else {
            remap[i] = obj.m_mesh.num_vertices() + (new_vert++);
            for (int k = 0; k < 3; ++k) am.vertices.push_back(obj1.m_x[V_ind(i)][k]);
        }
    }
//    std::for_each(am.vertices.begin(), am.vertices.end(), [](auto& i) { i /= 1e3; });
    am.face_label.resize(am.faces.size()/3 + obj1.m_mesh.num_faces(), 0);
    std::fill(am.face_label.begin(), am.face_label.begin() + am.faces.size()/3, 1);
    am.faces.reserve(am.faces.size() + 3*obj1.m_mesh.num_faces());
    for (auto f: obj1.m_mesh.faces()){
        auto vv = vert_around(obj1.m_mesh, f);
        for (int k = 0; k < 3; ++k)
            am.faces.push_back(remap[vv[k].idx()] + 1);
    }
    m = convert_to_Mesh(am, "", "f:watch");
    remark_boundary(m);
    Object3D obj_f(m);
//    auto is_movable = [] (Object3D& obj, const V_ind& v) -> bool { return !(obj.m_boundary[v] & 1); };
//    obj_f.set_is_movable_checker(is_movable);

    double H = 1, R = 15, R0 = 7, dP = 3;
    double mu_1 = 0.0122*1e3, mu_2 = 0.2747*1e3, gamma = 1.6090, kappa = 0.4019, teta = 1.5560;
    std::array<double, 3> Av = {cos(teta), sin(teta), 0};
    HGOModel1 hmf(H, mu_1, mu_2, gamma, kappa, Av, "../test_gen", false);
    hmf.prepareJacobianFunction(false);
    double P = 0;
    SimplePressureLoad Pr([&P](Object3D&, F_ind) { return P; });

//    thread t(start_debug_gui, argc, argv);
    World w;
//    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    int nfaces = obj_f.m_mesh.num_faces();
    auto oid = w.addObject3D(std::move(obj_f), 1);
    auto efid = w.addForce(std::move(hmf), oid);
    auto pid = w.addForce(std::move(Pr), oid);
    NewtonSolverStepAlgo nssa;
#ifdef USE_INMOST_SOLVERS
    nssa._solver.solver = std::move(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#else
    nssa._solver.solver = std::move(LinearSolver("eigen"));
#endif
    nssa._solver.solver.setVerbosityLevel(0);
    nssa.newton_algo = [](NewtonSolverStepAlgo& ns, World* w){
        ns.updateMeshData(w);
        int N = ns.getNodfs(w);
        ns.assembleSystem(w, N);
        double tau = 0.1;//*1e-4;//0.3;
        static double resid_init = -1;
        if (ns.m_it == 0) {
//            ns._solver.solver.SetParameterReal("absolute_tolerance", 1e-15);
//            ns._solver.solver.SetParameterReal("relative_tolerance", 1e-20);
//            double droptol = 2e-6;// /5;
            resid_init = ns.computeRhsNorm(1);
//            std::cout << "\tresid = " << resid_init
//                << " drop_tol = " << droptol
//                << std::endl;
//            ns._solver.solver.SetParameterReal("drop_tolerance", droptol);
//            ns._solver.solver.SetParameterReal("reuse_tolerance", 0.1 * droptol);
        } else {
            double eps = ns.computeRhsNorm(1) / resid_init;
//            double droptol = 1e-4/2;// /4;//8e-4;//4e-5;//1e-6;
            if (eps > 100){
                tau = 1e-1;
//                droptol = 5e-5/2;// /4;//1e-4;
            } else if (eps < 0.06){
                tau = 1;//*1e-3;
            } else {
                tau = 0.3;// * 1e-4;
            }
//            std::cout << "\teps = " << eps
//                << " drop_tol = " << droptol
//                << std::endl;
//            ns._solver.solver.SetParameterReal("drop_tolerance", droptol);
//            ns._solver.solver.SetParameterReal("reuse_tolerance", 0.1 * droptol);
        }
        ns.Solve();
        ns.applySolution(w, tau);
        ns.m_it++;
        return 0;
    };
    int freq = 1000, maxits = 1e6;
    double delta = 2e-3;//2e-5; // 0/100
    double deltaP = delta;
    w.setForceApplier(StaticForceApplierP(deltaP));
//    w.setStepSimulationAlgo(std::move(nssa));
//    int freq = 1, maxits = 100;

    int status = 0;
    double get_eps = -1;
    auto stopConditionGen = [&status, &eps = get_eps](auto &freq, auto &err, auto &time, auto &maxits) {
        return [&freq, &err, &time, &maxits, &status, &eps](StepSimInfo& info, World *w) -> bool {
            static double resid_init;
            int it = info.it;
            if (it == 1) resid_init = w->getWorldNorm("v:force");
//            if (it > 0) sleep(10);
            if (it > 0 && it % freq == 0) {
                double resid = w->getWorldNorm("v:force");
                eps = resid / resid_init;
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed()
                          << "\n";
                if (eps < err) {
                    std::cout << "Algorithm is converged: \n";
                    status = 0;
                    return true;
                }
                if (it > maxits || std::isnan(eps)) {
                    std::cout << "Algorithm is diverged or reached maximum of time - iters: \n";
                    status = 1;
                    return true;
                }
            }
            return false;
        };
    };
    double err = 1e-4;
    auto time = World3d::Timer();
    auto stopCondition = stopConditionGen(freq, err, time, maxits);

    HEDataGatherer gather(w.obj(oid).m_forces[efid].target<HyperElasticForceBase>()->f, true);

    StaticDefiniteGather::Traits traits;
    StaticDefiniteGather statgath;
    traits.setH(H).setMu(1e10).setDelta(6e-12).setErr(1e-4);
    traits.setSolveMethod(StaticDefiniteGather::Traits::SIMPLE_NEWTON).setMaxIters(15).setFreq(1);
#ifdef USE_INMOST_SOLVERS
    traits.setLinearSolver(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#else
    traits.setLinearSolver(LinearSolver("eigen"));
#endif
    statgath.setTraits(traits).setObject(w.obj(oid));
//    statgath.m_t.m_nssa._solver.solver.SetParameterReal("relative_tolerance", 1e-14);
    ResponseStatistics rs, rsd;

    Eigen::MatrixXd T(9*nfaces, 29), Td(9*nfaces, 29);
    for (int n = 1; n <= 23; ++n){
        P = n * dP;
        if (n  >= 24) {
            deltaP = delta / (4 * n);
            std::cout << "\nPressure = " << P << std::endl;
            if (n >= 2) {
                w.setStepSimulationAlgo([&ssa = nssa](World *w) { return ssa(w); });
                freq = 1, maxits = 10;
            }
            time.reset();
            w.Simulation(stopCondition);
            if (n > 1 && status != 0) {
                w.clearStepAlgo();
                freq = 1000, maxits = 5e4;
                double tmp_e = err;
                err = err / get_eps;
                w.Simulation(stopCondition);
                err = tmp_e;
            }
            w.obj(oid).save(res_path + "synth/synth_" + to_string(n) + ".txt");
        } else {
            Object3D tmp(res_path + "synth/synth_" + to_string(n) + ".txt");
            for (auto v: tmp.m_mesh.vertices())
                w.obj(oid).m_x[v] = tmp.m_x[v];

        }
        gather.gatherData(w.obj(oid), [&T = T.col(n-1).matrix(), nf = nfaces](const HEDataGatherer::CompDat& dat){
            for (int i = 0; i < 6; ++i) T[dat.f.idx() + i*nf] = dat.T[i];
            for (int i = 0; i < 3; ++i) T[dat.f.idx() + (i+6)*nf] = dat.Tee[i];
        });
        gather.gatherData(w.obj(oid), rs);

        statgath.m_t.setP(P);
        statgath.computeTensionedState();
        statgath.gatherResponseData(rsd).gatherMyData(
                [&T = Td.col(n-1).matrix(), nf = nfaces](const StaticDefiniteGather::ObjsData& dat)->void {
            std::pair<int, int> q3[] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}},
                                q2[] = {{0, 0}, {1, 1}, {1, 0}};
            if (dat.f.idx() < nf) {
                for (int i = 0; i < 6; ++i) T[dat.f.idx() + i * nf] = dat.T(q3[i].first, q3[i].second);
                for (int i = 0; i < 3; ++i) T[dat.f.idx() + (i + 6) * nf] = dat.Tee(q2[i].first, q2[i].second);
            }
        });
    }
    for_each(rs.stat.begin(), rs.stat.end(), [](auto& i){ i.resize(29); });
    for_each(rsd.stat.begin(), rsd.stat.end(), [](auto& i){ i.resize(29); });
    rs.save_to_file(res_path + "synth/gather_data_synth.tss");
    rsd.save_to_file(res_path + "synth/gather_data_synth_eval.tss");

    auto print = [](std::ostream& out, auto arr) -> std::ostream&{
        if (arr.empty()) return out;
        out << arr[0];
        for (int i = 1; i < arr.size(); ++i)
            out << ", " << arr[i];
        return out;
    };

    auto lbl = w.obj(oid).m_mesh.property_map<F_ind, int>("f:watch").first;
    ofstream f(res_path + "synth/synth_data.dist");
    for (int i = 0; i < nfaces; ++i){
//        if (lbl[F_ind(i)] != 1) continue;
        for (int n = 0; n < 29; ++n){
            print(f, rs.stat[i][n].ksi) << ", ";
            print(f, rs.stat[i][n].response) << ", ";
            for (int k = 0; k < 9; ++k) f << T(i + k*nfaces, n) << ", ";
            print(f, rsd.stat[i][n].ksi) << ", ";
            print(f, rsd.stat[i][n].response) << ", ";
            for (int k = 0; k < 9; ++k) f << Td(i + k*nfaces, n) << ", ";
            f << lbl[F_ind(i)]  << "\n";
        }
    }
    f.close();

#ifdef USE_INMOST_SOLVERS
    {
        auto& ob = obj;//w.obj(oid);
        using namespace INMOST;
        INMOST::Mesh* mesh = new INMOST::Mesh();
        vector<HandleType> new_nodes(ob.m_mesh.num_vertices());
        vector<HandleType> new_faces(ob.m_mesh.num_faces());
        for (auto v: ob.m_mesh.vertices()){
            std::array<Storage::real,3> xyz{ob.m_x0[v][0], ob.m_x0[v][1], ob.m_x0[v][2]};
            new_nodes[v.idx()] = mesh->CreateNode(xyz.data())->GetHandle();
        }
        for (auto f: ob.m_mesh.faces()){
            auto v = vert_around(ob.m_mesh, f);
            ElementArray<Node> subarr(mesh);
            subarr.reserve(3);
            for (int j = 0; j < 3; ++j)
                subarr.push_back(new_nodes[v[j].idx()]);
            const int nodesnum[3] = {0,1,2};
            const int sizes[1] = {3};

            new_faces[f.idx()] = mesh->CreateCell(subarr, nodesnum, sizes, 1).first.GetHandle();
        }
        const int ntag = 4;
        std::array<Tag, ntag+2> kk = {
                mesh->CreateTag("(Eexct-E)/Eexct", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("Eexct", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("Eexct-E", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("E", DATA_REAL, NODE, FACE, 3),
                mesh->CreateTag("T", DATA_REAL, NODE, FACE, 6),
                mesh->CreateTag("Texct", DATA_REAL, NODE, FACE, 6)
        };

        double mx = 0; int ix = -1;
        int n_epoch = 0;
        auto EE = [&T, nfaces, n = n_epoch](F_ind f, int k) { return T(f.idx() + (k + 6)*nfaces, n); };
        auto EEd = [&Td, nfaces, n = n_epoch](F_ind f, int k) { return Td(f.idx() + (k + 6)*nfaces, n); };
        auto TT = [&T, nfaces, n = n_epoch](F_ind f, int k) { return T(f.idx() + (k + 0)*nfaces, n); };
        auto TTd = [&Td, nfaces, n = n_epoch](F_ind f, int k) { return Td(f.idx() + (k + 0)*nfaces, n); };
        for (auto n = mesh->BeginNode(); n != mesh->EndNode(); ++n){
            auto ff = face_around(ob.m_mesh, V_ind(n->LocalID()));
            std::array<std::array<double, 3>, ntag> res{0};
            std::array<std::array<double, 6>, 2> resT{0};
            for (auto f: ff) {
                if ((EE(f, 2) - EEd(f, 2)) / EE(f, 2) > mx){
                    mx = (EE(f, 2) - EEd(f, 2)) / EE(f, 2);
                    ix = f.idx();
                }
                for (int k = 0; k < 3; ++k) {
                    res[0][k] += (EE(f, k) - EEd(f, k)) / EE(f, k);
                    res[1][k] += EE(f, k);
                    res[2][k] += EE(f, k) - EEd(f, k);
                    res[3][k] += EEd(f, k);
                }
                for (int k = 0; k < 6; ++k) {
                    resT[0][k] += TTd(f, k);
                    resT[1][k] += TT(f, k);
                }
            }
            for (int i = 0; i < 3 && !ff.empty(); ++i) {
                for (int j = 0; j < ntag; ++j)
                    res[j][i] /= ff.size();
            }
            for (int i = 0; i < 6 && !ff.empty(); ++i) {
                for (int j = 0; j < 2; ++j)
                    resT[j][i] /= ff.size();
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < ntag; ++j)
                    n->RealArray(kk[j])[i] = res[j][i];
            }
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 2; ++j) {
                    int r = i;
                    if (i == 4) r = 5;
                    if (i == 5) r = 4;
                    n->RealArray(kk[j + ntag])[i] = resT[j][r];
                }
            }
        }
        std::cout << "mx = "<< mx << " ix = " << ix
        << " TT(f, 2) =  " << EE(F_ind(ix), 2)  << " TTd(f, 2) = " << EEd(F_ind(ix), 2) << std::endl;
        mesh->SaveVTK(res_path + "/synth/tension.vtk");
        delete mesh;
    }
#endif

//    t.join();

    return 0;
}

int computeData(int argc, char* argv[]){
#ifndef USE_MATIO
    return 1;
#endif
    std::string path_pref = "../data-bench/membrane_kinematics/";
    std::string res_path = "../result/DataDrivenBench/";

    Eigen::MatrixXd
            X0 = read_mat(path_pref + "x_0.mat", false),
            Y0 = read_mat(path_pref + "y_0.mat", false),
            Z0 = read_mat(path_pref + "z_0.mat", false),
            Xn = read_mat(path_pref + "x_n.mat", false),
            Yn = read_mat(path_pref + "y_n.mat", false),
            Zn = read_mat(path_pref + "z_n.mat", false),
            TRI = read_mat(path_pref + "TRI.mat", false),
            F_d = read_mat(path_pref + "DEF_GRADIENT.mat", false),
            Sigma_d = read_mat(path_pref + "SIGMA.mat", false);

    AniMesh am;
    am.vertices.resize(3*X0.rows());
    am.faces.resize(3*TRI.rows());
    for (int i = 0; i < X0.rows(); ++i)
        am.vertices[3*i + 0] = X0(i, 0), am.vertices[3*i + 1] = Y0(i, 0), am.vertices[3*i + 2] = Z0(i, 0);
    for (int i = 0; i < TRI.rows(); ++i)
        for (int k = 0; k < 3; ++k)
            am.faces[3*i + k] = static_cast<int>(TRI(i, k));
    Mesh m = convert_to_Mesh(am);
    remark_boundary(m);
    Object3D obj(m);
    obj.save(res_path + "input.stl");
    obj.save(res_path + "input.txt");

    double H = 2.34, R = 15, R0 = 7, dP = 3;

    Object3D obj1(res_path + "input_outer.txt");
    vector<int> remap(obj1.m_mesh.num_vertices());
    int new_vert = 0;
    am.vertices.reserve(am.vertices.size() + 3*obj1.m_mesh.num_vertices());
    for (int i = 0; i < obj1.m_mesh.num_vertices(); ++i) {
        int ind = obj1.m_boundary[V_ind(i)] >> 2;
        if (obj1.m_boundary[V_ind(i)] != 1 &&
            !(ind > obj.m_mesh.num_vertices() || (obj1.m_x[V_ind(i)] - obj.m_x[V_ind(ind)]).squared_length() > 1e-5)){
            remap[i] = ind;
        } else {
            remap[i] = obj.m_mesh.num_vertices() + (new_vert++);
            for (int k = 0; k < 3; ++k) am.vertices.push_back(obj1.m_x[V_ind(i)][k]);
        }
    }
    am.faces.reserve(am.faces.size() + 3*obj1.m_mesh.num_faces());
    for (auto f: obj1.m_mesh.faces()){
        auto vv = vert_around(obj1.m_mesh, f);
        for (int k = 0; k < 3; ++k)
            am.faces.push_back(remap[vv[k].idx()] + 1);
    }
    m = convert_to_Mesh(am);
    remark_boundary(m);
    Object3D obj_f(m);
    obj_f.save(res_path + "input_full.stl");
    obj_f.save(res_path + "input_full.txt");

    ResponseStatistics rs;
    Eigen::MatrixXd F, Sigma;
    F.resize(F_d.rows(), F_d.cols()), Sigma.resize(Sigma_d.rows(), Sigma_d.cols());
    StaticDefiniteGather::Traits traits;
    StaticDefiniteGather statgath;
    traits.setH(H).setMu(1e10).setDelta(3e-12).setErr(1e-5);
    traits.setSolveMethod(StaticDefiniteGather::Traits::SIMPLE_NEWTON).setMaxIters(100).setFreq(1);
#ifdef USE_INMOST_SOLVERS
    traits.setLinearSolver(LinearSolver(INMOST::Solver::INNER_MPTILUC));
#else
    traits.setLinearSolver(LinearSolver("eigen"));
#endif
    traits.m_nssa._solver.solver.setVerbosityLevel(0);
    for (int n = 0; n < Xn.cols(); ++n){
        traits.setP(dP*(n + 1));
        for (auto v: obj.m_mesh.vertices())
            obj.m_next_x[v] = obj.m_x[v] = Point(Xn(v.idx(), n), Yn(v.idx(), n), Zn(v.idx(), n));
        statgath.setTraits(traits).setObject(obj);
        statgath.computeTensionedState().gatherResponseData(rs);
        statgath.gatherTension(Sigma.col(n).matrix());
//        Eigen::MatrixXd Watch(Sigma.rows(), 2);
//        Watch.col(0) = Sigma.col(n);
//        Watch.col(1) = Sigma_d.col(n);
//        for (int i = 0; i < obj.m_mesh.num_faces(); ++i) {
//            auto v  = vert_around(obj.m_mesh, F_ind(i));
//            if (!obj.is_movable(v[0]) || !obj.is_movable(v[1]) || !obj.is_movable(v[2]))
//                continue;
//            std::cout << Watch.row(i) << "\n";
//            std::cout << Watch.row(i+1*obj.m_mesh.num_faces()) << "\n";
//            std::cout << Watch.row(i+2*obj.m_mesh.num_faces()) << "\n" << std::endl;
//        }
        statgath.m_w.objs().begin()->second.first.save(res_path + "deform_" + to_string(n+1) + ".stl");
        gatherDefGrad(obj, F.col(n).matrix());
    }
    std::cout << "||F - F_d||_frob = " << (F - F_d).norm() << std::endl;
    std::cout << "||T - T_d||_frob = " << (Sigma - Sigma_d).norm() << std::endl;
    rs.save_to_file(res_path + "gather_data.tss");
    rs.save_ksi_distribution(res_path + "gather_xyz.dist");
    rs.save_distribution(res_path + "gather.dist");

    ofstream f(res_path + "tension.dist");
    for (int i = 0; i < obj.m_mesh.num_faces(); ++i){
        for (int n = 0; n < Xn.cols(); ++n) {
            int nf = obj.m_mesh.num_faces();
            f << Sigma(i + 0*nf, n);
            for (int k = 1; k < 3; ++k)
                f << ", " << Sigma(i + k*nf, n);
            f << "\n";
        }
    }
    f.close();

    return 0;
}

int BenchDataDriven(int argc, char* argv[]){
//    return testCircleStatDef(argc, argv);
    return testSphereStatDef(argc, argv);
//    sleep(1);
//    return testDataDriven(argc, argv);
//    testSigmaStability(argc, argv);
//    return processGeneratedState(argc, argv);
    std::string path_pref = "../data-bench/membrane_kinematics/";
    std::string res_path = "../result/DataDrivenBench/";
//    return computeData(argc, argv);
    return generateSyntheticData(argc, argv);

    double H = 2.34, R = 15, R0 = 7, dP = 3;
    ResponseStatistics rs(res_path + "gather_data.tss");

    thread t(start_debug_gui, argc, argv);
    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());

    Object3D obj(res_path + "input_full.txt");
    Object3D obj0(res_path + "input.txt");
    if (false){
        auto bcl = obj0.m_mesh.add_property_map<F_ind, int>("f:bnd_close").first;
        std::fill(bcl.begin(), bcl.end(), 0);
        for (auto f: obj0.m_mesh.faces()){
            auto v = vert_around(obj0.m_mesh, f);
            if (!obj0.is_movable(v[0]) || !obj0.is_movable(v[1]) || !obj0.is_movable(v[2]))
                bcl[f] = 1;
        }

        bool changed = false;
        int nits = 2000, nit = 1;
        std::vector<int> tmp_bcl(obj0.m_mesh.num_faces());
        do {
            std::copy(bcl.begin(), bcl.end(), tmp_bcl.begin());
            changed = false;
            ++nit;
            for (auto f: obj0.m_mesh.faces()) {
                if (bcl[f] != 0) continue;
                auto e = edge_around(obj0.m_mesh, f);
                std::array<F_ind, 3> ff;
                for (int i = 0; i < 3; ++i) {
                    auto fa = face_around(obj0.m_mesh, e[i]);
                    ff[i] = (fa[1].first != f) ? fa[1].first : fa[0].first;
                }
                auto q = {bcl[ff[0]], bcl[ff[1]], bcl[ff[2]]};
                int lvl = std::max(q);
                if (lvl > 0) {
                    tmp_bcl[f] = lvl + 1;
                    changed = true;
                }
            }
            std::copy(tmp_bcl.begin(), tmp_bcl.end(), bcl.begin());
        } while (changed && nit < nits);
        std::cout << "nit = " << --nit << std::endl;
        std::cout << "max_lvl = " << *std::max_element(bcl.begin(), bcl.end()) << std::endl;

        ofstream f(res_path + "gather.dist");
        auto print = [](std::ostream& out, auto arr) -> std::ostream&{
            if (arr.empty()) return out;
            out << arr[0];
            for (int i = 1; i < arr.size(); ++i)
                out << ", " << arr[i];
            return out;
        };
        for (int i = 0; i < rs.stat.size(); ++i){
            for (int j = 0; j < rs.stat[i].size(); ++j) {
                print(print(f, rs.stat[i][j].ksi) << ", ", rs.stat[i][j].response) << ", " << bcl[F_ind(i)];

                f << "\n";
            }
        }
        f.close();

        return 0;
    }

    int from_step = 2, k_step = 4;
    LinearZeroSearcher::Cloud stepDat; stepDat.reserve(rs.stat.size()*k_step);
    for (int i = 0; i < rs.stat.size(); ++i) {
        auto v = vert_around(obj0.m_mesh, F_ind(i));
        if (!obj0.is_movable(v[0]) || !obj0.is_movable(v[1]) || !obj0.is_movable(v[2]))
            continue;
        for (int j = from_step; j < from_step + k_step; ++j)
            stepDat.push_back(rs.stat[i][j]);
    }
    //filter data
    {
        std::vector<pair<int, int>> ids;
        ids.reserve(rs.stat.size() * rs.stat[0].size());
        for (int i = 0; i < rs.stat.size(); ++i)
            for (int j = 0; j < rs.stat[i].size(); ++j)
                ids.emplace_back(i, j);
        auto getDat = [&rs](pair<int, int> id) { return rs.stat[id.first][id.second]; };
        auto sq = [](auto i) { return i*i; };
        double median[3] = {0}, sd[3] = {0};
        double scl = 1 / 1.2;
        for (int k = 0; k < 3; ++k){
            auto median_it = ids.begin() + (ids.size() - 1) / 2;
            std::nth_element(ids.begin(), median_it , ids.end(), [k, &getDat](const pair<int, int>& id1, const pair<int, int>& id2) {
                return getDat(id1).response[k] < getDat(id2).response[k];
            });
            median[k] = getDat(*median_it).response[k];
            double mean = 0;
            for (auto& i: rs.stat) for (auto& j: i) mean += j.response[k];
            mean /= ids.size();
            for (auto& i: rs.stat) for (auto& j: i) sd[k] += sq(j.response[k] - mean);
            sd[k] /= (ids.size() - 1);
            sd[k] = sqrt(sd[k]);
        }
        for (auto& i: rs.stat)
            i.erase(std::remove_if(i.begin(), i.end(),
                                   [median, sd, scl](const auto& a) {
                bool rm = false;
                for (int k = 0; k < 3; ++k)
                    rm = rm || (fabs(a.response[k] - median[k]) > sd[k] * scl);
                return rm;}), i.end());
    }
    rs.combineAll();
    LinearZeroSearcher lzs;
    lzs.setNearZeroPointChooser([&stepDat](const LinearZeroSearcher::Cloud& cloud)-> LinearZeroSearcher::Cloud {
        LinearZeroSearcher::Cloud res;
//        double nzpc_c = 0.03;
//        auto mes = [](std::array<double, 3> p) { return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];};
//        for (const auto& i: cloud)
//            if (mes(i.ksi) < nzpc_c * nzpc_c)
//                res.push_back(i);
        res = stepDat;
        return res;
    });

    auto f = makeInterpolantFactory3D(KNearestSearcher3D(1), lzs);
    int nlzs = 0, nknn = 0;
    f.setStrategy([&nlzs, &nknn](decltype(f)& fact, std::array<double, 3> x) -> std::array<double, 3>{
        std::array<double, 3> res{0};
        //a is for too close to zero case and b is for too far to point cloud case
        double a = 0.001;
        double b = 0.1;
        KNearestSearcher3D& kn = *static_cast<KNearestSearcher3D*>(fact.factory[0].get());
        LinearZeroSearcher& lzs = *static_cast<LinearZeroSearcher*>(fact.factory[1].get());
        double l = fabs(x[0]) + fabs(x[1]) + fabs(x[2]);
        if (x[0] < 0 || x[1] < 0 || l < a){
            res = lzs(x); nlzs++;
            fact.m_interpolated = lzs.interpolated();
        } else
        {
            res = kn(x); nknn++;
            fact.m_interpolated = kn.interpolated();
            auto x0 = kn.lastX();
            double dl = std::max({fabs(x[0] - x0[0]), fabs(x[1] - x0[1]), fabs(x[2] - x0[2])});
            if (dl > b){
                res = lzs(x); nlzs++; nknn--;
                fact.m_interpolated = lzs.interpolated();
            }
        }

        return res;
    });
    World3d::DataDrivenElastic force(rs.stat[0], f);
    force.init();

    int freq = 100, maxits = 1e5; double err = 1e-4;
    double delta = 6e-6;
    double P = 10;
    auto id = w.addObject3D(move(obj), 1);
    auto pid = w.addForce(SimplePressureLoad(P), id);
    w.addForce(force);//NeoGookModelOpt(1e4, H));
    w.setForceApplier(StaticForceApplier(delta), id);
    for (int n = 10; n <= 29; ++n) {
        w.removeForce(id, pid);
        pid = w.addForce(SimplePressureLoad(P = dP*n), id);
        auto time = World3d::Timer();
        w.Simulation([&freq, &err, &time, &maxits, &nlzs, &nknn](StepSimInfo& info, World *w) -> bool {
            static double resid_init;
            int it = info.it;
            if (it == 1) resid_init = w->getWorldNorm("v:force");
            if (it > 0 && it % freq == 0) {
                double resid = w->getWorldNorm("v:force");
                double eps = resid / resid_init;
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed()
                          << " : nlzs / n_interp = " << (nlzs * 100.0) / (nknn + nlzs) << "%\n";
                nlzs = nknn = 0;
                if (eps < err) {
                    std::cout << "Algorithm is converged: \n";
                    return true;
                }
                if (it > maxits || std::isnan(eps)) {
                    std::cout << "Algorithm is diverged or reached maximum of time - iters: \n";
                    return true;
                }
            }
            return false;
        });
        w.objs().begin()->second.first.save(res_path + "reestablish_" + to_string(n) + ".stl");
        w.objs().begin()->second.first.save(res_path + "reestablish_" + to_string(n) + ".txt");
        break;
    }
    t.join();

    return 0;
}


