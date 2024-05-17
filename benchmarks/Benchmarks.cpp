//
// Created by alex on 21.09.2020.
//

#include "Benchmarks.h"
#include "NonLinearSolverCustom.h"

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

void computeShiftOperator(){
    using SX = casadi::SX;
    using Function = casadi::Function;
    SX P = SX::sym("P", 3, 3);
    SX X = SX::sym("X", 3, 3);
    SX Q = P + X;
    SX H = SX::sym("H");
    SX mu = SX::sym("mu");
    auto Pcol = SX::horzsplit(P, 1);
    auto Qcol = SX::horzsplit(Q, 1);
    SX Sp = 0.5 * SX::cross(Pcol[1] - Pcol[0], Pcol[2] - Pcol[0]);
    SX Sq = 0.5 * SX::cross(Qcol[1] - Qcol[0], Qcol[2] - Qcol[0]);
    SX Ap = SX::norm_2(Sp), Aq = SX::norm_2(Sq);
    SX np = Sp / Ap, nq = Sq / Aq;
    SX _Dp[3], _Dq[3];
    for (int i = 0; i < 3; ++i) {
        _Dp[i] = SX::cross(np, Pcol[(i + 2) % 3] - Pcol[(i + 1) % 3]);
        _Dq[i] = SX::cross(nq, Qcol[(i + 2) % 3] - Qcol[(i + 1) % 3]);
    }
    SX Dp = SX::horzcat({_Dp[0], _Dp[1], _Dp[2]}) / (2 * Ap);
    SX Dq = SX::horzcat({_Dq[0], _Dq[1], _Dq[2]}) / (2 * Ap);
    SX F = -Ap * H * mu * (Q * Dp.T() * Dp - SX::pow(Ap / Aq, 3) * Dq);
    SX Fr = SX::reshape(F, F.size1() * F.size2(), 1);
    SX Xr = SX::reshape(X, X.size1()*X.size2(), 1);
    SX dFdX = SX::jacobian(Fr, Xr);
    dFdX = SX::simplify(SX::substitute(dFdX, X, SX::zeros(3,3)));
    Function dF("dF", {P, H, mu}, {dFdX});
    double h = 0.25;
    vector<double> Pin = {0, 0, 0,    h, 0, 0,    h/2, sqrt(3.0)/2*h, 0};
    vector<double> Res(9*9);
    double Hin = 0.1;
    double muin = 1.2e6/3;
    const double* in[] = {Pin.data(), &Hin, &muin};
    double* out[] = {Res.data()};
    vector<casadi::DM> iin(3), oout(1);
    iin[0] = casadi::DM(3,3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            iin[0](i,j) = Pin[i + 3 * j];
    iin[1] = Hin;
    iin[2] = muin;

    dF.call(iin, oout);
//    std::cout << oout << std::endl;
    Eigen::Matrix<double, 9, 9> M, I;
    I.setIdentity();
    for (int q = 0, j = 0; j < 9; ++j)
        for (int i = 0; i < 9; ++i)
            M(i, j) = oout[0].has_nz(i, j) ? oout[0].get_nonzeros()[q++] : 0.0;
    std::cout << M <<std::endl;
    double Mdelta = 10e-6;
    for (int i = 1, cnt = 100; i < cnt; ++i){
        double delta = Mdelta / cnt * i;
//        Eigen::Matrix<double, 9, 9> R = /*I + */delta * M;
//        Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(R);
//        cout << "delta = " <<  delta << ": ||I + d dF/dX|| = " << svd.singularValues()[0] << endl;
        Eigen::Matrix<double, 9, 9> R = /*I + */delta * M;
//        Eigen::Matrix<double, 9, 9> RR = R.transpose() * R;
        Eigen::EigenSolver<Eigen::Matrix<double, 9, 9>> es(R);
        cout << "delta = " <<  delta << ": lambda(I + d dF/dX) = \n" << es.eigenvalues() << endl;

    }
    exit(0);
}

int BenchKyriacou(int argc, char **argv){
    thread t(start_debug_gui, argc, argv);

    double p = 1;
    double delta = 1.6e-2;//1.8e-6//3.8/2
    double err = 1e-9;
    double mesh_h = 0.025;//0.2495;
    string to_save = "../result";
    string prefix = "";
    int scenario = 10;

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << std::endl;
            std::cout << " Command line options: " << std::endl;
            std::cout << "  -f, --force    <Force value>" << std::endl;
            std::cout << "  -d, --delta    <Solver parameter>" << std::endl;
            std::cout << "  -e, --err      <Parameter of stop criteria>" << std::endl;
            std::cout << "  -h, --mesh_h   <Step of generated mesh>" << std::endl;
            std::cout << "  -t, --target   <Directory to save results>" << std::endl;
            std::cout << "  -p, --prefix   <Prefix to file name>" << std::endl;
            std::cout << "  -s, --scenario <Choosing of the running scenario: 0 - SVK, 1 - MR, 2 - MY, 3 - GOH>" << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            p = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--delta") == 0) {
            delta = atof(argv[i + 1]);
            i++;
            continue;
        }

        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--err") == 0) {
            err = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--mesh_h") == 0) {
            mesh_h = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0) {
            to_save = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) {
            prefix = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--scenario") == 0) {
            scenario = atoi(argv[i + 1]);
            i++;
            continue;
        }
    }
    std::cout << "Input parameters: " <<
              "\n\tforce = " << p <<
              "\n\tdelta = " << delta <<
              "\n\terr = " << err <<
              "\n\tmesh_h = " << mesh_h <<
              "\n\tscenario = " << scenario <<
              "\n\tdirectory = " << to_save << std::endl;
    {
        path path = to_save;
        auto s = status(path);
        if (!status_known(s) || s.type() != file_type::directory)
            create_directory(path);
    }

    string fname = to_save + "/" + prefix + "benchmarkKyriacou_" + to_string(p) + "_log.txt";
    //std::ofstream f(fname);
    auto& f = std::cout;

    f << "Input parameters: " <<
      "\n\tforce = " << p <<
      "\n\tdelta = " << delta <<
      "\n\terr = " << err <<
      "\n\tmesh_h = " << mesh_h <<
      "\n\tscenario = " << scenario <<
      "\n\tdirectory = " << to_save << std::endl;

    double lambda = 1.1;
    double L = 1.0;
    double a = L, b = L;
    static const int LEFT = 1 << 0, RIGHT = 1 << 1, BOTTOM = 1 << 2, TOP = 1 << 3, OXALIGN = 1 << 4, OYALIGN = 1 << 5;
    AniMesh am = generate_dual_aligned_rectangle(a, b, mesh_h);

    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = "rectangle";
    obj.set_is_movable_checker([] (Object3D& obj, const V_ind& v) -> bool
                               { return !(obj.m_boundary[v] & (LEFT | RIGHT | BOTTOM | TOP)); });
    for (auto v: obj.m_mesh.vertices()) {
        obj.m_x[v] = Point(lambda * obj.m_x[v][0], lambda * obj.m_x[v][1], 0);
        obj.m_next_x[v] = obj.m_x[v];
    }
    Force Pr;
    string dir = "../generated";
    Force elastic, bending;
    if (scenario == 1){
        Pr = SimplePressureLoad(p);
        double H = 1;
        double mu = 1, alpha = 0.0, beta = 1.0, teta = M_PI_4;
        std::array<double, 3> A = {cos(teta), sin(teta), 0};
        elastic = MooneyRivlinModel(H, mu, alpha, beta, A, dir, true);
        elastic.target<MooneyRivlinModel>()->prepareJacobianFunction();
//        bending = BendingForce(elastic.target<MooneyRivlinModel>()->f, true);
    } else if (scenario == 10){
        Pr = SimplePressureLoad(p);
        double H = 0.5;
        double mu = 2, alpha = 0.0, beta = 1.0, teta = M_PI_4;
        std::array<double, 3> A = {cos(teta), sin(teta), 0};
        elastic = MooneyRivlinModel(H, mu, alpha, beta, A, dir, true);
        elastic.target<HyperElasticForceBase>()->prepareJacobianFunction();
        bending = BendingForce(elastic.target<HyperElasticForceBase>()->f, true);
        bending.target<BendingForce>()->prepareJacobianFunction();
    }
        else throw std::runtime_error("Is not implemented");

    map<E_ind, BendingForce::BoundaryData::Entry> boundary;
    std::array<vector<E_ind>, 2> align;
    for (auto e: obj.m_mesh.edges()){
        auto f = face_around(obj.m_mesh, e);
        auto v = vert_around(obj.m_mesh, e);
        int amask[2] = {OXALIGN, OYALIGN};
        for (int i = 0; i < 2; ++i) {
            if ((obj.m_boundary[v[0]] & amask[i]) && (obj.m_boundary[v[1]] & amask[i])) {
                align[i].push_back(e);
                break;
            }
        }
        if (f[0].second && f[1].second)
            continue;
        int mask[] = {LEFT, RIGHT, BOTTOM, TOP};
        DReal vec[4][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}};//{{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}};
        for (int i = 0; i < 4; ++i) {
            if ((obj.m_boundary[v[0]] & mask[i]) && (obj.m_boundary[v[1]] & mask[i])) {
                boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, &vec[i][0]);
                break;
            }
        }
    }

    V_ind vc;
    for (auto v: obj.m_mesh.vertices()){
        if ((obj.m_boundary[v] & (OXALIGN | OYALIGN)) == (OXALIGN | OYALIGN)){
            vc = v;
            break;
        }
    }

    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto id = w.addObject3D(move(obj), 1);
    bending.target<BendingForce>()->set_boundary_data(w.obj(id), boundary);

    w.addForce(std::move(elastic));
    w.addForce(std::move(bending));
    w.addForce(std::move(Pr), id);
    w.setForceApplier(StaticForceApplierP(delta));
    int freq = 200;
//    if (scenario == 1) {
        w.setStepSimulationAlgo(std::move(NewtonSolverStepAlgo()));
        freq = 1;
//    }
    auto time = World3d::Timer();
    w.Simulation([&id, &time,  &err, &out = f, p0 = p, &to_save, &prefix, &vc, &freq](StepSimInfo& info, World* w)->bool{
        static double resid_init;
        int it = info.it;
        if (it == 1) resid_init = w->getWorldNorm();//w->getWorldShift();
        if (it > 0 && it % freq == 0){
            double resid = w->getWorldNorm();//w->getWorldShift();
            double eps = resid / resid_init;
            out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
            auto& xc = w->obj(id).m_x[vc];
            out << "\tcenter " << xc << std::endl;
            if (! (it % 400000)) {
                string fname = to_save + "/" + prefix + "benchmarkKyriacou_" + to_string(p0) + "_" + to_string(it) +".stl";
//                w->obj(id).save(fname);
            }
            if (eps < err || resid < err){
                out << "Algorithm is converged: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                string fname = to_save + "/" + prefix + "benchmarkKyriacou_" + to_string(p0) + "_converged.stl";
                w->obj(id).save(fname);
                return true;
            }
            if (eps > 300 || std::isnan(eps)){
                out << "Algorithm is diverged: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                string fname = to_save + "/" + prefix + "benchmarkKyriacou_" + to_string(p0) + "_diverged.stl";
                w->obj(id).save(fname);
                exit(-1);
            }
        }

        return false;
    });
    t.join();

    return 0;
}

int Benchmark4(int argc, char* argv[]){
    thread t(start_debug_gui, argc, argv);

    double force = 1000;
    double delta = 1.6e-6;//1.8e-6//3.8/2
    double err = 1e-5;
    double mesh_h = 0.5;//0.2495;
    string to_save = "../result";
    string prefix = "";
    int scenario = 1;

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << std::endl;
            std::cout << " Command line options: " << std::endl;
            std::cout << "  -f, --force    <Force value>" << std::endl;
            std::cout << "  -d, --delta    <Solver parameter>" << std::endl;
            std::cout << "  -e, --err      <Parameter of stop criteria>" << std::endl;
            std::cout << "  -h, --mesh_h   <Step of generated mesh>" << std::endl;
            std::cout << "  -t, --target   <Directory to save results>" << std::endl;
            std::cout << "  -p, --prefix   <Prefix to file name>" << std::endl;
            std::cout << "  -s, --scenario <Choosing of the running scenario: 0 - SVK, 1 - MR, 2 - MY, 3 - GOH>" << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            force = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--delta") == 0) {
            delta = atof(argv[i + 1]);
            i++;
            continue;
        }

        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--err") == 0) {
            err = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--mesh_h") == 0) {
            mesh_h = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0) {
            to_save = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) {
            prefix = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--scenario") == 0) {
            scenario = atoi(argv[i + 1]);
            i++;
            continue;
        }
    }
    std::cout << "Input parameters: " <<
              "\n\tforce = " << force <<
              "\n\tdelta = " << delta <<
              "\n\terr = " << err <<
              "\n\tmesh_h = " << mesh_h <<
              "\n\tscenario = " << scenario <<
              "\n\tdirectory = " << to_save << std::endl;
    {
        path path = to_save;
        auto s = status(path);
        if (!status_known(s) || s.type() != file_type::directory)
            create_directory(path);
    }

    string fname = to_save + "/" + prefix + "benchmark4_" + to_string(force) + "_log.txt";
    //std::ofstream f(fname);
    auto& f = std::cout;
//    if (mesh_h > 0.5) {
//        std::cout << "Too big mesh step!!! Stoped";
//        f << "Too big mesh step!!! Stoped";
//        exit(-1);
//    }

    f << "Input parameters: " <<
      "\n\tforce = " << force <<
      "\n\tdelta = " << delta <<
      "\n\terr = " << err <<
      "\n\tmesh_h = " << mesh_h <<
      "\n\tscenario = " << scenario <<
      "\n\tdirectory = " << to_save << std::endl;

    double a = 10, b = 10;
    static const int LEFT = 1 << 0, RIGHT = 1 << 1, BOTTOM = 1 << 2, TOP = 1 << 3, OXALIGN = 1 << 4, OYALIGN = 1 << 5;
    AniMesh am;
    if (true) {
         am = generate_dual_aligned_rectangle(a/100, b/100, mesh_h/100);
    } else  {
        am = generate_dual_aligned_rectangle(a, b, mesh_h);
    }
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = "rectangle";
    obj.set_is_movable_checker([] (Object3D& obj, const V_ind& v) -> bool
                               { return !(obj.m_boundary[v] & (LEFT | RIGHT | BOTTOM | TOP)); });
    SimplePressureLoad Pr(force);

    double H = 0.000174;//H = 0.05434;
    string dir = "../generated";
    Force elastic;
    Force bending;
    if (scenario == 0) {
        delta *= 100;
        double lambda = 51.63, mu = 1.59 * 1e4;
        elastic = SVKirchhoffModel(H, lambda, mu, dir, true);
        bending = BendingForce(elastic.target<SVKirchhoffModel>()->f, true);
    } else if (scenario == 1){
        delta = 4e-6*100;//10;
        err = 1e-10;
        double H = 0.0002;
        double mu = 2.45e4, alpha = 0.0245, beta = 0.0297;
        std::array<double, 3> A = {1, 0, 0};
        elastic = MooneyRivlinModel2(H, mu, alpha, beta, A, dir, true);
//        elastic.target<HyperElasticForceBase>()->prepareJacobianFunction();
        bending = BendingForce(elastic.target<HyperElasticForceBase>()->f, true);
//        bending.target<BendingForce>()->prepareJacobianFunction();
    } else if (scenario == 2){
        delta = 1e-3;
        err = 1e-7;
        double c0 = 5.95e6, c1 = 1.48e-3, c2 = 6.74e-4;
        std::array<double, 3> A = {1, 0, 0};
        elastic = MayYinModel(H, c0, c1, c2, A, dir, true);
        bending = BendingForce(elastic.target<MayYinModel>()->f, true);
    } else if (scenario == 3){
        delta = 1.0e-3;
        err = 1e-7;
//        double H = 0.00015;
        double c = 0.0511e6, k1 = 0.015e6, k2 = 0.0418, kappa = 0.05;
        std::array<double, 3> A = {1, 0, 0};
        elastic = HGOModel(H, c, k1, k2, kappa, A, dir, true);
        bending = BendingForce(elastic.target<HGOModel>()->f, true);
    } else if (scenario == 4){
        delta = 0.1e-4;
//        double mu = 2.45e4, alpha = 0.0297, beta = 0.347, teta = 0.347;
        H = 0.01;
        double mu = 1, alpha = 0.0245, beta = 0.0297, teta = 0.347;
        std::array<double, 3> A = {cos(teta), sin(teta), 0};
        elastic = MooneyRivlinModel(H, mu, alpha, beta, A, dir, true);
        bending = BendingForce(elastic.target<MooneyRivlinModel>()->f, true);
    }
    map<E_ind, BendingForce::BoundaryData::Entry> boundary;
    std::array<vector<E_ind>, 2> align;
    for (auto e: obj.m_mesh.edges()){
        auto f = face_around(obj.m_mesh, e);
        auto v = vert_around(obj.m_mesh, e);
        int amask[2] = {OXALIGN, OYALIGN};
        for (int i = 0; i < 2; ++i) {
            if ((obj.m_boundary[v[0]] & amask[i]) && (obj.m_boundary[v[1]] & amask[i])) {
                align[i].push_back(e);
                break;
            }
        }
        if (f[0].second && f[1].second)
            continue;
        int mask[] = {LEFT, RIGHT, BOTTOM, TOP};
        DReal vec[4][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}};//{{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}};
        for (int i = 0; i < 4; ++i) {
            if ((obj.m_boundary[v[0]] & mask[i]) && (obj.m_boundary[v[1]] & mask[i])) {
                boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, &vec[i][0]);
                break;
            }
        }
    }
    V_ind vc;
    for (auto v: obj.m_mesh.vertices()){
        if ((obj.m_boundary[v] & (OXALIGN | OYALIGN)) == (OXALIGN | OYALIGN)){
            vc = v;
            break;
        }
    }

    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto id = w.addObject3D(move(obj), 1);

    auto comp_aligns = [align, &obj = w.obj(id)](){
        std::array<double, 2> res = {0, 0};
        for (int i = 0; i < 2; ++i){
            for (auto e: align[i]) {
                auto v = vert_around(obj.m_mesh, e);
                res[i] += sqrt(Vector(obj.m_x[v[0]], obj.m_x[v[1]]).squared_length());
            }
        }
        return res;
    };

//    elastic.target<HyperElasticForceBase>()->prepareJacobianFunction();
//    bending.target<BendingForce>()->prepareJacobianFunction();
    bending.target<BendingForce>()->set_boundary_data(w.obj(id), boundary);
    w.addForce(std::move(elastic));
    w.addForce(std::move(bending));
    w.addForce(std::move(Pr), id);
    w.setForceApplier(StaticForceApplierP(delta));
    int freq = 200;
//    w.setStepSimulationAlgo(std::move(NewtonSolverStepAlgo()));
//    freq = 1;
    auto time = World3d::Timer();
    w.Simulation([&id, &time, &force, &err, &out = f, &to_save, &prefix, &delta, &comp_aligns, &vc, &freq](StepSimInfo& info, World* w)->bool{
        static double resid_init;
        int it = info.it;
        if (it == 1) resid_init = w->getWorldShift();
        if (it > 0 && it % freq == 0){
            double resid = w->getWorldShift();
            double eps = resid / resid_init;
            out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
            auto lens = comp_aligns();
            auto& xc = w->obj(id).m_x[vc];
            std::cout << "\taxis lens OX: " << lens[0] << "; OY: " << lens[1] << " :: center " << xc << std::endl;
            if (! (it % 400000)) {
                string fname = to_save + "/" + prefix + "benchmark4_" + to_string(force) + "_" + to_string(it) +".stl";
//                w->obj(id).save(fname);
            }
            if (eps < err || resid < err){
                out << "Algorithm is converged: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                string fname = to_save + "/" + prefix + "benchmark4_" + to_string(force) + "_converged.stl";
                w->obj(id).save(fname);
                return true;
            }
            if (eps > 300 || std::isnan(eps)){
                out << "Algorithm is diverged: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                string fname = to_save + "/" + prefix + "benchmark4_" + to_string(force) + "_diverged.stl";
                w->obj(id).save(fname);
                exit(-1);
            }
        }

        return false;
    });
    t.join();

    return 0;
}

#include "Solvers/NonLinearSolverKinsol.h"
#include "NSWorldWrapper.h"
#define TESTKINSOL
int BenchmarkBendAnnulus(int argc, char* argv[]){
    thread t(start_debug_gui, argc, argv);

    double mesh_h = 0.2;
    string to_save = "../result/benchmark2";
    string prefix = "";
    string test_name = "bBA";
    double force = 0.8;
    double err = 1e-5;

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << std::endl;
            std::cout << "Command line options: " << std::endl;
            std::cout << "-f, --force <Force value>" << std::endl;
            std::cout << "-d, --delta <Solver parameter>" << std::endl;
            std::cout << "-e, --err <Parameter of stop criteria>" << std::endl;
            std::cout << "-mh,--mesh_h <Step of generated mesh>" << std::endl;
            std::cout << "-t, --target <Directory to save results>" << std::endl;
            std::cout << "-p, --prefix <Prefix to file name>" << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            force = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--err") == 0) {
            err = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-mh") == 0 || strcmp(argv[i], "--mesh_h") == 0) {
            mesh_h = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0) {
            to_save = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) {
            prefix = argv[i + 1];
            i++;
            continue;
        }
    }
    std::cout << "Input parameters: " <<
              "\n\tforce = " << force <<
              "\n\terr = " << err <<
              "\n\tmesh_h = " << mesh_h <<
              "\n\tdirectory = " << to_save << std::endl;
    {
        path path = to_save;
        auto s = status(path);
        if (!status_known(s) || s.type() != file_type::directory)
            create_directory(path);
    }
    string fname = to_save + "/" + prefix + test_name + to_string(force) + "_" + to_string(mesh_h)+ "_log.txt";
    std::ofstream f(fname);

    f << "Input parameters: " <<
      "\n\tP/P_max = " << force/0.8 <<
      "\n\terr = " << err <<
      "\n\tmesh_h = " << mesh_h <<
      "\n\tdirectory = " << to_save << std::endl;

    double R1 = 6, R2 = 10;
    AniMesh am = generate_cutted_annulus(R1, R2, mesh_h);
    static const int CLAMPED = 4, PRESSURED = 2, FREE = 1;
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = "cut_annulus";
    obj.set_is_movable_checker([] (Object3D& obj, const V_ind& v) -> bool
                               { return !(obj.m_boundary[v] & CLAMPED); });
//    obj.save(to_save + "/" + test_name + "_init_mesh.stl");

    SimplePressureEdgeLoad load({0, 0, force*(R2 - R1)},
                                [](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] == PRESSURED); },
                                [](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] & PRESSURED); } );

    double E = 21e6, nu = 0;
    double lambda = E * nu / (1 - nu * nu), mu = E / (2 * (1 + nu));
    f << "SVK: lambda = " << lambda << ", mu = " << mu << std::endl;
    double H = 0.03;
    string dir = "../generated";
    SVKirchhoffNoPoissonModel svk(H, lambda, mu, dir, false);
    svk.prepareJacobianFunction(false);
    BendingForce svkbf(svk.f, false);
    svkbf.prepareJacobianFunction(false);
    map<E_ind, BendingForce::BoundaryData::Entry> boundary;
    for (auto e: obj.m_mesh.edges()){
        auto f = face_around(obj.m_mesh, e);
        if (f[0].second && f[1].second)
            continue;
        auto v = vert_around(obj.m_mesh, e);
        if ((obj.m_boundary[v[0]] & CLAMPED) && (obj.m_boundary[v[1]] & CLAMPED))
            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, {0, -1, 0});
        else if (obj.m_boundary[v[0]] != 0 && obj.m_boundary[v[1]] != 0)
            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
    }
    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto id = w.addObject3D(move(obj), 1);

    svkbf.set_boundary_data(w.obj(id), boundary);
    w.addForce(std::move(svk));
    w.addForce(std::move(svkbf));
    w.addForce(std::move(load), id);
    auto time = World3d::Timer();
    V_ind wv[2]; int ii = 0;
    for (auto& v: w.obj(id).m_mesh.vertices())
        if(w.obj(id).m_boundary[v] == (2|1))
            wv[ii++] = v;
    if  (w.obj(id).m_x[wv[0]].x() > w.obj(id).m_x[wv[1]].x()) swap(wv[0], wv[1]);
#ifndef TESTKINSOL
    NewtonSolverStepAlgo nssa;
#ifdef USE_INMOST_SOLVERS
    nssa._solver.solver = LinearSolver(INMOST::Solver::INNER_MPTILUC);
    nssa._solver.solver.SetParameterReal("drop_tolerance", 2e-5);
    nssa._solver.solver.SetParameterReal("reuse_tolerance", 0.1 * 2e-5);
    nssa._solver.solver.SetParameterReal("absolute_tolerance", 1e-20);
    nssa._solver.solver.SetParameterReal("relative_tolerance", 1e-11);
    nssa._solver.solver.setVerbosityLevel(1);
#else
    nssa._solver.solver = std::move(LinearSolver("eigen"));
#endif
//    nssa.newton_algo = [](NewtonSolverStepAlgo& ns, World* w){
//        ns.updateMeshData(w);
//        int N = ns.getNodfs(w);
//        ns.assembleSystem(w, N);
//        double tau = 0.3;
//        static double resid_init = -1;
//        if (ns.m_it == 0) {
//            ns._solver.solver.SetParameterReal("absolute_tolerance", 1e-18);
//            ns._solver.solver.SetParameterReal("relative_tolerance", 1e-11);
//            double droptol = 2e-6;
//            resid_init = ns.computeRhsNorm(1);
//            std::cout << "\teps = " << resid_init << " drop_tol = " << droptol << std::endl;
//            ns._solver.solver.SetParameterReal("drop_tolerance", droptol);
//            ns._solver.solver.SetParameterReal("reuse_tolerance", 0.1 * droptol);
//        } else {
//            double eps = ns.computeRhsNorm(1) / resid_init;
//            double droptol = 1e-4/2;//8e-4;//4e-5;//1e-6;
//            if (eps > 180){
//                tau = 0.3;
//                droptol = 5e-5/2;//1e-4;
//            } else if (eps < 0.6){
//                tau = 1;
//            } else {
//                tau = 0.9;
//            }
//            std::cout << "\teps = " << eps << " drop_tol = " << droptol << std::endl;
//            ns._solver.solver.SetParameterReal("drop_tolerance", droptol);
//            ns._solver.solver.SetParameterReal("reuse_tolerance", 0.1 * droptol);
//        }
//        ns.Solve();
//        ns.applySolution(w, tau);
//        ns.m_it++;
//        return 0;
//    };
    nssa.newton_algo = [](NewtonSolverStepAlgo& ns, World* w){
        ns.updateMeshData(w);
        int N = ns.getNodfs(w);
        ns.assembleSystem(w, N);
        double tau = 1.0;
        static double resid_init = -1;
        if (ns.m_it == 0) {
            resid_init = ns.computeRhsNorm(2);
            std::cout << "\tabs_init = " << resid_init << std::endl;
        } else {
            double abs_resid = ns.computeRhsNorm(2);
            double eps =  abs_resid / resid_init;
            std::cout << "\teps = " << eps  << " abs = " << abs_resid << std::endl;
        }
        ns.Solve();
        ns.applySolution(w, tau);
        ns.m_it++;
        return 0;
    };
    w.setStepSimulationAlgo(std::move(nssa));
    int freq = 1;
    w.Simulation([&id, &time, &force, &err, &out = f, wv, &to_save, &prefix, &freq, &test_name, &mesh_h](StepSimInfo& info, World* w)->bool{
        static double resid_init;
        int it = info.it;
        if (it == 1) resid_init = w->getWorldNorm("v:force", 2, 2, 2);
        if (it > 0 && it % freq == 0){
            double resid = w->getWorldNorm("v:force", 2, 2, 2);
            double eps = resid / resid_init;
            out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
            std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
            out << "watch { " << w->obj(id).m_x[wv[0]].z() << " " << w->obj(id).m_x[wv[1]].z() << " }" << std::endl;
            std::cout << "watch { " << w->obj(id).m_x[wv[0]].z() << " " << w->obj(id).m_x[wv[1]].z() << " }" << std::endl;
            if (eps < err || resid < err){
                out << "Algorithm is converged: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                out << "watch { " << w->obj(id).m_x[wv[0]].z() << " " << w->obj(id).m_x[wv[1]].z() << " }" << std::endl;
                string fname = to_save + "/" + prefix + test_name + to_string(force) + "_" + to_string(mesh_h) +"_converged.stl";
                w->obj(id).save(fname);
                string main_data = to_save + "/" + test_name + ".csv";
                ofstream fdata(main_data, std::ios_base::app);
                fdata << mesh_h << ", " << force/0.8 << ", " << w->obj(id).m_x[wv[0]].z() << ", " << w->obj(id).m_x[wv[1]].z()
                      << ", " <<  it << ", " << time.elapsed() << "\n";
                fdata.close();
                return true;
            }
            if (it > 600 || std::isnan(eps) || info.status < 0){
                out << "Algorithm is diverged or reached maximum of time - iters: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                out << "watch { " << w->obj(id).m_x[wv[0]].z() << " " << w->obj(id).m_x[wv[1]].z() << " }" << std::endl;
                string fname = to_save + "/" + prefix + test_name + to_string(force) + "_" + to_string(mesh_h) + "_diverged.stl";
                w->obj(id).save(fname);
                return true;
            }
        }

        return false;
    });
#else
    NSWorldWrapper nsww(w);
    NLProblem nlp;
    nlp.m_dofs = nsww.getNodfs();
    nlp.m_rhs = [&nsww](const double *x, double *b) -> int{
        return nsww.Residual(x, b);
    };
    nlp.m_jac = [&nsww](const double *x, SparseMatrix *sm) -> int{
        return nsww.Jacobian(x, sm);
    };

    NonLinearSolverKinsol nlsp;
    nlsp.setProblem(nlp);
#ifdef USE_INMOST
    LinearSolver ls(INMOST::Solver::INNER_MPTILUC);
#else
    LinearSolver ls;
#endif
    ls.setVerbosityLevel(1);
    double droptol = 2e-5;
    ls.SetParameterReal("drop_tolerance", droptol);
    ls.SetParameterReal("reuse_tolerance", 0.1 * droptol);
    ls.SetParameterReal("absolute_tolerance", 1e-12);
    ls.SetParameterReal("relative_tolerance", 1e-20);
    nlsp.SetLinearSolver(&ls);
    static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
    nlsp.SetInfoHandlerFn([](const char *module, const char *function, char *msg){
        std::cout << "[" << module << "] " << function << "\n   " << msg << std::endl;
    });
    nlsp.SetVerbosityLevel(1);
    nlsp.SetInitialGuess(nsww.getCurrentX());
    nlsp.setInterIterationMonitorFunc([&nsww, &w, &id, wv, &time](const double * x){
        nsww.setX(x);
        std::cout << "watch { " << w.obj(id).m_x[wv[0]].z() << " " << w.obj(id).m_x[wv[1]].z() << " }; time = " << time.elapsed() << std::endl;
        nsww.RenderIteration();
        return 0;
    });
    nlsp.SetFuncNormTol(err * 0.0367);
    nlsp.SetScaledStepTol(1e-7);
    nlsp.SetMaxNewtonStep(R2 / 100 * sqrt(nlp.m_dofs));
    bool slvFlag = false;

    nlsp.SetMAA(0);
    nlsp.SetSolveStrategy(NonLinearSolverKinsol::PICARD);
    nlsp.SetNumMaxIters(1);
    slvFlag = nlsp.Solve();
    nlsp.SetMAA(0);

    if (!slvFlag) {
        nlsp.SetNumMaxIters(200);
        nlsp.SetParameterInt("MaxSetupCalls", 1);
        nlsp.SetParameterInt("MaxSubSetupCalls", 1);
        nlsp.SetScaledStepTol(1e-9);
        nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
        nlsp.SetMaxBetaFails(40);
//    nlsp.SetEtaForm(NonLinearSolverKinsol::CONSTANT);
        slvFlag = nlsp.Solve();
    }
    if (!slvFlag) {
        nlsp.SetNumMaxIters(300);
        nlsp.SetParameterInt("MaxSetupCalls", 1);
        nlsp.SetParameterInt("MaxSubSetupCalls", 1);
        nlsp.SetSolveStrategy(NonLinearSolverKinsol::NONE);
        slvFlag = nlsp.Solve();
    }
    std::cout <<"\t#NonLinIts = "<<nlsp.GetNumNolinSolvIters() << " Residual = "<< nlsp.GetResidualNorm()
              << "\n\t#linIts = " << nlsp.GetNumLinIters() << " #funcEvals = " << nlsp.GetNumLinFuncEvals() << " #jacEvals = " << nlsp.GetNumJacEvals()
              << "\n\t#convFails = " << nlsp.GetNumLinConvFails() << " #betaCondFails = " << nlsp.GetNumBetaCondFails() << " #backtrackOps = " << nlsp.GetNumBacktrackOps() << "\n";
    std::cout << "Reason: " << nlsp.GetReason() << std::endl;
#endif

    t.join();

    return 0;
}
#undef TESTKINSOL

int Benchmark2(int argc, char* argv[]){
    double force = 4;
    double delta = 1.6e-6;//1.8e-6//3.8/2
    double err = 1e-5;
    double mesh_h = 0.012/2;//0.2495;//0.5;//0.2495;
    string to_save = "../result/benchmark2";
    string test_name = "bRС";
    string prefix = "6";
    string dir = "../generated6";
    double R = 1;

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << std::endl;
            std::cout << "Command line options: " << std::endl;
            std::cout << "-f, --force <Force value>" << std::endl;
            std::cout << "-d, --delta <Solver parameter>" << std::endl;
            std::cout << "-e, --err <Parameter of stop criteria>" << std::endl;
            std::cout << "-mh,--mesh_h <Step of generated mesh>" << std::endl;
            std::cout << "-t, --target <Directory to save results>" << std::endl;
            std::cout << "-p, --prefix <Prefix to file name>" << std::endl;
            std::cout << "-b, --width <Width of bulk>" << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            force = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--delta") == 0) {
            delta = atof(argv[i + 1]);
            i++;
            continue;
        }

        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--err") == 0) {
            err = atof(argv[i + 1]);
            i++;
            continue;
        }

        if (strcmp(argv[i], "-mh") == 0 || strcmp(argv[i], "--mesh_h") == 0) {
            mesh_h = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0) {
            to_save = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) {
            prefix = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--width") == 0) {
            R = atof(argv[i + 1]);
            i++;
            continue;
        }
    }
    std::cout << "Input parameters: " <<
              "\n\tP/Pmax = " << force/4 <<
              "\n\tdelta = " << delta <<
              "\n\terr = " << err <<
              "\n\tmesh_h = " << mesh_h <<
              "\n\tR = " << R <<
              "\n\tdirectory = " << to_save << std::endl;
    {
        path path = to_save;
        auto s = status(path);
        if (!status_known(s) || s.type() != file_type::directory)
            create_directory(path);
    }
    string fpref = to_save + "/" + prefix + test_name + "_" + to_string(force) + "_" + to_string(mesh_h) + "_" + to_string(static_cast<int>(R));
    string fname = fpref + "_log.txt";
    std::ofstream f(fname);

    AniMesh am = generate_half_sphere(R, mesh_h);
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = "half_sphere";
    obj.save(fpref + "_inited.stl");

    double E = 1.2e6, nu = 0;
    double lambda = E * nu / (1 - nu * nu), mu = E / (2 * (1 + nu));
    double H = 0.1;
    SVKirchhoffModel svk(H, lambda, mu, dir, false);
    BendingForce svkbf(svk.f, false);
    static const int CLAMPED = 2, FREE = 1;
    map<E_ind, BendingForce::BoundaryData::Entry> boundary;
    for (auto e: obj.m_mesh.edges()){
        auto f = face_around(obj.m_mesh, e);
        if (f[0].second && f[1].second)
            continue;
        auto v = vert_around(obj.m_mesh, e);
        if ((obj.m_boundary[v[0]] & CLAMPED) && (obj.m_boundary[v[1]] & CLAMPED))
            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, {0, -1, 0});
        else if ((obj.m_boundary[v[0]] & FREE) && (obj.m_boundary[v[1]] & FREE))
            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
    }

    World w;
//    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto id = w.addObject3D(move(obj), 1);

    svkbf.set_boundary_data(w.obj(id), boundary);
    w.addForce(std::move(svk));
    w.addForce(std::move(svkbf));
    w.setForceApplier(StaticForceApplierP(delta)); //устойчивость неогука до 6.0е-6
    int freq = 20000;
    w.obj(id).setDirichletCond([](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & CLAMPED) return 0;
        return 7;
    });

    w.Simulation([](StepSimInfo& info, World* w)->bool{
        static int nit = 0;
        return nit++ == 1;
    });
    auto k0 = w.obj(id).m_mesh.property_map<F_ind, std::array<DReal, 3>>("f:BendingForce:k0");
    LocalCartesian S(obj, false);
    static double n1[3] = {0}, n2[3] = {0};
    int cnt = 0;
    for (auto f: w.obj(id).m_mesh.faces()){
        auto v = vert_around(w.obj(id).m_mesh, f);
        auto vo = vert_opposite(w.obj(id).m_mesh, f);
        if (vo[0].second && vo[1].second && vo[2].second) {
            std::array<double, 3> dif = {1.0/R - k0.first[f][0], 1.0/R - k0.first[f][1],  k0.first[f][2]};
            for (int i = 0; i < 3; ++i){
                if (fabs(dif[i]) > n1[i]) n1[i] = fabs(dif[i]);
                n2[i] += fabs(dif[i]);
            }
            cnt++;
        }
        if (!vo[0].second && !vo[1].second && !vo[2].second)
            throw runtime_error("Mesh is bad");
    }
    for (int i = 0; i < 3; ++i) n2[i] /= cnt;
    std::cout << "Watch internal vals: mesh_h = " << mesh_h << std::endl;
    std::cout << "Max_deflection = " << n1[0] << " " << n1[1] << " " << n1[2] << std::endl;
    std::cout << "Mid deflection = " << n2[0] << " " << n2[1] << " " << n2[2] << std::endl;

#ifdef USE_INMOST_SOLVERS
    {
        auto& ob = w.obj(id);
        using namespace INMOST;
        INMOST::Mesh* mesh = new INMOST::Mesh();
        vector<HandleType> new_nodes(ob.m_mesh.num_vertices());
        vector<HandleType> new_faces(ob.m_mesh.num_faces());
        for (auto v: ob.m_mesh.vertices()){
            std::array<Storage::real,3> xyz{(Storage::real)ob.m_x[v][0], (Storage::real)ob.m_x[v][1], (Storage::real)ob.m_x[v][2]};
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
        Tag kk = mesh->CreateTag("dk", DATA_REAL, NODE, FACE, 3);
        for (auto n = mesh->BeginNode(); n != mesh->EndNode(); ++n){
            auto ff = face_around(ob.m_mesh, V_ind(n->LocalID()));
            std::array<double, 3> res = {0};
            for (auto f: ff){
                res[0] += 1.0/R - k0.first[f][0];
                res[1] += 1.0/R - k0.first[f][1];
                res[2] += k0.first[f][2];
            }
            for (int i = 0; i < 3 && !ff.empty(); ++i) res[i] /= ff.size();
            for (int i = 0; i < 3; ++i) n->RealArray(kk)[i] = res[i];
        }

        mesh->SaveVTK(fpref + "_kk.vtk");
        delete mesh;
    }
#endif

    return 0;
}

int OneAxisStretch(int argc, char* argv[]){
    double P = 1;
    double l = 1, b = 1, H = 1;
    double E = 2e11, nu = 0.3;
    double mesh_h = 0.05;
    string to_save = "../result/LinVsNonlin/benchs/oneaxis";
    string dir = "../generated";
    bool with_view = true;
    bool with_bnd = true;
//    double delta = 1e-5;
    const int SVK_m = 1, LIN_m = 0;
    int model_type = SVK_m;

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            std::cout << "Help message: " << std::endl;
            std::cout << "Command line options: " << std::endl;
            std::cout << "-P <Force value>" << std::endl;
            std::cout << "-l <Length of plate>" << std::endl;
            std::cout << "-b <Width of plate>" << std::endl;
            std::cout << "-H <Thickness of plate>" << std::endl;
            std::cout << "-E <Young modulus>" << std::endl;
            std::cout << "-nu <Poisson modulus>" << std::endl;
            std::cout << "-mh <Mesh step>" << std::endl;
//            std::cout << "-delta <Relaxate solver parameter>" << std::endl;
            std::cout << "-model <Choose SVK or LIN>" << std::endl;
            std::cout << "--no_view <Without view>" << std::endl;
            std::cout << "--no_bnd <Without bending force>" << std::endl;
            std::cout << "-t, --target <Directory to save results>" << std::endl;
            exit(0);
        } else if (std::string("-P") == argv[i]){
            if (i + 1 < argc) P = atof(argv[++i]);
            continue;
        } else if (std::string("-l") == argv[i]){
            if (i + 1 < argc) l = atof(argv[++i]);
            continue;
        } else if (std::string("-b") == argv[i]){
            if (i + 1 < argc) b = atof(argv[++i]);
            continue;
        } else if (std::string("-H") == argv[i]){
            if (i + 1 < argc) H = atof(argv[++i]);
            continue;
        } else if (std::string("-E") == argv[i]){
            if (i + 1 < argc) E = atof(argv[++i]);
            continue;
        } else if (std::string("-nu") == argv[i]){
            if (i + 1 < argc) nu = atof(argv[++i]);
            continue;
        } else if (std::string("-mh") == argv[i]){
            if (i + 1 < argc) mesh_h = atof(argv[++i]);
            continue;
//        } else if (std::string("-delta") == argv[i]){
//            if (i + 1 < argc) delta = atof(argv[++i]);
//            continue;
        } else if (std::string("-model") == argv[i]){
            if (i + 1 < argc){
                if (std::string("SVK") == argv[i+1])
                    model_type = SVK_m;
                else if (std::string("LIN") == argv[i+1])
                    model_type = LIN_m;
                else throw std::runtime_error("Unknown model name: \"" + std::string(argv[i+1]) + "\"");
                ++i;
            }
            continue;
        } else if (std::string("--no_view") == argv[i]){
            with_view = false;
            continue;
        } else if (std::string("--no_bnd") == argv[i]){
            with_bnd = false;
            continue;
        } else if (std::string("-t") == argv[i] || std::string("--target") == argv[i]){
            if (i + 1 < argc) to_save = argv[++i];
            continue;
        } else {
            throw std::runtime_error("Unknown option \"" + std::string(argv[i]) + "\"");
        }
    }

    thread* t = nullptr;
    if (with_view) {
        t = new thread(start_debug_gui, argc, argv);
    }

    AniMesh am = generate_rectangle(b/2, l/2, mesh_h, 1);
    for (int i = 0; i < am.vertices.size()/3; ++i) am.vertices[3*i+0] += b/4, am.vertices[3*i+1] += l/4;
    static const int OX = 1, OY = 2, AX = 4, AY = 8;
    for (int i = 0; i < am.vertices.size()/3; ++i) {
        am.vertex_label[i] = 0;
        if (fabs(am.vertices[3*i+1] - 0) < mesh_h*1e-5)  am.vertex_label[i] |= OX;
        if (fabs(am.vertices[3*i+0] - 0) < mesh_h*1e-5)  am.vertex_label[i] |= OY;
        if (fabs(am.vertices[3*i+1] - l/2) < mesh_h*1e-5)  am.vertex_label[i] |= AX;
        if (fabs(am.vertices[3*i+0] - b/2) < mesh_h*1e-5)  am.vertex_label[i] |= AY;
    }

    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};

    SimplePressureEdgeLoad load({0, P, 0},
                                [](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] & AX) != 0; },
                                [](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] & AX) != 0; });
    double lambda = E * nu / (1 - nu * nu), mu = E / (2 * (1 + nu));
    Force svk;
    switch (model_type) {
        case LIN_m:
            svk = LinearElasticModel1(H, lambda, mu, dir,  false);
            break;
        case SVK_m:
            svk = SVKirchhoffModel(H, lambda, mu, dir,  false);
            break;
    }
    svk.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
    BendingForce svkbf(svk.target<HyperElasticForceBase>()->f, false);
    svkbf.prepareJacobianFunction(false);

    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto id = w.addObject3D(move(obj), 1);
    w.obj(id).setDirichletCond([](Object3D& obj, const V_ind& v) -> unsigned char{
        unsigned char res = (0x7 & ((obj.m_boundary[v] & OX) ? 0x5 : 0x7)) & ((obj.m_boundary[v] & OY) ? 0x6 : 0x7);
        return res;
    });

    w.addForce(std::move(svk));
    if (with_bnd) w.addForce(std::move(svkbf));
    w.addForce(std::move(load), id);

    NSWorldWrapper nsww(w);
    NLProblem nlp = NSWorldWrapper::makeProblem(nsww);
    LinearSolver ls("inner_mptiluc");
    ls.setVerbosityLevel(1);
    ls.SetParameterReal("relative_tolerance", 1e-20);
    double droptol = 1e-4;
    ls.SetParameterReal("drop_tolerance", droptol);
    ls.SetParameterReal("reuse_tolerance", 0.1 * droptol);
    auto monitor_fcn = [&nsww](const double *x) {
        nsww.setX(x);
        nsww.RenderIteration();
        return 0;
    };

    World3d::Timer sym_time, com_time;
//    //initialize KINSOL solver
    NonLinearSolverKinsol nlsp(nlp);
    nlsp.SetLinearSolver(&ls);
    static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
    nlsp.SetVerbosityLevel(1);
    nlsp.SetInitialGuess(nsww.getCurrentX());
    nlsp.SetFuncNormTol(1e-20);
    nlsp.SetMaxNewtonStep(10*l * sqrt(nlp.m_dofs));
    nlsp.SetParameterInt("MaxSetupCalls", 1);
    nlsp.SetParameterInt("MaxSubSetupCalls", 1);
    nlsp.SetScaledStepTol(1e-5);
    nlsp.SetMaxBetaFails(10);
    nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
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
    bool slvFlag = nlsp.Solve();
    std::cout << result_str_nlsp(nlsp) << std::endl;

    double ax_dy = 0, ay_dx = 0;
    {
        int n_ax = 0, n_ay = 0;
        auto& obj = w.obj(id);
        for (auto v: obj.m_mesh.vertices()){
            if (obj.m_boundary[v] & AX) ax_dy += obj.m_x[v][1] - obj.m_x0[v][1], ++n_ax;
            if (obj.m_boundary[v] & AY) ay_dx += obj.m_x[v][0] - obj.m_x0[v][0], ++n_ay;
        }
        ax_dy /= n_ax, ay_dx /= n_ay;
    }

    std::cout << "Deflections (U, V, 0): U = " << ay_dx << " V = " << ax_dy << std::endl;
    w.obj(id).save(to_save + "/" + "res.vtk");


    if (with_view) {
        t->join();
        delete t;
    }
    return 0;
}

int Benchmark1(int argc, char* argv[]){
//    computeShiftOperator();


    double force = 4;
    double delta = 1.6e-6;//1.8e-6//3.8/2
    double err = 1e-5;
    double mesh_h = 0.1;//0.2495;//0.5;//0.2495;
    string to_save = "../result/benchmark1N";
    string test_name = "bRN";
    string prefix = "uni";
    string dir = "../generated";
    bool view = true;
    bool no_bnd = false;
    double a = 10;
    double b = 1;
    double H = 0.1;
    const int SVK_MODEL = 0, LIN_MODEL = 1;
    int model_type = SVK_MODEL;
    double E = 1.2e6, nu = 0;
    std::cout << "prefix = " << prefix << std::endl;

    for (int i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                std::cout << "Help message: " << std::endl;
                std::cout << "Command line options: " << std::endl;
                std::cout << "-f, --force <Force value>" << std::endl;
                std::cout << "-d, --delta <Solver parameter>" << std::endl;
                std::cout << "-e, --err <Parameter of stop criteria>" << std::endl;
                std::cout << "-mh,--mesh_h <Step of generated mesh>" << std::endl;
                std::cout << "-t, --target <Directory to save results>" << std::endl;
                std::cout << "-p, --prefix <Prefix to file name>" << std::endl;
                std::cout << "-b, --width <Width of bulk>" << std::endl;
            exit(0);
        }
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            force = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--delta") == 0) {
            delta = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "--etype") == 0) {
            if (argv[i + 1] == std::string("lin"))
                model_type = LIN_MODEL;
            else if (argv[i + 1] == std::string("svk"))
                model_type = SVK_MODEL;
            else
                throw std::runtime_error("Unknown model " + std::string(argv[i + 1]) + R"(. Available "lin", "svk")");
            i++;
            continue;
        }
        if (strcmp(argv[i], "--nobnd") == 0){
            no_bnd = true;
            continue;
        }
        if (strcmp(argv[i], "-pr") == 0 || strcmp(argv[i], "--params") == 0){
            E = atof(argv[i + 1]);
            nu = atof(argv[i + 2]);
            i += 2;
        }

        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--err") == 0) {
            err = atof(argv[i + 1]);
            i++;
            continue;
        }

        if (strcmp(argv[i], "-mh") == 0 || strcmp(argv[i], "--mesh_h") == 0) {
            mesh_h = atof(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--target") == 0) {
            to_save = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--prefix") == 0) {
            prefix = argv[i + 1];
            i++;
            continue;
        }
        if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--width") == 0) {
            b = atof(argv[i + 1]);
            i++;
            continue;
        }
    }
    std::cout << "Input parameters: " <<
        "\n\tP/Pmax = " << force/4 <<
        "\n\tdelta = " << delta <<
        "\n\terr = " << err <<
        "\n\tmesh_h = " << mesh_h <<
        "\n\t b = " << b <<
        "\n\tdirectory = " << to_save << std::endl;
    {
        path path = to_save;
        auto s = status(path);
        if (!status_known(s) || s.type() != file_type::directory)
            create_directory(path);
    }
    string fpref = to_save + "/" + prefix + test_name + "_" + to_string(force) + "_" + to_string(mesh_h) + "_" + to_string(static_cast<int>(b));
    string fname = fpref + "_log.txt";
    std::ofstream f(fname);
//    auto& f = std::cout;

#ifdef BEND_WATCH
    {
        std::vector<std::array<double, 4>>& KK = g_Kproxy.data;
        string RscriptFMT = "#This is generated code! You shouldn't modify it\n"
                            "rm(list = ls())\n"
                            "library(nleqslv); library(Carlson)\n"
                            "\n"
                            "P <- %lg; L <- %lg;\n"
                            "N <- %d\n"
                            "q <- P^0.5; r <- q/L;\n"
                            "fn <- function(x) { abs( elliptic_F(pi/2, x[1]^2) - elliptic_F(asin(1/(2^0.5 * x[1])), x[1]^2) - q ) }\n"
                            "z <- nleqslv(c(0.8), fn)\n"
                            "k <- z$x\n"
                            "phi0 <- asin(1/(2^0.5*k)); phi1 <- pi/2;\n"
                            "phi <- seq(from=phi0, to=phi1, by=(phi1-phi0)/N)\n"
                            "curvature <- function(x) { list(\n"
                            "  \"%s\"=Re(elliptic_F(x, k^2) - elliptic_F(phi0, k^2))/r, \n"
                            "  \"%s\"=2*k*r*cos(x),\n"
                            "  \"x\"=1/r*(-2*k*cos(x) + (2*(2*k^2 - 1))^0.5),\n"
                            "  \"z\"=1/r*Re(elliptic_F(x, k^2) - elliptic_F(phi0, k^2) - 2*(elliptic_E(x, k^2) - elliptic_E(phi0, k^2)))\n"
                            "  ) }\n"
                            "res <- curvature(phi)\n"
                            "write.csv(res, file = \"%s\", row.names = FALSE)\n";
        string Rscript;
        Rscript.resize(RscriptFMT.size() + 1024);
        const int kNSubSteps = 10000, kNSavedPrms = 4;
        string RSource = to_save + "/comp_curvature.R";
        string CurveRes = to_save + "/curvature.csv";
        remove(CurveRes);
        sprintf(Rscript.data(), RscriptFMT.c_str(), force, a, kNSubSteps, "s", "dtds", CurveRes.c_str());
        std::ofstream Rf(RSource);
        Rf << Rscript;
        Rf.close();
        int sysres = system(("Rscript " + RSource).c_str());
        if (!exists(CurveRes)) {
            throw std::runtime_error("In generated R code occured errors");
        } else {
            std::ifstream Rf(CurveRes);
            KK.reserve(kNSubSteps + 1);
            string s;
            char comma;
            Rf >> s;
            while (!Rf.eof()) {
                std::array<double, kNSavedPrms> dat = {0, 0};
                Rf >> dat[0] >> comma >> dat[1] >> comma >> dat[2] >> comma >> dat[3];
                dat[0] -= a/2;      // shift data to our mesh
                dat[2] -= a/2;
                if (!Rf.eof())
                    KK.push_back(dat);
            }
        }
    }
#endif
    f << "Input parameters: " <<
              "\n\tP/Pmax = " << force/4 <<
              "\n\tdelta = " << delta <<
              "\n\terr = " << err <<
              "\n\tmesh_h = " << mesh_h <<
              "\n\t b = " << b <<
              "\n\tdirectory = " << to_save << std::endl;

    bool rotate = false;
    bool shift = false;
    AniMesh am = generate_rectangle(a, b, mesh_h, 1);
    if (rotate){
        for (auto& lbl: am.vertex_label) lbl = 0;
        for (int i = 0; i < am.vertices.size() / 3; ++i){
            if (fabs(am.vertices[3*i+1] - b/2) < b*DBL_EPSILON) am.vertex_label[i] |= 2;
            if (fabs(am.vertices[3*i+1] + b/2) < b*DBL_EPSILON) am.vertex_label[i] |= 4;
            if (fabs(am.vertices[3*i+0] + a/2) < a*DBL_EPSILON || fabs(am.vertices[3*i+0] - a/2) < a*DBL_EPSILON) am.vertex_label[i] |= 1;
            swap(am.vertices[3*i+0], am.vertices[3*i+1]);
        }
        swap(a, b);
    }
    if (shift){
        std::array<double, 3> dv = {1, 8, 0};
        for (int i = 0; i < am.vertices.size() / 3; ++i){
            for (int l = 0; l < 3; ++l)
                am.vertices[3*i+l] += dv[l];
        }
    }
    static const int CLAMPED = 4, PRESSURED = 2, FREE = 1;
    Mesh m = convert_to_Mesh(am, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = "rectangle";
    obj.save(fpref + "_inited.vtk");
//#define RESIDUAL
//    //set right values to check residual
#ifdef RESIDUAL
    for (auto v: obj.m_mesh.vertices()){
        auto dat = g_Kproxy.getData(obj.m_x0[v][0]);
        obj.m_x[v] = obj.m_next_x[v] = Point(dat[2], obj.m_x0[v][1], dat[3]);
    }
#endif
    for (auto f: obj.m_mesh.faces()){
        auto f0 = *(obj.m_mesh.faces().begin());
        if (obj.m_normal[f][2] * obj.m_normal[f0][2] <= 0) {
            std::cout << "Wrong orientation!\n";
            exit(-1);
        }
    }
//    obj.name = "rectangle";
//    obj.save(fpref + "_inited.stl");
//    std::cout << "Saved in \"" + fpref + "_inited.stl\"" << std::endl;
//    exit(0);


    SimplePressureEdgeLoad load({0, 0, force*b},
                                [](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] == PRESSURED); },
                                [](const Object3D& obj, V_ind v) { return (obj.m_boundary[v] & PRESSURED); } );

//    SimpleNormalEdgeLoad load(force*b, PRESSURED);
    double lambda = E * nu / (1 - nu * nu), mu = E / (2 * (1 + nu));
    std::cout << "Elastic params: lambda = " << lambda << ", mu = " << mu << std::endl;

    Force svk;
    switch(model_type){
        case SVK_MODEL:
            svk = SVKirchhoffNoPoissonModel(H, lambda, mu, dir, false);
            break;
        case LIN_MODEL:
            svk = LinearElasticModel1(H, lambda, mu, dir,  false);
            break;
        default:
            throw std::runtime_error("Unknown model type");
    }
    svk.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
//    NeoGookModel svk(E/3, 5e-2, dir, true);
    BendingForce svkbf(svk.target<HyperElasticForceBase>()->f, false);
    svkbf.prepareJacobianFunction(false);
    map<E_ind, BendingForce::BoundaryData::Entry> boundary;
    for (auto e: obj.m_mesh.edges()){
        auto f = face_around(obj.m_mesh, e);
        if (f[0].second && f[1].second)
            continue;
        auto v = vert_around(obj.m_mesh, e);
        if ((obj.m_boundary[v[0]] & CLAMPED) && (obj.m_boundary[v[1]] & CLAMPED))
            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, {-1, 0, 0});
        else if (((obj.m_boundary[v[0]] & FREE) && (obj.m_boundary[v[1]] & FREE))
                || ((obj.m_boundary[v[0]] & PRESSURED) && (obj.m_boundary[v[1]] & PRESSURED)))
            boundary[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
    }

    thread* t = nullptr;
    if (view)
        t = new thread(start_debug_gui, argc, argv);

    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto id = w.addObject3D(move(obj), 1);

    svkbf.set_boundary_data(w.obj(id), boundary);
    w.addForce(std::move(svk));
    if (!no_bnd)
        w.addForce(std::move(svkbf));
    w.addForce(std::move(load), id);

    w.setForceApplier(StaticForceApplierP(delta)); //устойчивость неогука до 6.0е-6
    int freq = 20000;
    NewtonSolverStepAlgo nssa;

#ifdef USE_INMOST_SOLVERS
    nssa._solver.solver = LinearSolver(INMOST::Solver::INNER_MPTILUC);
#else
    nssa._solver.solver = std::move(LinearSolver("eigen"));
#endif
    nssa.newton_algo = [](NewtonSolverStepAlgo& ns, World* w){
        ns.updateMeshData(w);
        int N = ns.getNodfs(w);
        ns.assembleSystem(w, N);
        double tau = 0.3;//*1e-4;//0.3;
        static double resid_init = -1;
        if (ns.m_it == 0) {
            ns._solver.solver.SetParameterReal("absolute_tolerance", 1e-8);
            ns._solver.solver.SetParameterReal("relative_tolerance", 1e-20);
            double droptol = 2e-6;// /5;
            resid_init = ns.computeRhsNorm(1);
            std::cout << "\teps = " << resid_init << " drop_tol = " << droptol << std::endl;
            ns._solver.solver.SetParameterReal("drop_tolerance", droptol);
            ns._solver.solver.SetParameterReal("reuse_tolerance", 0.1 * droptol);
        } else {
            double eps = ns.computeRhsNorm(1) / resid_init;
            double droptol = 1e-4/2;// /4;//8e-4;//4e-5;//1e-6;
            if (eps > 100){
                tau = 0.3;//*1e-4;
                droptol = 5e-5/2;// /4;//1e-4;
            } else if (eps < 0.6){
                tau = 1;//*1e-3;
            } else {
                tau = 0.9;// * 1e-4;
            }
            std::cout << "\teps = " << eps << " drop_tol = " << droptol << std::endl;
            ns._solver.solver.SetParameterReal("drop_tolerance", droptol);
            ns._solver.solver.SetParameterReal("reuse_tolerance", 0.1 * droptol);
        }
        ns.Solve();
        ns.applySolution(w, tau);
        ns.m_it++;
        return 0;
    };

    w.setStepSimulationAlgo(std::move(nssa));
    freq = 1;

    V_ind wv;
    {
        Vector comp = CGAL::NULL_VECTOR;
        for (auto &v: w.obj(id).m_mesh.vertices())
            if (w.obj(id).m_boundary[v] == (PRESSURED | FREE))
                comp += (w.obj(id).m_x0[v] - CGAL::ORIGIN)/2;
        Point compP = CGAL::ORIGIN + comp;
        for (auto &v: w.obj(id).m_mesh.vertices())
            if ((w.obj(id).m_x0[v] - compP).squared_length() < 1e-7)
                wv = v;
    }

    w.obj(id).setDirichletCond([wv](Object3D& obj, const V_ind& v) -> unsigned char{
        if (obj.m_boundary[v] & CLAMPED) return 0;
//        if (obj.m_boundary[v] & FREE) return (1 << 0) | (1 << 2);
//        if (v == wv) return (1 << 0) | (1 << 2);
        return 7;
    });

#ifdef BEND_WATCH
    glob_k_it = -1;
    K_comp = w.obj(id).m_mesh.add_property_map<F_ind, std::array<double, 3*5>>("k_compare").first;
#endif
#ifndef RESIDUAL
    auto time = World3d::Timer();
    w.Simulation([&id, &time, &force, &err, &out = f, wv, &to_save, &fpref, &delta, &freq, &a, &b, &test_name, &mesh_h, &prefix](StepSimInfo& info, World* w)->bool{
        int it = info.it;
#ifdef BEND_WATCH
        glob_k_it++;
#endif
        int maxits = 200;
        static double resid_init;
        if (it == 1) resid_init = w->getWorldNorm("v:force", 1, 1, 1);
        if (it > 0 && it % freq == 0){
            double resid = w->getWorldNorm("v:force", 1, 1, 1);
            double eps = resid / resid_init;
            out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
            out << "watch { " << w->obj(id).m_x0[wv][0] - w->obj(id).m_x[wv][0] << " " << w->obj(id).m_x[wv][1] - w->obj(id).m_x0[wv][1] << " " << w->obj(id).m_x[wv][2] - w->obj(id).m_x0[wv][2] << " }" << std::endl;
            if (&out != &std::cout){
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << " ";
                std::cout << "watch { " << w->obj(id).m_x0[wv][0] - w->obj(id).m_x[wv][0] << " " << w->obj(id).m_x[wv][1] - w->obj(id).m_x0[wv][1] << " " << w->obj(id).m_x[wv][2] - w->obj(id).m_x0[wv][2] << " }" << std::endl;
            }
            if (eps < err || it >= maxits){
                out << ((eps < err) ? "Algorithm is converged: \n" : "Algorithm is achieved max iterations: \n");
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                out << "watch { " << w->obj(id).m_x0[wv][0] - w->obj(id).m_x[wv][0] << " " << w->obj(id).m_x[wv][1] - w->obj(id).m_x0[wv][1] << " " << w->obj(id).m_x[wv][2] - w->obj(id).m_x0[wv][2] << " }" << std::endl;
                if (&out != &std::cout) {
                    cout << ((eps < err) ? "Algorithm is converged: " : "Algorithm is achieved max iterations: ");
                    cout << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << " ";
                    cout << "watch { " << w->obj(id).m_x0[wv][0] - w->obj(id).m_x[wv][0] << " " << w->obj(id).m_x[wv][1] - w->obj(id).m_x0[wv][1] << " " << w->obj(id).m_x[wv][2] - w->obj(id).m_x0[wv][2] << " }" << std::endl;
                    string fname = fpref + "_converged.stl";
                    w->obj(id).save(fname);
                    string main_data = to_save + "/" + test_name + ".csv";
                    ofstream fdata(main_data, std::ios_base::app);
                    fdata << mesh_h << ", " << force / 4 << ", " << b << ", " << w->obj(id).m_x0[wv][0] - w->obj(id).m_x[wv][0] << ", "
                          << w->obj(id).m_x[wv][2] - w->obj(id).m_x0[wv][2]
                          << ", " << it << ", " << time.elapsed() << ", " + prefix + "\n";
                    fdata.close();
                }
                return true;
            }
            if (std::isnan(eps)){
                out << "Algorithm is diverged: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                out << "watch { " << w->obj(id).m_x0[wv][0] - w->obj(id).m_x[wv][0] << " " << w->obj(id).m_x[wv][1] - w->obj(id).m_x0[wv][1] << " " << w->obj(id).m_x[wv][2] - w->obj(id).m_x0[wv][2] << " }" << std::endl;
                if (&out != &std::cout) {
                    string fname = fpref + "_diverged.stl";
                    w->obj(id).save(fname);
                }
                return true;
            }
        }

        return false;
    });

    if (view) {
        t->join();
        delete t;
    }
#endif
#ifdef USE_INMOST_SOLVERS
#ifdef RESIDUAL
    auto res_tag = w.obj(id).m_mesh.add_property_map<V_ind, std::array<double, 3>>("v:residual").first;
    {
        NewtonSolverStepAlgo nssa;
        nssa._solver.solver = LinearSolver(INMOST::Solver::INNER_MPTILUC);
        nssa.newton_algo = [&res_tag, id](NewtonSolverStepAlgo &ns, World *w) {
            ns.updateMeshData(w);
            int N = ns.getNodfs(w);
            ns.assembleSystem(w, N);
            auto &solver = ns._solver;
            solver.x.resize(solver.rhs.size());
            for (auto v: w->obj(id).m_mesh.vertices()) {
                for (int i = 0; i < 3; ++i)
                    solver.x[v.idx() * 3 + i] = w->obj(id).m_x[v][i];
            }
            auto residual = LinearSolverInmost::SAXPY(1, ns._solver.sm, solver.x, -1, solver.rhs);
            for (auto v: w->obj(id).m_mesh.vertices()) {
                res_tag[v] = {residual[v.idx() * 3 + 0], residual[v.idx() * 3 + 1], residual[v.idx() * 3 + 2]};
            }

            ns.m_it++;
        };

        w.setStepSimulationAlgo(std::move(nssa));
        w.Simulation([](int it, World *w) -> bool {
#ifdef BEND_WATCH
            glob_k_it++;
#endif
            static int ret = 0;
            return (ret++ != 0);
        });
    }
#endif
#ifdef BEND_WATCH
    {
        auto& ob = w.obj(id);
        using namespace INMOST;
        INMOST::Mesh* mesh = new INMOST::Mesh();
        vector<HandleType> new_nodes(ob.m_mesh.num_vertices());
        vector<HandleType> new_faces(ob.m_mesh.num_faces());
        for (auto v: ob.m_mesh.vertices()){
            std::array<Storage::real,3> xyz{ob.m_x[v][0], ob.m_x[v][1], ob.m_x[v][2]};
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
        Tag k_exact = mesh->CreateTag("k_exact", DATA_REAL, NODE, FACE, 3);
        Tag k_comp = mesh->CreateTag("k_comp", DATA_REAL, NODE, FACE, 3);
        Tag dk = mesh->CreateTag("dk", DATA_REAL, NODE, FACE, 3);
#ifdef RESIDUAL
        Tag resid = mesh->CreateTag("resid", DATA_REAL, NODE, FACE, 3);
#endif
        for (auto n = mesh->BeginNode(); n != mesh->EndNode(); ++n){
#ifdef RESIDUAL
            for (int i = 0; i < 3; ++i) n->RealArray(resid)[i] = res_tag[V_ind(n->LocalID())][i];
#endif
            auto ff = face_around(ob.m_mesh, V_ind(n->LocalID()));
            std::array<double, 6> res = {0};
            for (auto f: ff){
                for (int i = 0; i < 6; ++i) res[i] += K_comp[f][i];
            }
            for (int i = 0; i < 6 && !ff.empty(); ++i) res[i] /= ff.size();
            for (int i = 0; i < 3; ++i) {
                n->RealArray(k_comp)[i] = res[i];
                n->RealArray(k_exact)[i] = res[i+3];
                n->RealArray(dk)[i] = res[i] - res[3+i];
            }
        }

        mesh->SaveVTK(fpref + "_kk.vtk");
        delete mesh;
    }
#endif
#endif

//    t.join();
    return 0;
}