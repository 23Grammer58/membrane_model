//
// Created by alex on 02.01.2021.
//

#include "Benchmarks.h"
#include "AorticSimulator/Messurer.h"

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

static bool readTag(Object3D& obj, string filename, string tag){
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

int BenchThreeLeaflet(int argc, char* argv[]){
//    thread t(start_debug_gui, argc, argv);
    string method = "NGK", Eval = "1";
    bool bend = true;
    for (int i = 1; i < argc; ++i){
        if (argv[i] == string("-h") || argv[i] == string("--help")){
            std::cout << "Help message: " << "\n";
            std::cout << " Command line options: " << "\n";
            std::cout << "  -m, --model   <Model formulation (ElasticModel, E_value, is_bend): {NGK, SVK, GNT}, {0.3, 1, 10}, {0, 1}>" << std::endl;
            exit(0);
        }
        if (argv[i] == string("-m") || argv[i] == string("--model")){
            if (i + 1 >= argc || argv[i+1][0] == '-') continue;
            method = argv[i+1];
            i++;
            if (i + 1 >= argc || argv[i+1][0] == '-') continue;
            Eval = argv[i+1];
            i++;
            if (i + 1 >= argc || argv[i+1][0] == '-') continue;
            bend = (argv[i+1] == string("1"));
        }
    }

    double P = 80 * 133.322;
    double R = 10;
    double mesh_h = 0.25;
    string to_save = "../result/bendInfluenceH/" + method + "-" + Eval + "/" + (bend ? "Bend/" : "NoBend/");
    string dir = "../generated";
    double margin = 0.05;
    double delta = 2.5e-7;
    double err = 0.076;//1e-5;
    int maxits = 3000000;
    int freq = 200;
    int saveFreq = bend ? 10000 : 5000;
    double H = 0.5, E = 1e6, nu = 0.5;
    if (Eval == "10") {
        //maxits = 600000;
        E = 1e7;
        delta = bend ? 2.5e-7/4 : 5e-7/10;
        err = 0.0036;
    }
    else if (Eval == "1") {
        //maxits = 100000;
        E = 1e6;
        delta = bend ? 2.5e-7 : 5e-7 / 2;
        err = 0.076;
    }
    else if (Eval == "0.1") {
        //maxits = 600000;
        E = 1e5;
        delta = bend ? 2.5e-7 : 5e-7;
        err = 0.076;
    }
    else if (Eval == "0.3") {
        //maxits = 200000;
        E = 3*1e5;
        delta = bend ? 2.5e-7/16 : 5e-7/2;
        err = 0.0036;
    }

    bool load_from_file = false;
    int ait = 0;
    path path = to_save;
    auto s = status(path);
    if (s.type() == file_type::directory){
        for (auto& p: directory_iterator(to_save)){
            if (p.status().type() == file_type::directory){
                string name = p.path().filename();
                auto pos = name.find_last_of('_');
                if (pos != name.npos && name.size() - pos > 1) {
                    int nit = stoi(name.substr(name.find_last_of('_') + 1));
                    if (nit > ait) ait = nit;
                }
            }
        }
        if (ait != 0) load_from_file = true;
    }

    string fname = to_save + "log.txt";
    std::ofstream f;
    if (load_from_file) f.open(fname, std::ios_base::app);
    else f.open(fname);
    f << "Input parameters: " <<
      "\n\tload_from_file = " << load_from_file <<
      "\n\tait = " << ait <<
      "\n\tForceType = " << method <<
      "\n\tE = " << E <<
      "\n\tP = " << P <<
      "\n\tR = " << R <<
      "\n\tH = " << H <<
      "\n\tnu = " << nu <<
      "\n\tmargin = " << margin <<
      "\n\terr = " << err <<
      "\n\tmaxits = " << maxits <<
      "\n\tfreq = " << freq <<
      "\n\tsaveFreq = " << saveFreq <<
      "\n\tbend = " << bend <<
      "\n\tdelta = " << delta <<
      "\n\tmesh_h = " << mesh_h <<
      "\n\tdirectory = " << to_save << std::endl;
    std::cout << "Input parameters: " <<
      "\n\tload_from_file = " << load_from_file <<
      "\n\tait = " << ait <<
      "\n\tForceType = " << method <<
      "\n\tE = " << E <<
      "\n\tP = " << P <<
      "\n\tR = " << R <<
      "\n\tH = " << H <<
      "\n\tnu = " << nu <<
      "\n\tmargin = " << margin <<
      "\n\terr = " << err <<
      "\n\tmaxits = " << maxits <<
      "\n\tfreq = " << freq <<
      "\n\tsaveFreq = " << saveFreq <<
      "\n\tbend = " << bend <<
      "\n\tdelta = " << delta <<
      "\n\tmesh_h = " << mesh_h <<
      "\n\tdirectory = " << to_save << std::endl;
    auto mt = World3d::Timer();

    std::array<Object3D, 3> obj;
    std::array<map<E_ind, BendingForce::BoundaryData::Entry>, 3> boundary;
    std::array<set<E_ind>, 3> free_b;
    static const int CLAMPED = 4, FREE = 1;

    if (load_from_file){
        string dir = to_save + "interim_" + to_string(ait) + "/";
        for (int i = 0; i < 3; ++i) {
            obj[i].read(dir + "leaf" + to_string(i) + ".txt");
            readTag(obj[i], dir + "leaf" + to_string(i) + "_m_x0.tag", "v:point");
            obj[i].name = "leaf" + to_string(i);
        }
    }
    else {
        std::array<AniMesh, 3> am;
        am[0] = generate_half_circle(R, mesh_h);
        if (bend) divideBorderElems(am[0]);
        am[1] = am[2] = am[0];
        double Rc = 6 * R / (2 * M_PI) * (1 + 1e-5);
        for (int i = 0; i < 3; ++i) {
            for (int n = 0; n < am[i].vertices.size() / 3; ++n) {
                double x0 = am[i].vertices[3 * n + 0], y0 = am[i].vertices[3 * n + 1];
                double phi = (x0 / R * M_PI + 2 * M_PI * i) / 3;
                am[i].vertices[3 * n + 0] = Rc * cos(phi);
                am[i].vertices[3 * n + 1] = Rc * sin(phi);
                am[i].vertices[3 * n + 2] = y0;
            }
        }
        f << "Generated mesh: " << mt.elapsed() << std::endl;

        for (int i = 0; i < 3; ++i) {
            obj[i].m_mesh = get_invert_Mesh(convert_to_Mesh(am[i], "v:boundary_lbl"));
            obj[i].reset_mesh();
            obj[i].name = "leaf" + to_string(i);
            obj[i].save(to_save + obj[i].name + "_init.stl");
        }
    }
    for (int i = 0; i < 3; ++i){
        for (auto e: obj[i].m_mesh.edges()) {
            auto f = face_around(obj[i].m_mesh, e);
            if (f[0].second && f[1].second)
                continue;
            auto v = vert_around(obj[i].m_mesh, e);
            if ((obj[i].m_boundary[v[0]] & CLAMPED) && (obj[i].m_boundary[v[1]] & CLAMPED)) {
                std::array<Point, 2> p = {obj[i].m_x0[v[0]], obj[i].m_x0[v[1]]};
                auto phi_branch = [i](double phi){ if (i == 2) phi += 2*M_PI; return phi; };
                std::array<double, 2> phi = {phi_branch(atan2(p[0].y(), p[0].x())), phi_branch(atan2(p[1].y(), p[1].x()))};
                std::array<Point_2, 2> x0 = {Point_2(R / M_PI * (3*phi[0] - 2 * M_PI * i), p[0].z()), Point_2(R / M_PI * (3*phi[1] - 2 * M_PI * i), p[1].z())};
                double _teta[2] = {atan2(x0[0].y(), x0[0].x()), atan2(x0[1].y(), x0[1].x())};
                for (int q = 0; q < 2; ++q) if (_teta[q] > 0) _teta[q] -= 2 * M_PI;
                double teta = (_teta[0] + _teta[1]) / 2;
                Point_2 n(cos(teta), sin(teta));
                double phim = (n[0] + 2 * i) * M_PI / 3;
                std::array<double, 3> norm = {-sin(phim) * n[0], cos(phim) * n[0], n[1]};
                boundary[i][e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, norm);
            }
            else if ((obj[i].m_boundary[v[0]] & FREE) && (obj[i].m_boundary[v[1]] & FREE)) {
                boundary[i][e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
                free_b[i].insert(e);
            }
        }
    }
    f << "Prepared objects: " << mt.elapsed() << std::endl;

    World w;
    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BulletCollisionManager>();
    w.setCollider(move(colMan));
//    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    std::array<ObjectID, 3> oid;
    for (int i = 0; i < 3; ++i){
        oid[i] = w.addObject3D(move(obj[i]), 1);
        Force Pr = SimplePressureLoad(P);
        double mu = E/3, Jm = 2.3, Lambda = E * nu / (1 - nu * nu), Mu = E / (2 * (1 + nu));
        Force elastic;
        if (method == "NGK")
            elastic = NeoGookModel(mu, H, dir, false);//i==0);
        else if (method == "SVK")
            elastic = SVKirchhoffModel(H, Lambda, Mu, dir, false);//i==0);
        else if (method == "GNT")
            elastic = GentModel(mu, H, Jm, dir, false);//i==0);
        else if (!bend && method == "MSM")
            elastic = SimpleMassSpringModel(E, H);
        else throw std::runtime_error("Force case \"" + method + "\" - is not implemented");
        Force bending;
        if (bend) {
            bending = BendingForce(elastic.target<HyperElasticForceBase>()->f, false);//i == 0);
            bending.target<BendingForce>()->set_boundary_data(w.obj(oid[i]), boundary[i]);
        }
        w.addForce(std::move(Pr), oid[i]);
        w.addForce(std::move(elastic), oid[i]);
        if (bend)
            w.addForce(std::move(bending), oid[i]);
        reinterpret_cast<BulletCollisionManager*>(w.getCollider())->set_margin(oid[i], margin);

        w.obj(oid[i]).setDirichletCond([](Object3D& obj, const V_ind& v) -> unsigned char{
            if (obj.m_boundary[v] & CLAMPED) return 0;
            return 7;
        });
    }

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
    static auto full_save = [](const Object3D& obj, string directory, string postfix){
        string new_dir = directory + postfix + "/";
        create_directory(new_dir);
        string fname = new_dir +  obj.name;
        obj.save(fname  + ".stl");
        obj.save(fname  + ".txt");
        saveTag(obj, fname  + "_m_x0.tag", "v:point");
    };

    f << "Generated forces and prepared world: " << mt.elapsed() <<
        "\nStart simulation" << std::endl;

    auto time = World3d::Timer();
//    for (int it  = 0; it < ait; ++it) {
//        if (bend && fabs(E - 1e6) < 1) {
//            if (it == 2400) delta /= 2;
//            if (it == 2400 + 2400 * 2) delta /= 2;
//            if (it == 45000) delta /= 2;
//            if (method == "GNT" && it == 60000) delta /= 2;
//            if (method == "GNT" && it == 155000) delta /= 2;
//        } else if (bend && fabs(E - 1e7) < 1) {
//            if (it == 3600) delta /= 2;
//        } else if (bend && fabs(E - 1e5) < 1) {
////            if (it == 1600) delta /= 2;
////            if (it == 2800) delta /= 2;
////            if (it == 15000) delta /= 2;
////            if (it == 22000) delta /= 2;
//        }
//    }

    StaticForceApplierP sfa(delta);
    sfa.addDamper(Damper(0.01*H, 0.8, 100, 100000).setReInitDxParams(1e7, 0.4));
    w.setForceApplier(sfa);
    w.Simulation([&oid, &time, &err, &out = f, &delta, &freq, &free_b, maxits, to_save, saveFreq, bend, E, ait](StepSimInfo& info, World* w)->bool{
        int it = info.it;
        it += ait;
        if (info.status < 0) {
            out << "Unrecoverable error in residual computation at it = " << it <<  std::endl;
            std::cout << "Unrecoverable error in residual computationat it = " << it << std::endl;
            return true;
        }
        auto comp_lengths = [&oid, &free_b](World* w) -> std::array<double, 3>{
            std::array<double, 3> res = {0, 0, 0};
            for (int i = 0; i < 3; ++i){
                auto& obj = w->obj(oid[i]);
                for (auto e: free_b[i]) {
                    auto v = vert_around(obj.m_mesh, e);
                    res[i] += sqrt((obj.m_x[v[1]] - obj.m_x[v[0]]).squared_length());
                }
            }
            return res;
        };
        static double resid_init = 6.26503e6;
//        if (ait > 0) {
//            std::cout << "Warning: unknown resid_init" << std::endl;
//            resid_init = 6.26503e6;//6.12804e6;
//        }
        static std::array<double, 3> len0 = {0, 0, 0};
        auto len1 = comp_lengths(w);
        std::array<double, 3> dll = {len1[0] - len0[0], len1[1] - len0[1], len1[2] - len0[2]};
        len0 = len1;
        double dl = fabs(dll[0]) + fabs(dll[1]) + fabs(dll[2]);
        if (it == 1) resid_init = w->getWorldNorm("v:force", 1, 1, 1);
//        if (bend && fabs(E - 1e6) < 1) {
//            if (it == 2400) delta /= 2;
//            if (it == 2400 + 2400 * 2) delta /= 2;
//          if (it == 45000) delta /= 2;
//        }
//        else if (bend && fabs(E - 1e7) < 1) {
//            if (it == 3600) delta /= 2;
//        }
//        else if (bend && fabs(E - 1e5) < 1) {
////            if (it == 1600) delta /= 2;
////            if (it == 2800) delta /= 2;
////            if (it == 15000) delta /= 2;
////            if (it == 22000) delta /= 2;
//        }
        if (it > ait && it % saveFreq == 0){
            for (int i = 0; i < 3; ++i){ full_save(w->obj(oid[i]), to_save, "interim_" + to_string(it)); }
        }
        if (it > ait && it % freq == 0){
            double resid = w->getWorldNorm("v:force", 1, 1, 1);
            double eps = resid / resid_init;
            out << "it " << it << ": eps = " << eps << " abs = " << resid << " dl = " << dl <<
                " len = {" << len1[0] << ", " << len1[1] << ", " << len1[2] << "}" <<" time = " << time.elapsed() << std::endl;
            if (&out != &std::cout){
                std::cout << "it " << it << ": eps = " << eps << " abs = " << resid << " dl = " << dl <<
                          " len = {" << len1[0] << ", " << len1[1] << ", " << len1[2] << "}" <<" time = " << time.elapsed() << "\n";
            }
            if (eps < err || it >= maxits /*|| dl < err*/){
                out << ((eps < err || dl < err) ? "Algorithm is converged: \n" : "Algorithm is achieved max iterations: \n");
                out << "it " << it << ": eps = " << eps << " abs = " << resid << " dl = " << dl <<
                      " len = {" << len1[0] << ", " << len1[1] << ", " << len1[2] << "}" <<" time = " << time.elapsed() << std::endl;
                if (&out != &std::cout) {
                    cout << ((eps < err|| dl < err) ? "Algorithm is converged: " : "Algorithm is achieved max iterations: ");
                    cout << "it " << it << ": eps = " << eps << " abs = " << resid << " dl = " << dl <<
                       " len = {" << len1[0] << ", " << len1[1] << ", " << len1[2] << "}" <<" time = " << time.elapsed() << "\n";
                    for (int i = 0; i < 3; ++i){ full_save(w->obj(oid[i]), to_save, "converged"); }
                }
                return true;
            }
            if (std::isnan(eps)){
                out << "Algorithm is diverged: \n";
                out << "it " << it << ": eps = " << eps << " abs = " << resid <<
                    " len = {" << len1[0] << ", " << len1[1] << ", " << len1[2] << "}" <<" time = " << time.elapsed() << std::endl;
                if (&out != &std::cout) {
                    cout << "Algorithm is diverged: " << std::endl;
                    for (int i = 0; i < 3; ++i){ full_save(w->obj(oid[i]), to_save, "diverged"); }
                }
                return true;
            }
        }

        return false;
    });

//    t.join();
    return 0;
}

int processThreeLeaflet(int argc, char* argv[]){
    vector<string> methods = {"SVK"};//{"NGK"};//{"GNT", "NGK", "SVK"};
    vector<string> Evals = {"0.3"};//{"1"};//{"0.3", "1", "10"};
    vector<bool> bends = {true};//{true};//{true, false};
    for (auto method: methods)
    for (auto Eval: Evals)
    for (auto bend: bends)
    {
        string to_save = "../result/bendInfluenceH/process/" + method + Eval + (bend ? "b" : "nb") + "/";
        create_directory(to_save);
        string to_read =
                "../result/bendInfluenceH/" + method + "-" + Eval + "/" + (bend ? "Bend/" : "NoBend/") + "converged/";
        Object3D dup;
        std::array<Object3D, 3> cusps;
        for (int i = 0; i < 3; ++i) {
            cusps[i].read(to_read + "leaf" + to_string(i) + ".txt");
            readTag(cusps[i], to_read + "leaf" + to_string(i) + "_m_x0.tag", "v:point");
            cusps[i].name = "leaf" + to_string(i);
            //convert to flat state
            for (auto v: cusps[i].m_mesh.vertices()) {
                auto phi_branch = [i](double phi) {
                    if (i == 2 && phi <= 0) phi += 2 * M_PI;
                    return phi;
                };
                Point p = cusps[i].m_x0[v];
                double phi = phi_branch(atan2(p.y(), p.x()));
                double R = 10;
                Point x0(R / M_PI * (3 * phi - 2 * M_PI * i), p.z(), 0);
                cusps[i].m_x0[v] = x0;
            }
        }
        Messurer mes(dup, cusps[0], cusps[1], cusps[2], Vector{0, 0, 1});
        mes.setMargin(0.2);
//        std::cout << "#F = " << cusps[0].m_mesh.num_faces() << " #E = " << cusps[0].m_mesh.num_edges() << " #V = "
//                  << cusps[0].m_mesh.num_vertices() << std::endl;
        std::cout << "method = " << method << " Eval = " << Eval << " bend = " << bend << std::endl;
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
            mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_init.stl", "v:point");
            mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_init_col.stl", "v:point",
                                colfilter);
            mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_init_nocol.stl", "v:point",
                                nocolfilter);
            mes.saveObjectShape(mes.m_colission_shapes[i], to_save + cusps[i].name + "_mes_real.stl", "v:x");
        }
        mes.computeHalfCoaptScans();
        auto distr = mes.computeHalfCoaptDistrib(500);
        mes.uniteCollidingCuspHalves();
        auto fdistr = mes.computeCoaptDistribution(1000);
        std::cout << "Hc = " << mes.computeHc() << " H = " << distr.getCoaptH() << " Closed = " << mes.isValveClosed() <<   std::endl;
        mes.computeBillowing().print();
        auto areas = mes.computeColArea();
        std::cout << "Coapt areas: \n\t0:" << areas[0] << "\n\t1:" << areas[1] << "\n\t2:" << areas[2] << "\n";
//        Messurer::saveHalfCoaptDistribCSV(distr, "../result/MeassureDebug/half_distrib.csv");
//        Messurer::saveCoaptDistribCSV(fdistr, "../result/MeassureDebug/full_distrib.csv");
        Messurer::saveHalfCoaptDistribCSV(distr, to_save + "half_distrib.csv");
    }

    return 0;
}

