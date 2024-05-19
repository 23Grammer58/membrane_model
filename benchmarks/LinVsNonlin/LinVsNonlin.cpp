//
// Created by alex on 23.01.2022.
//

#include "../BenchCommon.h"
#include "AVSim/AorticSimulator/Messurer.h"
#include "AVSim/Core/IO/Obj3dVTKSaver.h"
#include "AVSim/Core/NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverKinsol.h"
#include "AVSim/Core/NonLinearSolverCustom.h"
#include <Eigen/Dense>
#include "AVSim/Core/MeshGen/Helper.h"

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

static bool ObSaveTag(const Object3D& obj, string filename, string tag){
    ofstream ob(filename, ios::binary | ios::trunc);
    if (!ob) return false;
    auto m_x = obj.m_mesh.property_map<V_ind, Point>(tag).first;
    for (auto v: obj.m_mesh.vertices()){
        std::array<double, 3> P = {m_x[v][0], m_x[v][1], m_x[v][2]};
        ob.write(reinterpret_cast<const char *> ((P.data())), 3 * sizeof(double));
    }
    ob.close();
    return true;
}

static bool ObReadTag(Object3D& obj, string filename, string tag){
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
}

static map<E_ind, BendingForce::BoundaryData::Entry>
computeBndDataOnCilAligned(Object3D& leaf,
                           std::function<bool(Object3D&, E_ind, V_ind*/*[2]*/, V_ind)> choose_edge = [](Object3D& leaf, E_ind e, V_ind* v, V_ind vop){
                               const int FIXED = 2, FREE = 1|4|8;
                               return (leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED) && !(leaf.m_boundary[vop] & FREE);
                           },
                           double angle_align = 0)
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

static Mesh::Property_map<V_ind, Vector> addBdataTag(Object3D& leaf, const map<E_ind, BendingForce::BoundaryData::Entry>& bdata){
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

int LinVsNonlin(int argc, char* argv[]){
    double Ht = 0.3;
    double mu = 1e3; //kPa : 5e2, 1e3, 2e3, 3e3, 5e3
    double P = 90.0_mmHg / 1e3; //kPa
    int scenario = 0; // 0 - membrane, 1 - shape, 2 - shape + BC
    int material = 0; //0 - NHK, 1 - GNT, 2 - SVK, 3 - LIN
    int templ = 0;
    string init_config_prefix[] = {"_4", "_4_2", "_4_3", "_24", "_24_1"};
    double mesh_h[] = {0.5, 0.25, 0.125, 0.45, 0.45};
    static const char* material_types[] = {"NHK", "GNT", "SVK", "LIN"};
    static const double mu_def[] = {1.0e3/3, 2.0e3/3, 3.0e3/3, 5.0e3/3, 6.885e3/3, 9.0e3/3, 15.0e3/3};
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

    string  main_save = "../result/LinVsNonlin/",
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
            std::cout << "  -mt, --material     <Material type: 0 - NHK, 1 - GNT, 2 - SVK, 3 - LIN>" << "\n";
            std::cout << "  -tl, --template     <Which template: 0 - 0.5 mm, 1 - 0.25 mm, 2 - 0.125 mm, >2 - custom templates>" << "\n";
            std::cout << "  -tr, --target       <Directory in common directory to save results>" << "\n";
            std::cout << "  -md, --maindir      <Common directory to save results from different simulations>" << "\n";
            std::cout << "  -gd, --gendir       <Path to save generated codes>" << "\n";
            std::cout << "  -sp, --slvstop      <Solver stop criterion parameters: f_abs_err, maxnits>" << "\n";
            std::cout << "  -cs, --case         <Set default parameters for the case or define case by parameters>" << "\n";
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
        if (strcmp(argv[i], "-mt") == 0 || strcmp(argv[i], "--material") == 0) {
            try {
                if (i + 1 < argc) {
                    material = stoi(argv[++i]);
                    if (material < 0 || material > 3){
                        throw std::runtime_error("Wrong value of material, you should choose one of the following:\n\t0 - NHK\n\t1 - GNT\n\t2 - SVK\n\t3 - LIN\n");
                    }
                }
            } catch (std::exception& e){
                std::cout << "Waited number of material but error happens: '" << e.what() << "'\n"; --i;
            }
            continue;
        }
        if (strcmp(argv[i], "-tl") == 0 || strcmp(argv[i], "--template") == 0) {
            try {
                if (i + 1 < argc) {
                    templ = stoi(argv[++i]);
                    if (templ < 0 || templ > sizeof(init_config_prefix)/sizeof(string)){
                        throw std::runtime_error("Wrong value of template, you should choose one of the following:\n\t0 - 0.5mm\n\t1 - 0.25mm\n\t2 - 0.125mm\n>2 - custom leafs");
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
            bool read_next = i + 1 < argc && argv[i+1][0] != '-';
            try {
                int    mt = 0, //material: 0 - NHK, 1 - GNT, 2 - SVK, 3 - LIN
                       st = 0, //stiffness: 0 - 5e2, 1 - 1e3, 2 - 2e3, 3 - 3e3, 4 - 5e3
                       sc = 0; //scenario: 0 - membrane, 1 - shape + BC
                int nmt = 4, nst = sizeof(mu_def) / sizeof(double), nsc = 2;
                static const double* mu_arr = mu_def;//{ 5e2, 1e3, 2e3, 3e3, 5e3 };
                if (read_next) {
                    int cs = stoi(argv[++i]);
                    if (cs < 0 || cs >= nmt*nst*nsc){
                        throw std::runtime_error("Wrong value of case, case should be integer in [0, " + std::to_string(nmt*nst*nsc-1) + "]");
                    }
                    _case = cs;
                    sc = cs % nsc, mt = (cs / nsc) % nmt, st = (cs / nsc / nmt) % nst;
                    scenario = (sc > 0) ? 2 : 0;
                    material = mt;
                    mu = mu_arr[st];
                } else {
                    mt = material;
                    if (scenario == 2) sc = 1;
                    else if (scenario == 0) sc = 0;
                    else throw std::runtime_error("Can't define case by parameters");
                    st = -1;
                    for (int l = 0; l < nst && st < 0; ++l) if (abs(mu_arr[l] - mu) < 1e-2) st = l;
                    if (st < 0) throw std::runtime_error("Can't define case by parameters");
                    _case = sc + nsc * ( mt + nmt * ( st ) ) ;
                }
            } catch (std::exception& e){
                std::cout << "Waited number of case but error happens: '" << e.what() << "'\n";
                if (read_next) --i;
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
            "\n\tmaterial_type = " << material_types[material] <<
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
    leaf.name = "ModelLeaf";
    leaf.read("../result/ModelLeafBench/template/leaf_init" + init_config_prefix[templ] + ".txt");
    string init_aniso_crd = "v:initial";
    auto init_x = leaf.m_mesh.add_property_map<V_ind, Point>(init_aniso_crd).first;
    ObReadTag(leaf, "../result/ModelLeafBench/template/leaf_init_x0" + init_config_prefix[templ] + ".txt",
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
    Force Ef, bf;
    switch (material) {
//        {"NHK", "GNT", "SVK", "LIN"};
        case 0: {
            Ef = NeoGookModelOpt(mu, Ht);
            bf = BendingForce(NeoGookModel(mu, Ht, to_gen, regenerate).f, regenerate);
            break;
        }
        case 1: {
            Ef = GentModel(mu, Ht, 2.3, to_gen, regenerate);
            bf = BendingForce(Ef.target<HyperElasticForceBase>()->f, regenerate);
            Ef.target<HyperElasticForceBase>()->prepareJacobianFunction(regenerate); Ef.target<HyperElasticForceBase>()->f.setVerbosity(2);
            break;
        }
        case 2: {
            Ef = SVKirchhoffModel(Ht, 2*mu, mu, to_gen, regenerate ); // lambda = E * nu / (1 - nu * nu), mu = E / (2 * (1 + nu));
            bf = BendingForce(Ef.target<HyperElasticForceBase>()->f, regenerate);
            Ef.target<HyperElasticForceBase>()->prepareJacobianFunction(regenerate); Ef.target<HyperElasticForceBase>()->f.setVerbosity(2);
            break;
        }
        case 3: {
//            Ef = LinearElasticModel(Ht, 2*mu, mu);
            Ef = LinearElasticModel1(Ht, 2*mu, mu, to_gen, regenerate );
            Ef.target<HyperElasticForceBase>()->prepareJacobianFunction(regenerate); Ef.target<HyperElasticForceBase>()->f.setVerbosity(2);
            bf = BendingForce(Ef.target<HyperElasticForceBase>()->f, regenerate);
//            if (scenario != 0) throw std::runtime_error("Bending part for linear material is not implemented yet");
            break;
        }
        default: throw std::runtime_error("Faced wrong material number = " + to_string(material));
    }

    if (scenario != 0) { bf.target<BendingForce>()->prepareJacobianFunction(regenerate); bf.target<BendingForce>()->setVerbosity(2); }
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
        { //relative deformation
            auto& ob = w.obj(oids[0]);
            auto rel_stretch = ob.m_mesh.add_property_map<E_ind, double>("e:rel_stretch").first;
            for (auto e: ob.m_mesh.edges()) {
                auto v = vert_around(ob.m_mesh, e);
                rel_stretch[e] = sqrt( (ob.m_x[v[0]] - ob.m_x[v[1]]).squared_length() / (ob.m_x0[v[0]] - ob.m_x0[v[1]]).squared_length() );
            }
            auto rel_stretch_v = ob.m_mesh.add_property_map<V_ind, double>("v:rel_stretch").first;
            for (auto v: ob.m_mesh.vertices()){
                rel_stretch_v[v] = 0;
                auto ee = edge_around(ob.m_mesh, v);
                for (auto e: ee) rel_stretch_v[v] += rel_stretch[e];
                rel_stretch_v[v] /= ee.size();
            }
            auto H_F = ob.m_mesh.add_property_map<F_ind, double>("||H||_F").first;
            auto dR = ob.m_mesh.add_property_map<F_ind, double>("||R-R_0||_F").first;
            auto dU = ob.m_mesh.add_property_map<F_ind, double>("||U-U_0||_F").first;
            auto _Ap = set_Ap(ob);
            auto _Dv = set_D_vecs(ob);
            auto& x = ob.m_x;
            auto& x0 = ob.m_x0;
            for (auto f: ob.m_mesh.faces()){
                auto v = vert_around(ob.m_mesh, f);
                Eigen::Matrix<DReal, 3, 3> Q, P, D;
                Q <<    x [v[0]][0], x [v[1]][0], x [v[2]][0],
                        x [v[0]][1], x [v[1]][1], x [v[2]][1],
                        x [v[0]][2], x [v[1]][2], x [v[2]][2];
                P <<    x0[v[0]][0], x0[v[1]][0], x0[v[2]][0],
                        x0[v[0]][1], x0[v[1]][1], x0[v[2]][1],
                        x0[v[0]][2], x0[v[1]][2], x0[v[2]][2];
                D <<   _Dv[f][0][0],_Dv[f][1][0],_Dv[f][2][0],
                       _Dv[f][0][1],_Dv[f][1][1],_Dv[f][2][1],
                       _Dv[f][0][2],_Dv[f][1][2],_Dv[f][2][2];
                auto I = Eigen::Matrix<DReal, 3, 3>::Identity();
                using M3d = Eigen::Matrix<DReal, 3, 3>;
                M3d F0 = P*D.transpose(), F = Q*D.transpose();
                M3d H = F-F0;
                H_F[f] = H.norm();
                auto Fsvd = F.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV);
                auto R = Fsvd.matrixU()*Fsvd.singularValues().asDiagonal()*Fsvd.matrixU().adjoint();
                M3d U;
                {
                    auto lU = Fsvd.matrixU(), lV = Fsvd.matrixV();
                    U = lU.col(0)*lV.col(0).adjoint() + lU.col(1)*lV.col(1).adjoint();
                    if (lU.col(2).adjoint()*lV.col(2) > 0)
                        U += lU.col(2)*lV.col(2).adjoint();
                    else
                        U -= lU.col(2)*lV.col(2).adjoint();
                }
                auto F0svd = F0.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV);
                auto R0 = F0svd.matrixU()*F0svd.singularValues().asDiagonal()*F0svd.matrixU().adjoint();
                M3d U0;
                {
                    auto& U = U0;
                    auto lU = F0svd.matrixU(), lV = F0svd.matrixV();
                    U = lU.col(0)*lV.col(0).adjoint() + lU.col(1)*lV.col(1).adjoint();
                    if (lU.col(2).adjoint()*lV.col(2) > 0)
                        U += lU.col(2)*lV.col(2).adjoint();
                    else
                        U -= lU.col(2)*lV.col(2).adjoint();
                }
                dR[f] = (R-R0).norm();
                dU[f] = (U-U0).norm();
            }
        }
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
                        "scenario; Ht; P; mu; material; mesh_h; f_abs_err; maxnits; delta_init; sim_status\n";
        }
        else
            comtable.open(main_save + "comres.csv", std::ios_base::app);

        comtable << "\"" << to_save << "\"" << ";" << H << ";" << Hc << ";" << (areas[0].colArea + areas[1].colArea + areas[2].colArea)/3 << ";"
                 << bill.getBillowing(0) << ";" << l_free << ";"
                 << scenario << ";" << Ht << ";" << P << ";" << mu << ";" << material_types[material] << ";"
                 << mesh_h[templ] << ";" << f_err_lim << ";" << maxnit << ";" << delta_relax << ";" << (do_solve ? slvFlag : 2) << std::endl;

        comtable.close();
    }

    if (view){
        view_window->join();
        delete view_window;
    }

    return 0;
}

