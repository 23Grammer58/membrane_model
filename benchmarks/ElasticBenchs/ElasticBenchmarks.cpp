//
// Created by alex on 01.10.2020.
//

#include "../BenchCommon.h"
#include <gsl/gsl_poly.h>
#include "nlohmann/json.hpp"

#if __cplusplus >= 201703L
#include <filesystem>
#else
#include <experimental/filesystem>
#endif

struct Lendl{};
struct Logger{
    ostream& out = cout;
    ofstream log;
    bool log_b = false;
    bool verb = true;

    void open( const std::string &filename,
               ios_base::openmode mode = ios_base::out ){
        log.open(filename, mode);
        log << scientific;
        log.precision(14);
        log_b = true;
    }
    void close(){
        if (log_b)
            log.close();
        log_b = false;
    }
    template <typename T>
    Logger& operator<< (const T& value){
        if (log_b)
            log << value;
        if (verb)
            out << value;
        return *this;
    }
    template <typename T>
    void to_log (const T& value){
        if (log_b)
            log << value;
    }
    Logger& operator<< (const Lendl& value){
        if (log_b)
            log << endl;
        if (verb)
            out << endl;
        return *this;
    }
};
Lendl lendl;
Logger Log;

enum ElasticModels{
    EMOD_MSM = 1,
    EMOD_SVK,
    EMOD_NGK,
    EMOD_GNT,
    EMOD_MNY
};

struct InputParams{
    using json = nlohmann::json;

    double R0 = 10;
    double l = 10;
    double phi = M_PI_2;
    double h_mesh = 0.4;
    int MeshType = 1;
    ElasticModels EModelType = EMOD_GNT;
    double thickness = 0.5;
    double eprms[10] = {1000000, 2.3, 0, 0};
    double rel_allow_shift = 0.02 * 4;
    double rel_max_shift = 0.05 * 4;
    double P = 80; //mm Hg
    double delta = 5e-7;
    double err = 5e-4;
    double maxits = 100000;
    int compute_print_freq = 200;
    string dir = "../bin/";
    string prefix = "gent_";
    string logfile = "logfile.txt";
    int scenario = 0;
    json init_data;
private:
    static json get_json(string filename){
        std::ifstream js;
        js.open(filename);
        if (! js.good()) {
            throw std::runtime_error("Unsuccessful attempt to open file \"" + filename + "\"");
        }
        json conf;
        js >> conf;
        js.close();

        return std::move(conf);
    }
    static ElasticModels getModel(string ModelName){
        transform(ModelName.begin(), ModelName.end(), ModelName.begin(), ::tolower);
        map<string, ElasticModels> models = {
                {"msm", EMOD_MSM},
                {"mass-spring", EMOD_MSM},
                {"mass-spring model", EMOD_MSM},
                {"trqs", EMOD_SVK},
                {"st venan", EMOD_SVK},
                {"neogook", EMOD_NGK},
                {"neogokean", EMOD_NGK},
                {"neo-gookean", EMOD_NGK},
                {"reinforcing", EMOD_GNT},
                {"gent", EMOD_GNT},
                {"compressible neogook", EMOD_GNT},
                {"mooney", EMOD_MNY}
        };
        if (!models.count(ModelName)) return EMOD_MSM;
        else return models[ModelName];
    }
    static string getModelName(ElasticModels m){
        static map<ElasticModels, string> names = {
                {EMOD_MSM, "msm"},
                {EMOD_SVK, "st-venan"},
                {EMOD_NGK, "neogook"},
                {EMOD_GNT, "gent"},
                {EMOD_MNY, "mooney"}};
        return names[m];
    }
    void processConfigFile(const char* filename){
        json conf = get_json(filename);
        R0 = conf["R"];
        l = conf["l"];
        double _phi = conf["phi"];
        phi = _phi / 180 * M_PI;
        h_mesh = conf["mesh_step"];
        MeshType = conf["mesh_type"];
        EModelType = getModel(conf["elastic_model"]);
        prefix = getModelName(EModelType);
        thickness = conf["thickness"];
        vector<double> v = conf["model_params"].get<vector<double >>();
        for (int i = 0; i < 10 && i < v.size(); ++i)
            eprms[i] = v[i];
        rel_allow_shift = conf["rel_allow_shift"];
        rel_max_shift = conf["rel_max_shift"];
        P = conf["P"];
        delta = conf["delta"];
        err = conf["velocity_decay"];
        maxits = conf["max_count_of_iterations"];
        compute_print_freq = conf["frequency_of_printing"];
        dir = conf["save_directory"];
        if (conf.count("file_prefix"))
            prefix += conf["file_prefix"];
        logfile = conf["logfile_name"];
        if (conf.count("scenario"))
            scenario = conf["scenario"];
        init_data = std::move(conf);
    }
public:
    static InputParams processMainArgs(int argc, char* argv[]){
        InputParams ip = InputParams();
        int i = 0;
//        if (argc == 1) goto helpMessage;
        for (i = 1; i < argc; i++) {
            //Print help message and exit
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                helpMessage:
                std::cout << "Help message: " << std::endl;
                std::cout << "Command line options: " << std::endl;
                std::cout << "-c, --config <Main configuration parameters file name>" << std::endl;
                exit(0);
            }
            if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--config") == 0) {
                std::cout << "Main configuration file found: " << argv[i + 1] << std::endl;
                ip.processConfigFile(argv[i + 1]);
                i++;
                continue;
            }
        }
        return ip;
    }
};

static double evaluate_radius(int meshtype, Object3D& obj, double& max_div, double& mid_div);
static double evaluate_radius(int meshtype, Object3D& obj);
static double predicted_radius(double P /*Hg mm*/, double R0 /*mm*/, double T /*Pa*/, ElasticModels mod, int mesh_type,  double* prms);
static function<int(Object3D& )> getForceApplier(InputParams& prms);
static function<bool(Object3D& obj, const V_ind& v)> getMovable(InputParams& prms);
static int generate_config_files(string dir);

int ElasticBenchmark(int argc, char* argv[]){
    thread t(start_debug_gui, argc, argv);

    InputParams params = InputParams::processMainArgs(argc, argv);
    {
#if __cplusplus >= 201703L
        using namespace std::filesystem;
#else
        using namespace std::experimental::filesystem;
#endif
        path path = params.dir;
        auto s = status(path);
        if (!status_known(s) || s.type() != file_type::directory)
            create_directory(path);
    }
    generate_config_files(params.dir);
    if (params.logfile != "") Log.open(params.dir + params.prefix + params.logfile);
    Log.to_log(params.init_data); Log << lendl;

    AniMesh mesh;
    string mname = "";
    if (params.MeshType == 0) {
        mesh = generate_eigth_sphere_part(params.R0, params.h_mesh);
        mname = "sph";
    }
    else if (params.MeshType == 1) {
        mesh = generate_cilinder_part(params.R0, params.l, params.phi, params.h_mesh);
        mname = "cil";
    }
    else if (params.MeshType == 2) {
        mesh = generate_circle(params.R0, params.h_mesh);
        mname = "circle";
    }
    else if (params.MeshType == 3) {
        mesh = generate_annulus(params.l, params.R0, params.h_mesh);
        mname = "annulus";
    }
    else {
        Log << "ERROR: MeshType = " << params.MeshType << " is not defined" << lendl;
        return -1;
    }

    Mesh m = convert_to_Mesh(mesh, "v:boundary_lbl");
    Object3D obj{m};
    obj.name = mname;
    string gendir = "../generated";
    obj.set_is_movable_checker(getMovable(params));

    Force Fe;
    if (params.scenario <= 0) {
        switch (params.EModelType) {
            case EMOD_MSM:
                Fe = SimpleMassSpringModel(params.eprms[0], params.thickness);
                break;
            case EMOD_SVK:
                Fe = SVKirchhoffModel(params.thickness, params.eprms[0], params.eprms[1],  gendir);
                break;
            case EMOD_NGK:
                Fe = NeoGookModel(params.eprms[0]/3, params.thickness, gendir);
                break;
            case EMOD_GNT:
                Fe = GentModel(params.eprms[0]/3, params.thickness, params.eprms[1],  gendir);
                break;
            case EMOD_MNY:
                Fe = MooneyRivlinModel(params.thickness, params.eprms[0], params.eprms[1],  gendir);
                break;
            default:
                std::runtime_error("It's not implemented");
        }
    } else std::runtime_error("It's not implemented");

    double Pf = 0;
    SimplePressureLoad Fp([&Pf](Object3D&, F_ind){ return Pf; });

    World w;
    w.setRenderer(std::make_unique<World3d::DefaultRenderer>());
    auto id = w.addObject3D(move(obj), 1);

    w.addForce(std::move(Fe));
    w.addForce(std::move(Fp));
    w.setForceApplier(getForceApplier(params));
    auto time = World3d::Timer();
    auto stopCondition = [maxits = params.maxits, err = params.err, freq = params.compute_print_freq, &time](StepSimInfo& info, World* w)->bool{
        static double resid_init;
        int it = info.it;
        if (it == 1) resid_init = w->getWorldShift();
        if (it > 0 && it % freq == 0) {
            double resid = w->getWorldShift();
            double eps = resid / resid_init;
            Log << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
            if (eps < err){
                Log << "Algorithm is converged: \n";
                Log << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                return true;
            }
            if (eps > 3 || std::isnan(eps) || it > maxits){
                if (eps > 3 || std::isnan(eps))
                    Log << "Algorithm is diverged: \n";
                else if (it > maxits)
                    Log << "Algorithm reached maximum of iterations = " << maxits << ": \n";
                Log << "it " << it << ": eps = " << eps << " abs = " << resid << " time = " << time.elapsed() << "\n";
                return true;
            }
        }
        return false;
    };

    if (params.MeshType == 0 || params.MeshType == 1) {
        ofstream for_R(params.dir + "/res/" + params.prefix + obj.name + "_R(P)_ext_" + to_string(params.scenario) + ".csv", std::ios::trunc);
        double max_div = 0, mid_div = 0;
        double R_init = evaluate_radius(params.MeshType, w.obj(id), max_div, mid_div);
        for_R << 0 << ";" << R_init << ";" << mid_div << ";" << max_div << ";" << params.R0 << endl;
        int st = 1;
        for (int i = st, count = ((int)(round(params.P) + 0.1))*2; i <= count; ++i){
            double P = params.P * i / count;
            Pf = -133.3 * P;
            double R_exact = predicted_radius(P, params.R0, params.thickness, params.EModelType, params.MeshType,
                                              params.eprms);
            double kk = 133.3 * P * params.R0 / (params.eprms[0] / 3 * params.thickness);
            Log << "P * R / ( mu * H) = " << kk << lendl;
            w.Simulation(stopCondition);
            double R_model = evaluate_radius(params.MeshType, w.obj(id), max_div, mid_div);
            for_R << kk << ";" << R_model << ";" << mid_div << ";" << max_div << ";" << R_exact << endl;
            Log << "Evaluate radius = " << R_model << lendl;
            Log << "Exact radius = " << R_exact << lendl;
            Log << "|R_exact - R_eval| = " << fabs(R_exact - R_model) << lendl;
        }
        for_R.close();
    }

    return 0;
}

static function<bool(Object3D& obj, const V_ind& v)> getMovable(InputParams& prms){
    if (prms.MeshType == 0 || prms.MeshType == 1){
        return [](Object3D& obj, const V_ind& v) { return true; };
    } else if (prms.MeshType == 2){
        return [](Object3D& obj, const V_ind& v) { return !(obj.m_boundary[v] % 2); };
    } else if (prms.MeshType == 3){
        return [](Object3D& obj, const V_ind& v) { return !(obj.m_boundary[v] & 3); };
    } else
        return [](Object3D& obj, const V_ind& v) { return !(obj.m_boundary[v] % 2); };
}

static function<int(Object3D& )> getForceApplier(InputParams& prms){
    if (prms.MeshType == 0){
        return [&_delta = prms.delta](Object3D& obj) -> int{
            for (auto i: obj.m_mesh.vertices()){
                if (obj._is_movable(obj, i)){
                    obj.m_next_x[i] = obj.m_x[i] + _delta * obj.m_F[i];
                    int mask = obj.m_boundary[i];
                    obj.m_next_x[i] = Point((mask & 1) ? 0 : obj.m_next_x[i][0],
                                            (mask & 2) ? 0 : obj.m_next_x[i][1],
                                            (mask & 4) ? 0 : obj.m_next_x[i][2]);
                }
            }
            return 0;
        };
    }else if (prms.MeshType == 1){
        return [&_delta = prms.delta, phi = prms.phi](Object3D& obj) -> int{
            for (auto i: obj.m_mesh.vertices()){
                if (obj._is_movable(obj, i)){
                    obj.m_next_x[i] = obj.m_x[i] + _delta * obj.m_F[i];
                    int mask = obj.m_boundary[i];
                    if (mask) {
                        obj.m_next_x[i] = Point(obj.m_next_x[i][0],
                                                (mask & 1) ? obj.m_x0[i][1] : obj.m_next_x[i][1],
                                                (mask & 2) ? obj.m_x0[i][2] : obj.m_next_x[i][2]);
                    }
                    if (mask & 4) {
                        Vector n = Vector(-sin(phi), cos(phi), 0);
                        obj.m_next_x[i] -= ((obj.m_next_x[i] - obj.m_x0[i]) * n) * n;
                    }
                }
            }
            return 0;
        };
    }else if (prms.MeshType == 2){
        return [&_delta = prms.delta](Object3D& obj) -> int{
            for (auto i: obj.m_mesh.vertices()){
                if (obj._is_movable(obj, i)){
                    obj.m_next_x[i] = obj.m_x[i] + _delta * obj.m_F[i];
                    int mask = obj.m_boundary[i];
                    if (mask == 6){
                        obj.m_next_x[i] = Point(0, 0, obj.m_next_x[i][2]);
                    }
                    mask /= 2;
                    if (mask > 0) {
                        obj.m_next_x[i] = Point(obj.m_next_x[i][0],
                                                obj.m_x0[i][1],
                                                obj.m_next_x[i][2]);
                    }
                }
            }
            return 0;
        };
    }else if (prms.MeshType == 3){
        return [&_delta = prms.delta](Object3D& obj) -> int{
            for (auto i: obj.m_mesh.vertices()){
                if (obj._is_movable(obj, i)){
                    obj.m_next_x[i] = obj.m_x[i] + _delta * obj.m_F[i];
                    int mask = obj.m_boundary[i];
                    mask /= 4;
                    if (mask > 0) {
                        obj.m_next_x[i] = Point(obj.m_next_x[i][0],
                                                obj.m_x0[i][1],
                                                obj.m_next_x[i][2]);
                    }
                }
            }
            return 0;
        };
    } else return StaticForceApplier(prms.delta);
}

static double solve_poly_equation(const int n, double* a){
    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (n);
    double z[2 * (n - 1)];
    gsl_poly_complex_solve (a, n, w, z);
    gsl_poly_complex_workspace_free (w);
    int i = 0;
    int real = 0;
    double rreal = z[1];
    double eps = 1.0e-7;
    int solved = 0;
    for (i = 0; i < n - 1; ++i)
        if ((0 < z[2*i] && z[2*i] < 1 + eps)) {
            real = i;
            rreal = z[2*i + 1];
            if (fabs(rreal) < eps)
                solved = 1;
            break;
        }

    for (; i < n - 1; ++i)
        if ((0 < z[2*i] && z[2*i] < 1 + eps) && (fabs(z[2*i+1]) < fabs(rreal)) && (fabs(z[2*i+1]) < eps)){
            real = i;
            rreal = z[2*i + 1];
            solved = 2;
        }
    if (!solved || fabs(z[real * 2 + 1]) > eps)
#ifndef _RADIUS_TESTER
        throw std::runtime_error("Failed to solve the equation");
#else
    return 0;
#endif
    return (z[2*real] < 1) ? z[2*real] : 1;
}

double pred_lambda_neogook_cilinder(double pr_hmu){
    if (pr_hmu >= 1)
#ifndef _RADIUS_TESTER
        throw std::runtime_error("Failed to solve the equation, cilinder will grow to infininity");
#else
    return 0.0;
#endif
    return 1 / pow(1 - pr_hmu, 0.25);
}

static double pred_lambda_gent_sphere(double pr_hmu, double Jm){
    double k = pr_hmu / 2;
    const int n = 10;
    double a[n] = {2*k, 0, -k*(Jm + 3), Jm, 0, 0, k, 0, 0, -Jm};
    return 1/solve_poly_equation(n, a);
}

static double pred_lambda_gent_cilinder(double pr_hmu, double Jm){
    double k = pr_hmu;
    if (fabs(k) < DBL_EPSILON) return 1.0;
    double Jm_k = Jm / k;
    double a = Jm_k - (Jm + 2), b = 1, c = -Jm_k;
    array<double, 3> x;
    int n = gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
    sort(x.begin(), x.begin() + n);
    double z = 0; bool finded = false;
    for (int i = 0; i < n; ++i)
        if (x[i] > 1) {
            z = x[i];
            finded = true;
            break;
        }
    if (!finded)
#ifndef _RADIUS_TESTER
        throw std::runtime_error("Failed to solve the equation, there are no root > 1");
#else
    return 0.0;
#endif

    return sqrt(z);
}

static double pred_lambda_neogook_sphere(double pr_hmu){
    double k = pr_hmu / 2;
    const int n = 8;
    double a[n] = { -k, 1, 0, 0, 0, 0, 0, -1 };

    return 1/solve_poly_equation(n, a);
}

static double predicted_radius(double P /*Hg mm*/, double R0 /*mm*/, double T /*Pa*/, ElasticModels mod, int mesh_type,  double* prms){
    double mu = prms[0] / 3;
    P *= 133.3;
    double k = P * R0 / (mu * T);
    if (mod == EMOD_NGK && mesh_type == 0){
        return R0 * pred_lambda_neogook_sphere(k);
    }
    if (mod == EMOD_NGK && mesh_type == 1){
        return R0 * pred_lambda_neogook_cilinder(k);
    }
    if (mod == EMOD_GNT && mesh_type == 0){
        return R0 * pred_lambda_gent_sphere(k, prms[1]);
    }
    if (mod == EMOD_GNT && mesh_type == 1){
        return R0 * pred_lambda_gent_cilinder(k, prms[1]);
    }
    return -1e20;
}

static double test_spheriable(Object3D& obj, double& _max_div, double& _mid_div){
    double r = 0;
    for (auto v: obj.m_mesh.vertices())
        r += sqrt((obj.m_x[v] - CGAL::ORIGIN).squared_length());
    r /= obj.m_mesh.num_vertices();
    double mid_div = 0, max_div = 0;
    for (auto v: obj.m_mesh.vertices()){
        double div = fabs(sqrt((obj.m_x[v] - CGAL::ORIGIN).squared_length()) - r);
        mid_div += div;
        max_div = (max_div >= div) ? max_div : div;
    }
    mid_div /= obj.m_mesh.num_vertices();
    Log << "Evaluate radius = " << r << lendl;
    Log << "Maximal deviation from radius is " << max_div << lendl;
    Log << "Middle deviation from radius is " << mid_div << lendl;

    _max_div = max_div;
    _mid_div = mid_div;

    return r;
}

static double test_cilindrable(Object3D& obj, double& _max_div, double& _mid_div){
    double r = 0;
    auto length = [](const Point& p){
        Vector v = p - CGAL::ORIGIN;
        return sqrt(v.x() * v.x() + v.y() * v.y());
    };
    for (auto v: obj.m_mesh.vertices())
        r += length(obj.m_x[v]);
    r /= obj.m_mesh.num_vertices();
    double mid_div = 0, max_div = 0;
    for (auto v: obj.m_mesh.vertices())
    {
        double div = fabs(length(obj.m_x[v]) - r);
        mid_div += div;
        max_div = (max_div >= div) ? max_div : div;
    }
    mid_div /= obj.m_mesh.num_vertices();
    Log << "Evaluate radius = " << r << lendl;
    Log << "Maximal deviation from radius is " << max_div << lendl;
    Log << "Middle deviation from radius is " << mid_div << lendl;

    _max_div = max_div;
    _mid_div = mid_div;

    return r;
}

static double evaluate_radius(int meshtype, Object3D& obj, double& max_div, double& mid_div){
    if (meshtype == 0) return test_spheriable(obj, max_div, mid_div);
    if (meshtype == 1) return test_cilindrable(obj, max_div, mid_div);
    return 0;
}

static double evaluate_radius(int meshtype, Object3D& obj){
    double max_div, mid_div;
    evaluate_radius(meshtype, obj, max_div, mid_div);

    return 0;
}

static int generate_config_files(string dir){
    map<pair<int, ElasticModels>, double> brds;
    brds.insert({{0, EMOD_NGK}, 0.6});
    brds.insert({{1, EMOD_NGK}, 0.99});
    brds.insert({{0, EMOD_GNT}, 1.98});
    brds.insert({{1, EMOD_GNT}, 1.98});
    array<double, 1> mesh_steps = {/*0.1, 0.2,*/ 0.4/*, 0.8*/};
    array<int, 2> mesh_types = {0, 1};
    array<pair<ElasticModels, string>, 2> e_models = {pair<ElasticModels, string>{EMOD_NGK, "neogook"}, pair<ElasticModels, string>{EMOD_GNT, "gent"}};
    vector<tuple<double, int, string, double>> variants;
    variants.reserve(4 * 2 * 2 * 10);
    for(auto& i: e_models)
        for (auto& j: mesh_types){
            double brd = brds[{j, i.first}];
            for (int k = 1; k <= 10; ++k)
                for (auto l: mesh_steps)
                    variants.push_back(tuple<double, int, string, double>{l, j, i.second, (brd * k) / 10 + brd});
        }
    for (int i = 0; i < variants.size(); ++i){
        //if ((i + 1) % 4 != 1) continue;
        auto& v = variants[i];
        using json = nlohmann::json;
        json js;
        js["R"] = 10;
        js["l"] = 7;
        js["phi"] = 90;
        js["mesh_step"] = get<0>(v);
        js["mesh_type"] = get<1>(v);
        js["elastic_model"] = get<2>(v);
        js["thickness"] = 0.5;
        js["model_params"] = {1e6, 2.3};
        js["rel_allow_shift"] = 0.08;
        js["rel_max_shift"] = 0.2;
        double P = get<3>(v) * 1e6 / 3 * 0.5 / (10 * 133.3);
        js["P"] = P;
        double _P = ((int)(P * 100)) / 100;
        js["delta"] = 5e-7;
        js["velocity_decay"] = 1e-5;
        js["max_count_of_iterations"] = 40000000;
        js["frequency_of_printing"] = 250;
        js["save_directory"] = "../bin/";
        js["file_prefix"] = to_string(get<1>(v)) + "_" + get<2>(v) + "_" + to_string((int)_P*100) + "_" + to_string((int)(get<0>(v)*10));
        js["logfile_name"] = "logfile.txt";
        std::ofstream ff;
        ff.open(dir + "config_" + to_string(170 + i+1) + ".json");
        ff << js;
        ff.close();
    }
    return 0;
}

