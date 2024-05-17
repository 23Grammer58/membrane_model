#include "VirtualSuturer.h"
#include <boost/program_options.hpp>
#include "TestVirtualSuturer.h"

using namespace std;
using namespace World3d;

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

/*
 *      *        * |
 *      *        * | - linear part with height Hc
 *      *        * |
 *       *      * |
 *         *__*   | - elliptic part with height Hb
 */
std::function<std::array<double, 3>(double t)> get_line1(double R, double phi0, double Hc, double Hb){
    double ae = Hb, be = R*phi0 / 2;
    double l = 2*Hc + 2*(M_PI*ae*be + (ae-be)*(ae-be))/(ae + be);
    return [R, phi0, Hc, Hb, t1 = Hc / l, t2 = 1 - Hc / l](double t)->std::array<double, 3>{
        if (t < t1)
            return {R*cos(-phi0/2), R*sin(-phi0/2), Hc * (1 - t / t1)};
        else if (t < t2) {
            double teta = M_PI * (t - t1) / (t2 - t1);
            double phi = -phi0/2 * cos(teta);
            return {R*cos(phi), R*sin(phi), -Hb*sin(teta)};
        }
        else
            return {R*cos(+phi0/2), R*sin(+phi0/2), Hc * (t - t2) / t1};
    };
}

std::function<std::array<double, 3>(double t)> get_line2(double R, double phi0, double Hc, double Hb,
                                                         double a_bot = 0.8, double phi_bot = 3*M_PI_4){
    assert(a_bot < 0.95);
    assert(phi_bot >= M_PI_2 && phi_bot < M_PI);
    double ae = Hb, be = R*phi0 / 2;
    auto flat = [d = R*phi0, Hc, Hb, a_bot, phi_bot](double t)->std::array<double, 2>{
        if (t < 0.25){
            double l = 1 - t*4;
            return {-d/2, Hc*l};
        } else if (t < 0.5){
            double l = (t - 0.25)*4;
            double teta = phi_bot * l;
            return {-(1 - a_bot)*d/2 - a_bot*d/2*cos(teta), -Hb*sin(teta)};
        } else if (t < 0.75){
            double l = (t - 0.5)*4;
            double  x0 = -(1 - a_bot)*d/2-a_bot*d/2*cos(phi_bot) - d/2,
                    y0 = -Hb*sin(phi_bot),
                    k = Hb/(a_bot*d/2)*tan(phi_bot - M_PI_2);
            double x1 = 0, y1 = y0 - k*x0;
            if (!(y1 <= 0)) {
                double k_inv = (a_bot*d/2) / Hb * tan(M_PI - phi_bot);
                x1 = x0 - y0*k_inv, y1 = 0;
            }
            //create bezier curve
            std::array<double, 2> p[5] = {
                    {x0, y0},
                    {x0 + (x1 - x0)/3, y0 + (y1-y0)/3},
                    {x0 + 2*(x1 - x0)/3, -Hb},
                    {0, y0 + 2*(y1-y0)/3},
                    {0, 0}
            };
            for (auto & i : p) i[0] += d/2;
            std::array<double, 2> res = {0, 0};
            double r = 1 - l;
            double l2 = l*l, r2 = r*r;
            double l3 = l2*l, l4 = l2*l2;
            double r3 = r2*r, r4 = r2*r2;
            double a[5] = {r4, 4*r3*l, 6*r2*l2, 4*r*l3, l4};
            for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 2; ++j) res[j] += a[i]*p[i][j];
            return res;
        } else {
            double l = (t - 0.75)*4;
            return {d/2, Hc*l};
        }
    };
    return [R, flat](double t)->std::array<double, 3>{
        auto p = flat(t);
        return {R*cos(p[0]/R), R*sin(p[0]/R), p[1]};
    };
}

OzakiTempl MeshD::getOzaki() const {
    OzakiTempl templ;
    templ.s = 0;
    templ.m = 0;
    templ.w = 0;
    templ.D = D;
    templ.h = H;
    templ.beta = H_f;
    templ.alpha = (L_f - D)/2;
    return templ;
}
std::pair<bool, double> MeshD::check_Hf_condition() const {
    if (H_f < -H) return {false, -H - H_f};
    else return {true, 0.0};
}
std::pair<bool, double> MeshD::check_D_condition() const {
    double d2 = H*H + (L_f/2)*(L_f/2);
    if (d2 < (D/2 * D/2)) return {false, sqrt(d2) - D/2*(1+DBL_EPSILON)};
    else return {true, 0.0};
}
///return err code and approximate change to be applyeid to recover problem to unacceptable variable
std::pair<MeshD::ErrCode, double> MeshD::check_state() const {
    if (L_f <= 0 || D <= 0) return {ERR_UNRECOVERED, 0};
    if (H < 0) return {ERR_H, -H*(1+DBL_EPSILON)};
    auto r = check_Hf_condition();
    if (!r.first) return {ERR_H_f, r.second};
    r = check_D_condition();
    if (!r.first) return {ERR_D, r.second};
    return {OK, 0};
}

double MeshD::get_sutured_boundary_length() const {
    double d2 = H*H + (L_f/2)*(L_f/2);
    double Ll = sqrt(d2 - (D/2 * D/2));
    double phi = M_PI - (atan(L_f/ (2 * H)) + acos(D / (2*sqrt(d2))));
    return D*phi + 2*Ll;
}
double MeshD::get_free_boundary_length() const {
    return 2*sqrt(H_f*H_f + (L_f/2)*(L_f/2));
}
std::ostream& operator<<(std::ostream& out, const MeshD& t){
    return out << "{\n"
                  "\tmesh_size = " << t.mesh_size << "\n"
                  "\tMeshD{\n"
                  "\t\tD = " << t.D << "\n"
                  "\t\tL_f = " << t.L_f << "\n"
                  "\t\tH = " << t.H << "\n"
                  "\t\tH_f = " << t.H_f << "\n"
                  "\t}\n"
                  "}\n";
}
std::ostream& operator<<(std::ostream& out, const PhysProperties& t){
    return out << "{\n"
                  "\tHt = " << t.Ht << " mm\n"
                  "\tE = " << t.E << " kPa\n"
                  "\tP = " << t.P << " mm Hg\n"
                  "}\n";
}
void LeafMessures::print(std::ostream& out){
    out << "Lcoapt = " << length_coaptation << ", Hcentral = " << central_coaptation << " Heff = " << effective_height << "\n"
        << "isClosed = " << form_closed_valve <<  " closure_degree = " << closure_degree << "\n"
        << "billowing = " << billowing << " " << coaptation_area << "\n"
        << "l_free = " << free_edge_init_length << " l_free_init = " << free_edge_res_length << "\n";
}

TestVirtualSutureParams::TestVirtualSutureParams(){
    props.Ht = 0.5;
    props.E = 1000;
    props.P = 80;
    for (int i = 0; i < 3; ++i) {
        mhd[i].mesh_size = 1.0;
        mhd[i].D = 19; mhd[i].L_f = 3.35*2 + 19; mhd[i].H = 11.3; mhd[i].H_f = 2.1;
        phi_relations[i] = 1.0;
    }
} 
std::array<double, 6> TestVirtualSutureParams::getLeafCommisureAngles() const{
    std::array<double, 6> phi;
    DReal dphi = dist_suture / R;
    DReal fphi = 2*M_PI - 3*dphi;
    DReal c = phi_relations[0] + phi_relations[1] + phi_relations[2];
    auto phi_rel = phi_relations;
    for (int i = 0; i < 3; ++i) phi_rel[i] *= fphi / c;
    phi[0] = -phi_rel[0]/2; phi[1] = phi_rel[0]/2;
    phi[2] = phi[1] + dphi; phi[3] = phi[2] + phi_rel[1];
    phi[4] = phi[3] + dphi; phi[5] = phi[4] + phi_rel[2];
    return phi;
}
std::array<Point, 3> TestVirtualSutureParams::getCommissurePoints() const{
    std::array<double, 6> phi = getLeafCommisureAngles();
    std::array<Point, 3> comissure_points;
    for (int i = 0; i < 3; ++i){
        DReal phi0 = phi[2*((i+1)%3)+1];
        if (phi0 < 0) phi0 += 2*M_PI;
        DReal phi1 = phi[2*((i+2)%3)];
        if (phi1 < 0) phi1 += 2*M_PI;
        DReal lphi = (phi0 + phi1)/2;
        comissure_points[i] = Point(R*cos(lphi), R*sin(lphi), Hc);
    }
    return comissure_points;
}
void TestVirtualSutureParams::setQuasiOptWidth(){
    auto A = getCommissurePoints();
    std::array<double, 3> l;
    for (int i = 0; i < 3; ++i) l[i] = sqrt((A[(i+1)%3] - A[(i+2)%3]).squared_length());
    for (int i = 0; i < 3; ++i) mhd[i].L_f = l[i] * 2 / sqrt(3);// + dist_suture/2;
}
void TestVirtualSutureParams::parseInputArgs(int argc, char* argv[]){
    namespace po = boost::program_options;
    TestVirtualSutureParams::argc = argc, TestVirtualSutureParams::argv = argv;
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
                ("thickness,t", po::value<double>(&props.Ht)->default_value(0.5), "Leafs thickness, mm")
                ("young_modulus,E", po::value<double>(&props.E)->default_value(1000), "Young modulus, kPa")
                ("pressure,p", po::value<double>(&props.P)->default_value(80), "Diastolic pressure, Hg mm")
                ("mesh_step", po::value<double>()->default_value(1), "Step of generated mesh, mm")
                ("num_leafs,n", po::value<int>(&nleafs)->default_value(3), "Number of leafs")
                ("template_geometry0", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{19, 3.35*2 + 19, 11.3, 2.1}), "Parameters of template geometry {D, L_f, H, H_f} mm")
                ("template_geometry1", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{19, 3.35*2 + 19, 11.3, 2.1}), "Parameters of template geometry {D, L_f, H, H_f} mm")
                ("template_geometry2", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{19, 3.35*2 + 19, 11.3, 2.1}), "Parameters of template geometry {D, L_f, H, H_f} mm")
                ("root_directory", po::value<std::string>(&root_directory)->default_value(""), "Root directory for saving and reading of results")
                ("save_directory", po::value<std::string>(&result_directory)->default_value("default"), "Path relative root directory to save steps")
                ("view", "Create ImGui debug view window")
                ("save_steps", "Save sewing steps")
                ("save_collisions", "Save collision shapes")
                ("use_bezier", "Use bezier contacts instead bezier")
                ("use_aortic", "Use constraint from aorta")
                ("use_bending", "Turn on bending forces")
                ("use_clamped_bc", "Set clamped BC on sutured boundary")
                ("regenerate_expressions", "Regenerate expressions")
                ("use_opt_lwidth", "Use quasi-optimal leaflet width for the cilindric aorta")
                ("bnd.types", po::value<std::vector<int>>()->multitoken()->default_value(std::vector<int>{1, 1, 1}), "Types[3] of sew line (1 - simple, 2 - complex)")
                ("bnd.R", po::value<double>(&R)->default_value(12), "Radius of cilinder, mm")
                ("bnd.phi", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{1, 1, 1}), "Relations[3] of angle size of aortic suture lines")
                ("bnd.Hc", po::value<double>(&Hc)->default_value(1.5), "Vertical allowance of the line, mm")
                ("bnd.Hb", po::value<double>(&Hb)->default_value(9), "Deflection depth of the line, mm")
                ("dist_suture", po::value<double>(&dist_suture)->default_value(4), "Distance on cilinder aorta between sutured lines, mm")
                ("verbosity", po::value<int>(&verbose_level)->default_value(2), "Verbosity level")
                ("collision_thickness", po::value<double>(&collision_Ht_coef)->default_value(1.0), "Scaled final collision thickness")
                ("collision_shift", po::value<double>(&collision_Ht_shift)->default_value(0.0), "Scaled shift of collision objects")
                ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config);

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << cmdline_options << "\n";
            exit(0);
        }
        if (vm.count("version")){
            cout << "CilindricSuturer version 0.1.0\n";
            exit(0);
        }
        {
            auto& r = vm["bnd.types"].as<std::vector<int>>();
            auto& _phi = vm["bnd.phi"].as<std::vector<double>>();
            std::copy(r.data(), r.data() + 3, bnd_types.data());
            std::copy(_phi.data(), _phi.data() + 3, phi_relations.data());
        }
        if (vm.count("bnd.types"))
        with_gui_window = vm.count("view");
        save_steps = vm.count("save_steps");
        save_col_shapes = vm.count("save_collisions");
        with_bezier_contact = (vm.count("use_bezier") > 0);
        with_clamped_bc = (vm.count("use_clamped_bc") > 0);
        with_aorta = vm.count("use_aortic");
        with_bending_forces = vm.count("use_bending");
        use_quasi_optimal_leaf_width = (vm.count("use_opt_lwidth") > 0);
        regenerate_expressions = vm.count("regenerate_expressions");
        DReal mh = vm["mesh_step"].as<double>();
        for (int i = 0; i < 3; ++i) mhd[i].mesh_size = mh;
        for (int i = 0; i < 3; ++i){
            auto& tg = vm["template_geometry" + to_string(i)].as<std::vector<double>>();
            auto& ot = mhd[i];
#define IN_T(X, Y) if (tg.size() > X) ot.Y = tg[X]
        IN_T(0, D); IN_T(1, L_f); IN_T(2, H); IN_T(3, H_f);
#undef IN_T
        }
        if (nleafs <= 0 || nleafs > 3)
            throw std::runtime_error("Wrong number of leafs, expected value in range {1, 2, 3}");
        
        if (save_steps || save_col_shapes){
            path p(root_directory);
            if (!exists(status(p)) || !is_directory(status(p)))
                throw std::runtime_error("root_directory = \"" + root_directory + "\" is not exists");
            if (save_steps || save_col_shapes){
                p = root_directory + "/" + result_directory;
                if (!exists(status(p))) create_directory(p);
                else if (!is_directory(status(p)))
                    throw std::runtime_error("save_directory = \"" + p.string() + "\" but it's exists as nondirectory");
            }
        }
        if (!root_directory.empty()) result_directory = root_directory + "/" + result_directory;
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        exit(-1);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        exit(-2);
    }
}

std::map<int, LeafMessures> run_suture_simulation(TestVirtualSutureParams p){
    std::array<MeshD, 3> mhd = p.mhd;
    PhysProperties props = p.props;
    int nleafs = p.nleafs;
    std::array<DReal, 3> phi_relations = p.phi_relations;
    bool with_gui_window = p.with_gui_window;
    std::string root_directory = p.root_directory;
    std::string result_directory = p.result_directory;
    double R = p.R;
    double dist_suture = p.dist_suture;
    double Hc = p.Hc, Hb = p.Hb;
    std::array<int, 3> bnd_types = p.bnd_types;
    bool save_steps = p.save_steps;
    bool is_bezier = p.with_bezier_contact;
    bool with_aorta = p.with_aorta;
    bool with_bending_forces = p.with_bending_forces;
    bool without_clamped_bc = !p.with_clamped_bc;
    std::vector<int> ileafs = p.ileafs;
    bool save_col_shapes = p.save_col_shapes;
    std::array<DReal, 6> phi = p.getLeafCommisureAngles();
    int verbose_level = p.verbose_level;

    auto print_input = [&](){
        std::cout << "geom0:" << mhd[0] << "geom1:" << mhd[1] << "geom2:" << mhd[2]
                  << "PhysProperties:" << props << "\n"
                  << "with_gui_window = " << (with_gui_window ? "true" : "false") << "\n"
                  << "with_save_steps = " << (save_steps ? "true" : "false") << "\n"
                  << "root_directory = " << root_directory << "\n"
                  << "result_directory = " << result_directory << "\n"
                  << "{\n"
                     "\tnleafs = " << nleafs << "\n"   
                     "\tR = " << R << "\n"
                     "\tdist_suture = " << dist_suture << "\n"
                     "\tHc = " << Hc << "\n"
                     "\tHb = " << Hb << "\n"
                     "\tbnd.types = {" << bnd_types[0] << ", " << bnd_types[1] << ", " << bnd_types[2] << "}\n"
                     "\tphi_sizes = {" << phi_relations[0] << ", " << phi_relations[1] << ", " << phi_relations[2] << "}\n"
                     "\tphi_pos = {[" << phi[0] << ", " << phi[1] << "], [" << phi[2] << ", " << phi[3] << "], [" << phi[4] << ", " << phi[5] << "]}\n" 
                     "}\n";
    };
    if (verbose_level >= 0) print_input();
    std::array<std::function<Point(double t)>, 3> lines;
    for (int i = 0; i < nleafs; ++i){
        std::function<std::array<double, 3>(double t)> line;
        if (bnd_types[i] == 1)
            line = get_line1(R, phi[2*i+1] - phi[2*i], Hc, Hb);
        else if (bnd_types[i] == 2)
            line = get_line2(R, phi[2*i+1] - phi[2*i], Hc, Hb);
        else
            throw std::runtime_error("Faced unknown sew line type");
        DReal dphi = (phi[2*i+1] + phi[2*i])/2;

        lines[i] = [line, dphi](double t)->Point {
            auto x = line(t);
            auto phi_x = atan2(x[1], x[0]);
            auto R = hypot(x[0], x[1]);
            x[0] = R*cos(dphi + phi_x);
            x[1] = R*sin(dphi + phi_x);
            return Point(x[0], x[1], x[2]);
        };
        std::cout << "lines[" << i << "]: " << lines[i](0) << " " << lines[i](0.5) << " " << lines[i](1) << std::endl;
    }

    //specify suture parameters
    VirtualSuturerCil vs;
    vs  .setViewCtx(with_gui_window, p.argc, p.argv);
    if (with_aorta)
        vs.setCilindricAortic();
    if (save_steps)
        vs.setSaveDirectory(result_directory);
    {
        std::array<Point, 3> comissure_points = p.getCommissurePoints();
        std::array<DReal, 3> l;
        for (int i = 0; i < 3; ++i){
            l[i] = mhd[i].L_f;
            // DReal phi0 = phi[2*((i+1)%3)+1];
            // if (phi0 < 0) phi0 += 2*M_PI;
            // DReal phi1 = phi[2*((i+2)%3)];
            // if (phi1 < 0) phi1 += 2*M_PI;
            // DReal lphi = (phi0 + phi1)/2;
            // comissure_points[i] = Point(R*cos(lphi), R*sin(lphi), Hc);
        }
        vs.setCommissurePoints(comissure_points).setLeafLengths(l);
    }
    std::sort(ileafs.begin(), ileafs.end());
    std::array<AniMesh, 3> am;
    {
        int j = -1;
        for (int i = 0; i < nleafs && j < 0; ++i){
            if (mhd[i].D > 0 && mhd[i].L_f > 0 && mhd[i].H >= 0) j = i;
        }
        if (j < 0) throw std::runtime_error("Mesh geometry was not initialized");
        for (int i = 0; i < nleafs; ++i){
            if (mhd[i].D > 0 && mhd[i].L_f > 0 && mhd[i].H >= 0) am[i] = generate_ozaki(mhd[i].getOzaki(), mhd[i].mesh_size);
            else {
                double sc = mhd[i].L_f / mhd[j].L_f;
                auto it = std::lower_bound(ileafs.begin(), ileafs.end(), i);
                if (it != ileafs.end() && *it == i){
                    mhd[i].D = sc * mhd[j].D;
                    mhd[i].H = sc * mhd[j].H;
                    mhd[i].H_f = sc * mhd[j].H_f;
                    am[i] = generate_ozaki(mhd[i].getOzaki(), mhd[i].mesh_size);
                } else {
                    am[i] = am[j];
                    std::for_each(am[i].vertices.begin(), am[i].vertices.end(), [sc](auto& i){ i *= sc; });
                }
            }
        }
    }
    for (int i = 0; i < nleafs; ++i){
        auto mesh = convert_to_Mesh(am[i], "v:boundary_lbl");
        auto& lf = vs.pushLeaf(mesh, lines[i])
            .setCilindricClampField()
            .setPressure(operator""_mmHg(props.P)/1e3)
            .setYoungModulus(props.E)
            .setThickness(props.Ht)
            .setSutureBFieldTagName();
        if (is_bezier)
            lf.setContactSurfType(VirtualSuturerCil::LeafInfo::BEZIER_CILINDRIC);
        else 
            lf.setContactSurfType(VirtualSuturerCil::LeafInfo::PLANE);
    }
    if (without_clamped_bc){
        for (auto& f: vs.m_p.m_fsc) for (auto& b: f.m_bend) b.use_clamped_bc = false;
        for (auto& b: vs.m_p.m_csc.m_bend) b.use_clamped_bc = false;
        for (auto& b: vs.m_p.m_asc.m_bend) b.use_clamped_bc = false;
    }
    if (!with_bending_forces){
        vs.m_p.m_fsc.clear();
        std::for_each(vs.m_p.m_csc.m_bend.begin(), vs.m_p.m_csc.m_bend.end(), [](auto& i){ i.cE = 0.0; });
        std::for_each(vs.m_p.m_asc.m_bend.begin(), vs.m_p.m_asc.m_bend.end(), [](auto& i){ i.cE = 0.0; });
    }
    if (true) {
        double csigma = 3e-3;
        double cHt = 0.05;
        vs.m_p.m_csc.m_damps.csigma = csigma;
        vs.m_p.m_csc.m_contact.resize(0);
        vs.m_p.m_csc.m_solve_ctx.maxits = 30;
        auto push_contact = [&c = vs.m_p.m_csc.m_contact](double cP, double cH){ 
            VirtualSuturerCil::SDFForceCtx ctx; ctx.cP = cP, ctx.cHt = cH;
            ctx.cshift = 1;
            c.push_back({ctx, ctx, ctx}); 
        };
        push_contact(0.05, cHt);
        push_contact(0.1, cHt);
        push_contact(0.2, cHt);
        push_contact(0.3, cHt);
        push_contact(0.4, cHt);
        push_contact(0.5, cHt);
        push_contact(0.8, cHt);
        push_contact(1.0, cHt);
        for (int i = 0; i < 5; ++i) push_contact(2.0, cHt);
        push_contact(4.0, cHt);
    }
    {
        VirtualSuturerCil::CustomCxt cxt;
        cxt.m_elast.first = vs.m_p.m_csc.m_elast, cxt.m_elast.second = true;
        cxt.m_bend.first = vs.m_p.m_csc.m_bend, cxt.m_bend.second = true;
        cxt.m_pres.first = vs.m_p.m_csc.m_pres, cxt.m_pres.second = true;
        cxt.m_free_edge.second = false;
        cxt.m_aortic.second = false;
        // for (auto& pctx: cxt.m_pres.first) pctx.cP = 2;
        cxt.m_contact.first = vs.m_p.m_csc.m_contact.back(); cxt.m_contact.second = true;
        for (int i = 0; i < 3; ++i) {
            cxt.m_contact.first[i].cP = 1.2;
            cxt.m_contact.first[i].cHt = 0.5;
            cxt.m_contact.first[i].cshift = 0.9;
        }
        cxt.m_damps.second = true; cxt.m_damps.first.csigma = 5e-2;
        cxt.m_solve_ctx.maxits = 30;
        for (int i = 0; i < 15; ++i)
            vs.m_p.m_custom.push_back(cxt);
    }
    {
        std::vector<VirtualSuturerCil::SolverCtx*> sctx; 
        sctx.push_back(&vs.m_p.m_fmc.m_solve_ctx);
        for (int i = 0; i < vs.m_p.m_fsc.size(); ++i)
            sctx.push_back(&vs.m_p.m_fsc[i].m_solve_ctx);
        sctx.push_back(&vs.m_p.m_csc.m_solve_ctx); 
        sctx.push_back(&vs.m_p.m_asc.m_solve_ctx); 
        for (int i = 0; i < vs.m_p.m_custom.size(); ++i)
            sctx.push_back(&vs.m_p.m_custom[i].m_solve_ctx); 
        for (auto ctx: sctx) ctx->verbose_level = verbose_level;     
    }
    {
        std::vector<VirtualSuturerCil::ContactForce*> cctx;
        cctx.push_back(&vs.m_p.m_asc.m_contact);
        for (auto& c: vs.m_p.m_custom) if (c.m_contact.second) cctx.push_back(&c.m_contact.first);
        for (auto& c: cctx) for (int i = 0; i < 3; ++i) {
            (*c)[i].cHt = p.collision_Ht_coef;
            (*c)[i].cshift = p.collision_Ht_shift;
            (*c)[i].cP = 2.0;
        }
    }
    
    if (p.regenerate_expressions){
        vs.m_regenerate_elastic_force = true;
        vs.m_regenerate_bending_force = true;
    }
    //perform calculations
    vs.setup();
    for (int i = 0, j = 0; i < 3; ++i){
        while (j < ileafs.size() && ileafs[j] < i) ++j;
        if (j < ileafs.size() && ileafs[j] == i)
            vs.suture_leaf(i);
        else {
            auto bnd = vs.geometrical_shifting_leaf(i);
            bnd.addShift(1);
            ///we move nonsutured leaflets to avoid it's action on coaptation evaluating
            Vector n = CGAL::cross_product(vs.m_commissure.m_n, vs.m_commissure.m_A[(i+2)%3] - vs.m_commissure.m_A[(i+1)%3]);
            n /= sqrt(n.squared_length());
            for (auto v: vs.m_objs[i].m_mesh.vertices())
                vs.m_objs[i].m_x[v] += -2*vs.m_commissure.m_R*n;
        }    
    }
    //evaluate surgical parameters
    Object3D dup;
    Messurer mes(dup, vs.m_objs[0], vs.m_objs[1], vs.m_objs[2], vs.m_commissure.m_n);
    mes.setMargin(0.5*props.Ht);
    Point C; //< triple contact point
    {
        auto pd = CilindricContactSurface::computeDilationPointFromLengths(
                                {vs.m_commissure.m_A[0], vs.m_commissure.m_A[1], vs.m_commissure.m_A[2]},
                                {vs.m_linfo[0].m_l, vs.m_linfo[1].m_l, vs.m_linfo[2].m_l});
        C = CGAL::ORIGIN + pd.C;
    }
    auto specialCase = [&vs, &ileafs, C, H = 0.5*props.Ht](const Messurer::ObjectShape& from, const Messurer::ObjectShape& to, std::vector<Messurer::CollisionInfo::CPD>& info)->bool{
        int from_i = -1, to_i = -1;
        for (int i = 0; i < 3; ++i){
            if (&vs.m_objs[i] == from.obj) from_i = i;
            if (&vs.m_objs[i] == to.obj) to_i = i;
        }
        auto it = std::lower_bound(ileafs.begin(), ileafs.end(), from_i);
        if (it == ileafs.end() || *it != from_i) return false;
        
        Point O = vs.m_commissure.m_A[3 - (from_i + to_i)];
        Vector n = CGAL::cross_product(vs.m_commissure.m_n, O - C);
        n /= sqrt(n.squared_length());
        n *= (n * (vs.m_commissure.m_A[from_i] - C) > 0) ? 1 : -1;
        
        bool changed = false;
        for (int v = 0; v < from.mesh.points.size(); ++v){
            Point p = from.mesh.points[v];
            double d = (p - C)*n;
            if (d > -H) {
                Messurer::CollisionInfo::CPD c;
                c.v = V_ind(v);
                c.f = F_ind();
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
            Point O = vs.m_commissure.m_A[3 - (i + j)];
            Vector n = CGAL::cross_product(vs.m_commissure.m_n, O - C);
            n /= sqrt(n.squared_length());
            n *= (n * (vs.m_commissure.m_A[i] - C) > 0) ? 1 : -1;
            mes.m_planes(i, j) = {Messurer::Plane_3(C,  n), true};
            mes.m_planes(j, i) = {Messurer::Plane_3(C, -n), true};
        }    
    mes.computeHalfCoaptScans();
    std::map<int, LeafMessures> lmes;
    auto distr = mes.computeHalfCoaptDistrib(500);
    mes.uniteCollidingCuspHalves(); 
    auto fdistr = mes.computeCoaptDistribution({1000, 1000, 1000});
    auto Hcentr = mes.computeHc();
    auto H = distr.getCoaptH();
    auto coaptStat = mes.computeCoaptStatus();
    bool isClosed = ((coaptStat & (1 << 6)) > 0);
    auto bill = mes.computeBillowing();
    auto areas = mes.computeColArea(); 

    auto make_colfilter = [&mes](int ileaf){
        return [&mes, &col = mes.m_colMap[ileaf]](F_ind f)->bool{
            auto lbl = col.face_lbl[col.remap[f]];
            bool res = true;
            for (int k = 0; k < 3; ++k)
                res &= ((lbl & (7 << 3 * k)) > 0);
            return res;
        };  
    }; 
    auto make_nocolfilter = [&make_colfilter](int i){ return [filter = make_colfilter(i)](F_ind f){ return !filter(f); }; };

    if (ileafs.size() == 1){
        Object3D leafs[3];
        leafs[2] = Object3D(convert_to_Mesh(am[ileafs[0]], "v:boundary_lbl"));
        for (auto v: leafs[2].m_mesh.vertices()) leafs[2].m_x[v] = vs.m_objs[ileafs[0]].m_x[v] + Vector(0.5, 0, 0);
        
        for (int i = 0; i < 2; ++i){
            leafs[i] = leafs[2];
            for (auto v: leafs[i].m_mesh.vertices()) {
                double x = leafs[2].m_x[v][0], y = leafs[2].m_x[v][1], z = leafs[2].m_x[v][2];
                double phi = atan2(y, x);
                double R = std::hypot(x, y);
                leafs[i].m_x[v] = Point(R*cos(2*M_PI/3*(i+1) + phi), R*sin(2*M_PI/3*(i+1) + phi), z);
            }
        }
        World w;
        w.addObject3D(std::move(leafs[2]), 1);
        w.addObject3D(std::move(leafs[0]), 1);
        w.addObject3D(std::move(leafs[1]), 1);
        unique_ptr<CollisionManagerBase> colMan = std::make_unique<BridsonCollisionManager>();
        w.setCollider(std::move(colMan));
        BridsonCollisionManager* pcol = reinterpret_cast<BridsonCollisionManager*>(w.getCollider());
        DReal P = p.props.P;
        pcol->set_Kspring(P * 250);
        pcol->set_Dt(0.1);
        pcol->setCollisionStatus(BridsonCollisionManager::SelfCollision | BridsonCollisionManager::DynamicVsDynamicCollission);
        pcol->setAsInitialNonCollidingState();
        pcol->set_MaxRelaxIts(0);
        pcol->set_RelaxEps(0.5);
        for (auto& o: w.objs()){
            auto id = o.first;
            DReal Ht = 0.05;

            CollisionThicknessCompressed<CollisionThicknessConst> cthickness{CollisionThicknessConst(Ht)};
            CollisionThickness Htf(cthickness);
            pcol->set_thickness(id, Htf);
        }
        for (auto v: leafs[2].m_mesh.vertices()) leafs[2].m_x[v] = leafs[2].m_next_x[v] = leafs[2].m_x[v] - Vector(1.1, 0, 0);
        for (int i = 0; i < 2; ++i){
            for (auto v: leafs[i].m_mesh.vertices()) {
                double x = leafs[2].m_next_x[v][0], y = leafs[2].m_next_x[v][1], z = leafs[2].m_next_x[v][2];
                double phi = atan2(y, x);
                double R = std::hypot(x, y);
                leafs[i].m_x[v] = leafs[i].m_next_x[v] = Point(R*cos(2*M_PI/3*(i+1) + phi), R*sin(2*M_PI/3*(i+1) + phi), z); 
            }
        }
        
        pcol->findCollisions();
        pcol->solveCollisions();
        for (int i = 0; i < 3; ++i) leafs[i].save(result_directory +"/"+ "tempres" + std::to_string(i) + ".vtk");
    }

    for (int ileaf: ileafs){
        if (ileaf < 0 || ileaf >=3 ) continue;
        static const int FIXED = 2, FREE = 1|4|8;
        auto& leaf = vs.m_objs[ileaf];
        double l_free = 0, l_free_init = 0;
        V_ind FE_centr;
        {
            for (auto e: leaf.m_mesh.edges()) {
                auto f = face_around(leaf.m_mesh, e);
                if (f[0].second && f[1].second)
                    continue;
                auto v = vert_around(leaf.m_mesh, e);
                auto vop = vert_opposite(leaf.m_mesh, e);
                if ((leaf.m_boundary[v[0]] & FREE) && (leaf.m_boundary[v[1]] & FREE)) {
                    double loc = sqrt((leaf.m_x[v[0]] - leaf.m_x[v[1]]).squared_length());
                    double loc_init = sqrt((leaf.m_x0[v[0]] - leaf.m_x0[v[1]]).squared_length());
                    l_free += loc, l_free_init += loc_init;
                }
            }
            for (auto v: leaf.m_mesh.vertices()) if (leaf.m_boundary[v] == 5) FE_centr = v;
        }
        LeafMessures lm;
        lm.half_distrib = distr.data[ileaf];
        lm.union_distrib = fdistr[ileaf];
        lm.central_coaptation = Hcentr;
        lm.length_coaptation = std::max(distr.getCoaptH(ileaf, (ileaf+1)%3), distr.getCoaptH(ileaf, (ileaf+2)%3));
        lm.billowing = bill.getBillowing(ileaf);
        auto bil_plane = CGAL::plane_from_points<CGAL::Simple_cartesian<DReal>>(bill.plane[0], bill.plane[1], bill.plane[2]);
        lm.effective_height = sqrt((leaf.m_x[FE_centr] - bil_plane.projection(leaf.m_x[FE_centr])).squared_length());
        lm.coaptation_area = areas[ileaf];
        lm.form_closed_valve = isClosed;
        lm.free_edge_init_length = l_free_init, lm.free_edge_res_length = l_free;
        {
            double d = std::numeric_limits<double>::max();
            for (auto v: leaf.m_mesh.vertices()){
                Vector cx = C - leaf.m_x[v];
                cx = cx - (cx * vs.m_commissure.m_n)*vs.m_commissure.m_n;
                double ld = sqrt(cx.squared_length());
                if (ld < d) d = ld;
            }
            lm.closure_degree = d;
        }
        lmes.insert({ileaf, std::move(lm)});
        if (save_col_shapes && !result_directory.empty()){
            lmes[ileaf].print();
            leaf.save(result_directory +"/"+ mes.getCusps()[ileaf]->name + "res.vtk");
            Messurer::saveObjectShape(mes.m_colission_shapes[ileaf], result_directory +"/"+ mes.getCusps()[ileaf]->name + "_mes_init.stl", "v:point");
            Messurer::saveObjectShape(mes.m_colission_shapes[ileaf], result_directory +"/"+ mes.getCusps()[ileaf]->name + "_mes_init_col.stl",
                                      "v:point",make_colfilter(ileaf));
            Messurer::saveObjectShape(mes.m_colission_shapes[ileaf], result_directory +"/"+ mes.getCusps()[ileaf]->name + "_mes_init_nocol.stl",
                                      "v:point",make_nocolfilter(ileaf));
            Messurer::saveObjectShape(mes.m_colission_shapes[ileaf], result_directory +"/"+ mes.getCusps()[ileaf]->name + "_mes_real_nocol.stl",
                                      "v:x", make_nocolfilter(ileaf));
            Messurer::saveObjectShape(mes.m_colission_shapes[ileaf], result_directory +"/"+ mes.getCusps()[ileaf]->name + "_mes_real_col.stl",
                                      "v:x", make_colfilter(ileaf));
            Object3D(mes.m_hcScan(ileaf, (ileaf+1)%3).to_Mesh()).save(
                        result_directory  + "/half_scan_" + std::to_string(ileaf) + "_" + std::to_string((ileaf+1)%3) + ".stl");   
            Object3D(mes.m_hcScan(ileaf, (ileaf+2)%3).to_Mesh()).save(
                        result_directory  + "/half_scan_" + std::to_string(ileaf) + "_" + std::to_string((ileaf+2)%3) + ".stl"); 
            Object3D(mes.m_coaptScan[ileaf].to_Mesh()).save(result_directory + "/unite_scan_" + std::to_string(ileaf) + ".stl");
            Messurer::saveCoaptDistribCSV(fdistr, result_directory + "/full_distrib.csv");
        }

    }

    return lmes;
}