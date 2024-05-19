#include "examples.h"

static void save_in_all_formats(MetaTriMesh& mtm, const std::string name, bool verbose = true){
    VTK_File().write(name + "_txt.vtk", mtm);
    VTK_File().write(name + "_bin.vtk", mtm, true);
    TMH_File().write_txt(name + "_txt.tmh", mtm);
    TMH_File().write_bin(name + "_bin.tmh", mtm);
    if (verbose)
        std::cout << "Save mesh in\n" 
            << "\t" << (name + "_txt.vtk") << "\n"
            << "\t" << (name + "_bin.vtk") << "\n"
            << "\t" << (name + "_txt.tmh") << "\n"
            << "\t" << (name + "_bin.tmh") << "\n";
}

static std::pair<bool, double> read_double(const std::string& s){
    std::istringstream iss(s);
    double f = NAN;
    iss >> f; 
    return {iss.eof() && !iss.fail(), f};    
}

static bool read_bool(const std::string& s){
    std::string s0 = s;
    std::transform(s0.begin(), s0.end(), s0.begin(), [](char c){ return std::tolower(c); });
    if (s0 == "true") return true;
    else if (s0 == "false") return false;
    int r = stoi(s0);
    if (r == 0) return false;
    if (r == 1) return true;
    throw std::runtime_error("Can't convert \"" + s + "\" to boolean");
    return false;
}

std::ostream& stretch_rectangle::print_param_list(std::ostream& out) const{
    out << "a [mm], b [mm], ht [mm], mesh_h [mm], dx [mm], dy [mm], theta_f [rad], theta_s [rad] is_2d [bool] is_fixed [bool]";
    return out;
}
std::ostream& stretch_rectangle::print_state(std::ostream& out) const {
    out << "a = " << a << " [mm], b = " << b << " [mm], ht = " << ht << " [mm], mesh_h = " 
        << mesh_h << " [mm], dx = " << dx <<" [mm], dy = " << dy << " [mm], theta_f = " 
        << theta_f << " [rad], theta_s = " << theta_s << " [rad]" 
        << ", is_2d = " << (is_2d ? "true" : "false") << ", is_fixed = " << (is_fixed ? "true" : "false");
    return out;
}
stretch_rectangle& stretch_rectangle::init(int argc, char* argv[]){
    double* p[]{&a, &b, &ht, &mesh_h, &dx, &dy, &theta_f, &theta_s};
    for (int i = 0; i < argc && i < sizeof(p) / sizeof(double*); ++i){
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
            std::cout << "Usage ./prog stretch_rectangle "; print_param_list() << std::endl;
            exit(0);
        }
        auto r = read_double(argv[i]);
        if (!r.first) break;
        *(p[i]) = r.second;
    }
    bool* bp[]{&is_2d, &is_fixed};
    for (int i = sizeof(p) / sizeof(double*), j = 0; i < argc && j < sizeof(bp) / sizeof(bool*); ++i, ++j){
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
            std::cout << "Usage ./prog stretch_cross_piece "; print_param_list() << std::endl;
            exit(0);
        }
        *(bp[j]) = read_bool(argv[i]);
    }
    return *this;
}
int stretch_rectangle::generate() const {
    AniMesh am = generate_rectangle(a, b, mesh_h);
    Mesh m = convert_to_Mesh(am, "v:bnd");
    auto lbl = m.property_map<V_ind, int>("v:bnd").first;
    auto x = m.add_property_map<V_ind, Point>("v:x").first;
    auto cl = m.add_property_map<E_ind, Vector>("be:clamped_direction", Vector(0, 0, 0)).first;
    auto htt = m.add_property_map<F_ind, double>("f:thickness", ht);
    for (auto e: m.edges()){
        auto v = vert_around(m, e);
        if ((lbl[v[0]] & 4) && (lbl[v[1]] & 4)){
            cl[e] = Vector(-1, 0, 0);
        } else if ((lbl[v[0]] & 2) && (lbl[v[1]] & 2)){
            cl[e] = Vector(+1, 0, 0);
        } else if ((lbl[v[0]] & 1) && (lbl[v[1]] & 1)){
            cl[e] = Vector(0, (m.point(v[0])[1] > 0 ? 1 : -1), 0);
        }
    }
    for (auto v: m.vertices()) {
        x[v] = m.point(v);
        if (lbl[v] != 0){ 
            if (lbl[v] & 2)
                x[v] += Vector(dx, 0, 0);
            if ((lbl[v] & 1) && m.point(v)[1] > 0)
                x[v] += Vector(0, dy, 0);  
            if (is_fixed)  
                lbl[v] = 0xF;
            else {        
                if ((lbl[v] & 1) && (lbl[v] & 4) && m.point(v)[1] < 0)
                    lbl[v] = 0xF;
                else if (lbl[v] & (4|2))
                    lbl[v] = 0xD;
                else 
                    lbl[v] = 0xE;
            }
        } else 
            lbl[v] = is_2d ? 0x4 : 0;
    }
    Vector f{cos(theta_f), sin(theta_f), 0}, s{cos(theta_s), sin(theta_s), 0};
    auto fib_f = m.add_property_map<F_ind, Vector>("f:fiber_f", f);
    auto fib_s = m.add_property_map<F_ind, Vector>("f:fiber_s", s);
    std::string name = "stretch_rectangle";
    MetaTriMesh mtm(m);
    mtm.readTagDataFromMesh<V_ind>(m, "v:bnd", "v:bnd");
    mtm.readTagDataFromMesh<V_ind>(m, "v:x", "v:x");
    mtm.readTagDataFromMesh<E_ind>(m, "be:clamped_direction", "be:clamped_direction");
    mtm.readTagDataFromMesh<F_ind>(m, "f:thickness", "f:thickness");
    mtm.readTagDataFromMesh<F_ind>(m, "f:fiber_f", "f:fiber_f");
    mtm.readTagDataFromMesh<F_ind>(m, "f:fiber_s", "f:fiber_s");
    save_in_all_formats(mtm, name);

    return 0;
}

std::ostream& stretch_cross_piece::print_param_list(std::ostream& out) const{
    out << "a [mm], ac [mm], r [mm], ht [mm], mesh_h [mm], dx0 [mm], dy0 [mm], dx1 [mm], dy1 [mm], theta_f [rad], theta_s [rad], fixed [bool], r1 [mm], r2 [mm], r3 [mm]";
    return out;
}

std::ostream& stretch_cross_piece::print_state(std::ostream& out) const{
    out << "a = " << a << " [mm], ac = " << ac << " [mm], r = " << r << " [mm], ht = " << ht << " [mm], mesh_h = " << mesh_h << " [mm], "
        << "dx0 = " << dx0 << " [mm], dy0 = " << dy0 << " [mm], dx1 = " << dx1 << " [mm], dy1 = " << dy1 << " [mm], "
        << "theta_f = " << theta_f << " [rad], theta_s = " << theta_s << " [rad], fixed = " << is_fixed << " [bool]"
        << ", r1 = " << r1 << " [mm], r2 = " << r2 << " [mm], r3 = " << r3 << " [mm]";
    return out;
}
stretch_cross_piece& stretch_cross_piece::init(int argc, char* argv[]){
    double* p[]{&a, &ac, &r, &ht, &mesh_h, &dx0, &dy0, &dx1, &dy1, &theta_f, &theta_s};
    for (int i = 0; i < argc && i < sizeof(p) / sizeof(double*); ++i){
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
            std::cout << "Usage ./prog stretch_cross_piece "; print_param_list() << std::endl;
            exit(0);
        }
        auto r = read_double(argv[i]);
        if (!r.first) break;
        *(p[i]) = r.second;
    }
    for (int i = sizeof(p) / sizeof(double*), j = 0; i < argc && j < 1; ++i, ++j){
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
            std::cout << "Usage ./prog stretch_cross_piece "; print_param_list() << std::endl;
            exit(0);
        }
        is_fixed = read_bool(argv[i]);
    }
    double rdef = r;
    double* p1[]{&r1, &r2, &r3};
    for (int i =  sizeof(p) / sizeof(double*) + 1, j=0; i < argc && j < sizeof(p1) / sizeof(double*); ++i, ++j){
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
            std::cout << "Usage ./prog stretch_cross_piece "; print_param_list() << std::endl;
            exit(0);
        }
        auto r = read_double(argv[i]);
        if (!r.first) break;
        *(p1[j]) = rdef = r.second;
    }
    for (int i = 2; i >= 0 && *(p1[i]) < 0; --i) 
        *(p1[i]) = rdef;
    return *this;
}

int stretch_cross_piece::generate() const{
    AniMesh am = generate_cross_piece({a, a}, {ac, ac}, {r, r1, r2, r3}, mesh_h);
    Mesh m = convert_to_Mesh(am, "v:bnd", "f:region");
    auto htt = m.add_property_map<F_ind, double>("f:thickness", ht);
    auto x = m.add_property_map<V_ind, Point>("v:x").first;
    auto lbl = m.property_map<V_ind, int>("v:bnd").first;
    auto cl = m.add_property_map<E_ind, Vector>("be:clamped_direction", Vector(0, 0, 0)).first;
    for (auto e: m.edges()){
        auto v = vert_around(m, e);
        if ((lbl[v[0]] & 1) && (lbl[v[1]] & 1))
            cl[e] = Vector(+1, 0, 0);
        else if ((lbl[v[0]] & 2) && (lbl[v[1]] & 2)) 
            cl[e] = Vector(0, +1, 0); 
        if ((lbl[v[0]] & 4) && (lbl[v[1]] & 4))
            cl[e] = Vector(-1, 0, 0);
        else if ((lbl[v[0]] & 8) && (lbl[v[1]] & 8)) 
            cl[e] = Vector(0, -1, 0);      
    }
    for (auto v: m.vertices()) {
        x[v] = m.point(v);
        if (lbl[v] != 0){ 
            if (lbl[v] & 1) {
                x[v] += Vector(dx0, 0, 0);
                lbl[v] = is_fixed ? 0xF : (0x8 | 0x4 | 0x1);
            } else if (lbl[v] & 2) {
                x[v] += Vector(0, dy0, 0); 
                lbl[v] = is_fixed ? 0xF : (0x8 | 0x4 | 0x2);  
            } else if (lbl[v] & 4) {
                x[v] -= Vector(dx1, 0, 0);
                lbl[v] = is_fixed ? 0xF : (0x8 | 0x4 | 0x1);
            } else if (lbl[v] & 8) {
                x[v] -= Vector(0, dy1, 0);
                lbl[v] = is_fixed ? 0xF : (0x8 | 0x4 | 0x2);
            } else 
                lbl[v] = 0x4;
        } else 
            lbl[v] = 0x4;
    }
    Vector f{cos(theta_f), sin(theta_f), 0}, s{cos(theta_s), sin(theta_s), 0};
    auto fib_f = m.add_property_map<F_ind, Vector>("f:fiber_f", f);
    auto fib_s = m.add_property_map<F_ind, Vector>("f:fiber_s", s);
    std::string name = "stretch_cross";
    MetaTriMesh mtm(m);
    mtm.readTagDataFromMesh<F_ind>(m, "f:thickness", "f:thickness");
    mtm.readTagDataFromMesh<V_ind>(m, "v:x", "v:x");
    mtm.readTagDataFromMesh<V_ind>(m, "v:bnd", "v:bnd");
    mtm.readTagDataFromMesh<E_ind>(m, "be:clamped_direction", "be:clamped_direction");
    mtm.readTagDataFromMesh<F_ind>(m, "f:fiber_f", "f:fiber_f");
    mtm.readTagDataFromMesh<F_ind>(m, "f:fiber_s", "f:fiber_s");
    mtm.readTagDataFromMesh<F_ind>(m, "f:region", "f:region");

    save_in_all_formats(mtm, name);

    return 0;
}

std::ostream& inflation_circle::print_param_list(std::ostream& out) const {
    out << "R [mm], ht [mm], mesh_h [mm], theta_f [rad], theta_s [rad], R_c [mm], K [mm^-1]";
    return out;
}
std::ostream& inflation_circle::print_state(std::ostream& out) const {
    out << "R = " << R << " [mm], ht = " << ht << " [mm], mesh_h = " << mesh_h << " [mm], theta_f = " << theta_f << " [rad], theta_s = " << theta_s << " [rad]";
    if (R_c > 0 && R_c < R) out << ", R_c = " << R_c << " [mm]";
    if (K != 0) out << ", K = " << K << " [mm^-1]";
    return out;
}
inflation_circle& inflation_circle::init(int argc, char* argv[]){
    double* p[]{&R, &ht, &mesh_h, &theta_f, &theta_s, &R_c, &K};
    for (int i = 0; i < argc && i < sizeof(p) / sizeof(double*); ++i){
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){
            std::cout << "Usage ./prog inflation_circle "; print_param_list() << std::endl;
            exit(0);
        }
        auto r = read_double(argv[i]);
        if (!r.first) break;
        *(p[i]) = r.second;
    }
    return *this;
}
int inflation_circle::generate() const {
    AniMesh am;
    bool with_coaxial = !(R_c <= 0 || R_c >= R); //(K <= 0.0 && (R_c <= 0 || R_c >= R))
    std::string region_lbl = "";
    if (!with_coaxial)
        am = generate_circle(R, mesh_h);
    else {
        am = generate_circle_with_coaxial(R, R_c / R, mesh_h, K * R);   
        region_lbl = "f:region";
    } 
    Mesh m = convert_to_Mesh(am, "v:bnd", region_lbl);
    if (with_coaxial && abs(1. - K * R) < 1e-9)
        for (auto v: m.vertices()){
            Point p = m.point(v);
            std::array<double, 3> pr{p[0], p[1], p[2]};
            double len = std::hypot(pr[0], pr[1], pr[2]);
            m.point(v) = Point(pr[0]/len, pr[1]/len, pr[2]/len);
        }
    auto lbl = m.property_map<V_ind, int>("v:bnd").first;
    auto htt = m.add_property_map<F_ind, double>("f:thickness", ht);
    auto cl = m.add_property_map<E_ind, Vector>("be:clamped_direction", Vector(0, 0, 0)).first;
    for (auto e: m.edges()){
        auto v = vert_around(m, e);
        if ((lbl[v[0]] & 1) && (lbl[v[1]] & 1)){
            Vector l = ((m.point(v[0]) - CGAL::ORIGIN) + (m.point(v[1]) - CGAL::ORIGIN)) / 2;
            double phi = atan2(l[1], l[0]);
            cl[e] = Vector(cos(phi), sin(phi), 0);
        }
    }
    for (auto v: m.vertices()) {
        if (lbl[v] & 1)
            lbl[v] = 0xF;
        else    
            lbl[v] = 0; 
    }
    Vector f{cos(theta_f), sin(theta_f), 0}, s{cos(theta_s), sin(theta_s), 0};
    auto fib_f = m.add_property_map<F_ind, Vector>("f:fiber_f", f);
    auto fib_s = m.add_property_map<F_ind, Vector>("f:fiber_s", s);
    std::string name = "inflation_circle";
    MetaTriMesh mtm(m);
    mtm.readTagDataFromMesh<V_ind>(m, "v:bnd", "v:bnd");
    mtm.readTagDataFromMesh<E_ind>(m, "be:clamped_direction", "be:clamped_direction");
    mtm.readTagDataFromMesh<F_ind>(m, "f:thickness", "f:thickness");
    mtm.readTagDataFromMesh<F_ind>(m, "f:fiber_f", "f:fiber_f");
    mtm.readTagDataFromMesh<F_ind>(m, "f:fiber_s", "f:fiber_s");
    if (with_coaxial)
        mtm.readTagDataFromMesh<F_ind>(m, "f:region", "f:region");
    save_in_all_formats(mtm, name);

    return 0;
}