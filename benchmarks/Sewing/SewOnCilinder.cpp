//
// Created by Liogky Alexey on 20.05.2022.
//

#include <cmath>
#include <functional>
#include <array>
#include "AVSim/LeafTemplates/TemplatesCollection.h"
#include "AVSim/Core/TriangularMeshHelpers.h"
#include "AVSim/AorticSimulator/AVSimulator.h"
#include "AVSim/Core/MeshGen/StretchFigure.h"
#include "AVSim/Core/Forces/Forces.h"
#include "AVSim/Core/ForceAppliers/ForceAppliers.h"
#include "AVSim/Core/NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverKinsol.h"
#include <boost/program_options.hpp>
#include <AVSim/Core/Collision/BridsonCollisionManager.h>
#include "AVSim/AorticSimulator/VirtSuture/ContactSurface.h"
#include <AVSim/Core/Forces/ContactForce.h>
#include <AVSim/Core/Forces/SDFForce.h>
#include <regex>
#include <AVSim/Core/Collision/CollisionThickness.h>

using namespace std;
using namespace World3d;

#if __cplusplus >= 201703L
using namespace std::filesystem;
#else
using namespace std::experimental::filesystem;
#endif

struct CilDist: public SignedDistanceField{
    Vector C;
    Vector n;
    DReal R;
    SDF operator()(const Vector& x) const override;
    CilDist(Vector _C, Vector _n, DReal _R): C{_C}, n{_n}, R{_R} {}
};
SignedDistanceField::SDF CilDist::operator()(const Vector& x) const{
    SDF d;
    Vector rv = (x - C) - ((x - C)*n)*n;
    DReal r = sqrt(rv.squared_length());
    d.sdf = R - r;
    d.grad_sdf = -rv / r;
    Vector f = CGAL::cross_product(n, d.grad_sdf);
    for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j)
        d.grad_grad_sdf(i, j) = f[i]*f[j] / r;
    return d;
}

static bool ObSaveTag(const Mesh& m, string filename, string tag){
    ofstream ob(filename, ios::binary | ios::trunc);
    if (!ob) return false;
    auto m_x = m.property_map<V_ind, Point>(tag).first;
    for (auto& v: m.vertices()){
        std::array<double, 3> P = {m_x[v][0], m_x[v][1], m_x[v][2]};
        ob.write(reinterpret_cast<const char *> ((P.data())), 3 * sizeof(double));
    }
    ob.close();
    return true;
}

static bool ObSaveTag(Object3D& obj, string filename, string tag){
    return ObSaveTag(obj.m_mesh, filename, tag);
}

static bool ObReadTag(Mesh& m, string filename, string tag){
    ifstream ob(filename, std::ios::binary);
    if (!ob) return false;
    auto m_x = m.add_property_map<V_ind, Point>(tag);
    for (auto& v: m.vertices()){
        std::array<double, 3> P;
        ob.read(reinterpret_cast<char *> ((P.data())), 3 * sizeof(double));
        m_x.first[v] = Point(P[0], P[1], P[2]);
    }
    ob.close();
    return true;
}

static bool ObReadTag(Object3D& obj, string filename, string tag){
    return ObReadTag(obj.m_mesh, filename, tag);
}

/*
 *      *        * |
 *      *        * | - linear part with height Hc
 *      *        * |
 *       *      *  |
 *         *__*    | - elliptic part with height Hb
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

bool gaussableCDF(std::vector<double>& val, double J0, double l_l){
    if (J0 < 0.33555 || J0 > 0.97) return true;
    //p(x) = 1/J0 - c * e^((x-0.5)^2 / w^2)
    //p(0.5) = 1/(3*J0)
    //int p(x) dx from 0 to 1 = 1
    double c = 2 / (3 * J0);
    double w = 2.5;
    double sqrt_pi = sqrt(M_PI);
    {
        double q = (1 - J0) * 3 / sqrt_pi;
        double m = 0, M = 10, tol = 1e-7;
        while (M - m > tol){
            double x = (m+M)/2;
            double f = x * erf(1/x);
            if (f > q) M = x;
            else m = x;
        }
        w = (M + m) / 4;
    }
    double erf_1d2w = erf(1/(2*w));
    for (auto& x: val)
        x = x / J0 - c * sqrt_pi / 2 * w * (erf((2*x - 1) / (2*w)) + erf_1d2w);

    return true;
}

bool onethirdCDF(std::vector<double>& val, double J0, double l_l){
    if (J0 < 0.34 || J0 > 0.97) return true;
    double d = 3*(1 - J0)/4;
    for (auto& x: val)
        if (x < 0.5 - d) x = 1/J0*x;
        else if (x < 0.5 + d) x = (1 - 2*d + x) / (3 * J0);
        else x = (x - 4*d / 3) / J0;
    return true;
}

void rotateToNewBasis(Object3D* leaf, Vector f0, Vector f1){
    f0 /= sqrt(f0.squared_length());
    f1 /= sqrt(f1.squared_length());

    Point p[3];
    for (auto v: leaf->m_mesh.vertices()){
        if (leaf->m_boundary[v] == (2 | 1) )
            p[0] = leaf->m_x[v];
        else if (leaf->m_boundary[v] == (2 | 4) )
            p[1] = leaf->m_x[v];
    }
    Vector e0 = p[1] - p[0];
    e0 /= sqrt(e0.squared_length());
    //here we invert mesh orientation
    Vector e1(e0[1], -e0[0], 0);
    for (auto v: leaf->m_mesh.vertices()){
        leaf->m_x[v] = Point((leaf->m_x[v]-CGAL::ORIGIN)*e0, (leaf->m_x[v]-CGAL::ORIGIN)*e1, 0);
    }

    p[2] = *leaf->m_x.begin();
    for (auto v: leaf->m_mesh.vertices()){
        if (leaf->m_boundary[v] == (2 | 1) )
            p[0] = leaf->m_x[v];
        else if (leaf->m_boundary[v] == (2 | 4) )
            p[1] = leaf->m_x[v];
        if (p[2].y() > leaf->m_x[v].y())
            p[2] = leaf->m_x[v];
    }

    Point c = p[0] + (p[1] - p[0])/2;
    for (auto v: leaf->m_mesh.vertices()){
        leaf->m_x[v] = CGAL::ORIGIN + (leaf->m_x[v][0] - c[0])* f0 + (leaf->m_x[v][1] - c[1]) * f1;;
    }
}

struct MeshD{
    double mesh_size; //mm

    //D, h, w, alpha, beta, s, m
    OzakiTempl templ;
    friend std::ostream& operator<<(std::ostream& out, const MeshD& t);
};
std::ostream& operator<<(std::ostream& out, const MeshD& t){
    return out << "{\n"
                  "\tmesh_size = " << t.mesh_size << "\n"
                  "\tOzakiTemplate{\n"
                  "\t\tD = " << t.templ.D << "\n"
                  "\t\th = " << t.templ.h << "\n"
                  "\t\tw = " << t.templ.w << "\n"
                  "\t\talpha = " << t.templ.alpha << "\n"
                  "\t\tbeta = " << t.templ.beta << "\n"
                  "\t\ts = " << t.templ.s << "\n"
                  "\t\tm = " << t.templ.m << "\n"
                  "\t}\n"
                  "}\n";
}

struct MaterialD{
    double thickness; //mm
    int MatType;
    double E; //kPa
    double nu;
    std::vector<double> other;
    bool with_bending;
    int MatType_bend = -1;
    double E_bend = -1, nu_bend = -1, thickness_bend = -1;
    double density = 1; //mg / mm^3
    friend std::ostream& operator<<(std::ostream& out, const MaterialD& t);
};
std::ostream& operator<<(std::ostream& out, const MaterialD& t){
    out << "{\n"
           "\tmat_type = " << t.MatType << "\n"
           "\tthickness = " << t.thickness << " mm\n"
           "\tE = " << t.E << " kPa\n"
           "\tnu = " << t.nu << "\n"
           "\trho = " << t.density << " mg/mm^3\n";
    if (t.other.size()){
        out << "\tother = { " << t.other[0];
        for (int i = 1; i < t.other.size(); ++i)
            out << ", " << t.other[i];
        out << " }\n";
    }
    out << "\twith_bending = " << (t.with_bending ? "true" : "false") << "\n";
    if (t.with_bending){
        out << "\tthickness_bend = " << t.thickness_bend << " mm\n"
               "\tmat_type_bend = " << t.MatType << "\n"
               "\tE_bend = " << t.E_bend << " kPa\n"
               "\tnu_bend = " << t.nu_bend << "\n";
    }
    return out << "}\n";
}

struct RelaxSolverD{
    double delta;
    double err;
    int maxits;
    int pfreq;
    friend std::ostream& operator<<(std::ostream& out, const RelaxSolverD& t);
};
std::ostream& operator<<(std::ostream& out, const RelaxSolverD& t){
    return out << "{\n"
                  "\tdelta = " << t.delta << "\n"
                  "\terr = " << t.err << "\n"
                  "\tmaxits = " << t.maxits << "\n"
                  "\tpfreq = " << t.pfreq << "\n"
                  "}\n";
}
struct ExternalForcesD{
    double pressure; //mm Hg

    bool with_free_edge_penalty;
    std::vector<double> free_edge_penalty_params;

    bool with_mass_force;
    std::array<double, 3> mass_force_dir;
    std::vector<double> mass_force_params;

    bool with_collision_walls;

    friend std::ostream& operator<<(std::ostream& out, const ExternalForcesD& t);
};

std::ostream& operator<<(std::ostream& out, const ExternalForcesD& t){
    auto print_array = [&out](std::string name, auto& dat){
        out << "\t" << name << " = { ";
        if (dat.size()) out << dat[0];
        for (int i = 1; i < dat.size(); ++i)
            out << ", " << dat[i];
        out << " }\n";
    };
    out << "{\n"
           "\tpressure = " << t.pressure << " mm Hg\n";
    out << "\twith_free_edge_penalty = " << (t.with_free_edge_penalty ? "true" : "false") << "\n";
    if (t.with_free_edge_penalty) print_array("free_edge_penalty_params", t.free_edge_penalty_params);
    out << "\twith_mass_force = " << (t.with_mass_force ? "true" : "false") << "\n";
    if (t.with_mass_force) print_array("mass_force_dir", t.mass_force_dir);
    if (t.with_mass_force) print_array("mass_force_params", t.mass_force_params);
    out << "\twith_collision_walls = " << (t.with_collision_walls ? "true" : "false") << "\n";
    return out << "}\n";
}

map<E_ind, BendingForce::BoundaryData::Entry>
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

Mesh::Property_map<V_ind, Vector> addBdataTag(Object3D& leaf, const map<E_ind, BendingForce::BoundaryData::Entry>& bdata){
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

namespace po = boost::program_options;

// struct SewingParam{
//     bool is_bezier = true;
//     MeshD mhd;
//     MaterialD mtd;
//     RelaxSolverD rsd;
//     ExternalForcesD efd;
//     bool with_gui_window = false;
//     bool save_steps  = false;
//     int read_step = -1;
//     std::string root_directory;
//     std::string result_directory;
//     std::string read_directory;
//     double Kcil = 0; //init cilidric curvature Rcil = 1/Kcil
//     double R = 12;
//     double init_shift = 1, shift_coef = 0.1;
//     std::array<DReal, 3> phi_all, leaf_lengths;
//     std::array<DReal, 6> phi;
//     DReal dphi = 0.2;
//     double Hc = 1.5, Hb = R;
//     double add_dist = 0;
//     int bnd_type = 1;
//     struct CommissureGeom{
//         std::vector<Point> m_A; ///<commissure points
//         Point m_O; ///< circumcenter
//         Point m_M; ///< midcenter
//         Vector m_n; ///< blood flow direction
//         DReal m_R; ///< circumradius
//         DReal m_area; ///< area of commissure triangle
//         std::array<DReal, 3> m_a; ///< sides of commissure triangle

//         void setup(){
//             auto& A = m_A;
//             m_n = CGAL::cross_product(A[1] - A[0], A[2] - A[0]);
//             m_area = sqrt(m_n.squared_length()) / 2;
//             m_n /= (2*m_area);
//             std::array<DReal, 3> a;
//             for (int i = 0; i < 3; ++i) a[i] = (A[(i+1)%3] - A[(i+2)%3]).squared_length();
//             std::array<DReal, 3> l;
//             DReal s = 0;
//             for (int i = 0; i < 3; ++i){
//                 l[i] = a[i] * (a[(i+1)%3] + a[(i+2)%3] - a[i]);
//                 s += l[i];
//             }
//             for (int i = 0; i < 3; ++i) m_a[i] = sqrt(a[i]);
//             for (int i = 0; i < 3; ++i) l[i] /= s;
//             m_O = CGAL::ORIGIN + l[0] * (A[0]-CGAL::ORIGIN) + l[1] * (A[1]-CGAL::ORIGIN) + l[2] * (A[2]-CGAL::ORIGIN);
//             m_M = CGAL::ORIGIN + ((A[0]-CGAL::ORIGIN) + (A[1]-CGAL::ORIGIN) + (A[2]-CGAL::ORIGIN))/3;
//             m_R = sqrt((A[0] - m_O).squared_length());
//         }
//     };
//     CommissureGeom m_commissure;
//     // std::array<Point, 3> m_commissure_points;
//     // Vector m_commissure_normal;
//     DReal phi0;
//     std::function<std::array<double, 3>(double t)> line;

//     SewingParam(){
//         mhd.templ.s = 0;
//         mhd.templ.m = 0;
//         mhd.templ.w = 0;
//         mhd.templ.h = 11.3;
//         mhd.templ.D = 9.5;
//         mhd.templ.alpha = 12.85 - 9.5;
//         mhd.templ.beta = 2.1;
//     }
//     void printInput(std::ostream& out = std::cout){
//         out << mhd << mtd << rsd << efd;
//         out << "with_gui_window = " << (with_gui_window ? "true" : "false") << "\n";
//         out << "with_save_steps = " << (save_steps ? "true" : "false") << "\n";
//         out << "read_step = " << (read_step < 0 ? "false" : to_string(read_step)) << "\n";
//         out << "root_directory = " << root_directory << "\n"
//                   << "result_directory = " << result_directory << "\n"
//                   << "read_directory = " << read_directory << "\n";
//         out << "{\n"
//                      "\tR = " << R << "\n"
//                      "\tphi_sizes = {" << phi_all[0] << ", " << phi_all[1] << ", " << phi_all[2] << "}\n"
//                      "\tphi_pos = {[" << phi[0] << ", " << phi[1] << "], [" << phi[2] << ", " << phi[3] << "], [" << phi[4] << ", " << phi[5] << "]}\n" 
//                      "\tHc = " << Hc << "\n"
//                      "\tHb = " << Hb << "\n"
//                      "}\n";
//         out << "add_dist = " << add_dist << "\n";
//     }
//     void parseArgs(int argc, char* argv[]){
//         try {
//             // Declare a group of options that will be
//             // allowed only on command line
//             po::options_description generic("Generic options");
//             generic.add_options()
//                     ("version,v", "print version string")
//                     ("help,h", "produce help message")
//                     ;
//             // Declare a group of options that will be
//             // allowed both on command line and in
//             // config file
//             po::options_description config("Configuration");
//             config.add_options()
//                     ("thickness,t", po::value<double>(&mtd.thickness)->default_value(0.5), "Leafs thickness, mm")
//                     ("material_model,m", po::value<int>(&mtd.MatType)->default_value(0), "Choosing of the material model: 0 - SVKNoPoisson, 1 - SVK, 2 - NHK, 3 - LIN")
//                     ("young_modulus,E", po::value<double>(&mtd.E)->default_value(1000), "Young modulus, kPa")
//                     ("poisson_modulus,N", po::value<double>(&mtd.nu)->default_value(0), "Poisson modulus")
//                     ("init_cilindric", po::value<double>(&Kcil)->default_value(0), "Initial cilidric curvature Rcil = 1/Kcil, mm^{-1}")
//                     ("with_bending,b", "Take into account bending force")
//                     ("with_planes", "Use planes insted bezier curves")
//                     ("tb", po::value<double>(&mtd.thickness_bend)->default_value(-1), "Leafs thickness for bending force, mm")
//                     ("mb", po::value<int>(&mtd.MatType_bend)->default_value(-1), "Choosing of the material model for bending force: 0 - SVKNoPoisson, 1 - SVK, 2 - NHK")
//                     ("Eb", po::value<double>(&mtd.E_bend)->default_value(-1), "Young modulus for bending force, kPa")
//                     ("Nb", po::value<double>(&mtd.nu_bend)->default_value(-1), "Poisson modulus for bending force")
//                     ("mesh_step", po::value<double>(&mhd.mesh_size)->default_value(1), "Step of generated mesh, mm")
//                     ("template_geometry,g", po::value<std::vector<double>>()->multitoken()->default_value({19, 11.3, 0, 3.35, 2.1, 0, 0}), "Parameters of template geometry {D, h, w, alpha, beta, s, m}")
//                     ("delta,d", po::value<double>(&rsd.delta)->default_value(5e-8), "Set relaxation solver iteration parameter")
//                     ("err,e", po::value<double>(&rsd.err)->default_value(1e-7), "Set stop criterion based on relative error decrease")
//                     ("maxit,i", po::value<int>(&rsd.maxits)->default_value(1e5), "Set stop criterion based on maximum number of iterations")
//                     ("pfreq,q", po::value<int>(&rsd.pfreq)->default_value(1e3), "Frequency of printing of state")
//                     ("pressure,p", po::value<double>(&efd.pressure)->default_value(80), "Diastolic pressure, Hg mm")
//                     ("with_free_edge_penalty", po::value<std::vector<double>>(&efd.free_edge_penalty_params), "Turn on edge penalty")
//                     ("with_mass_force", po::value<std::vector<double>>()->multitoken(), "Turn on mass force arg{direction[3], params[*]]}")
//                     ("with_collision_walls", "Turn on collision walls")
//                     ("root_directory", po::value<std::string>(&root_directory)->default_value("../../../result/Sewing"), "Root directory for saving and reading of results")
//                     ("save_directory", po::value<std::string>(&result_directory)->default_value("default"), "Path relative root directory to save steps")
//                     ("read_directory", po::value<std::string>(&read_directory)->default_value("default"), "Path relative root directory to read start step")
//                     ("view", "Create ImGui debug view window")
//                     ("save_steps", "Save sewing steps")
//                     ("read_step", po::value<int>(&read_step)->default_value(-1), "Start from specific previously saved step")
//                     ("bnd.type", po::value<int>(&bnd_type)->default_value(1), "Type of sew line (1 - simple, 2 - complex)")
//                     ("bnd.R", po::value<double>(&R)->default_value(12), "Radius of cilinder, mm")
//                     ("bnd.phi", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{1, 1, 1}), "Relations[3] of angle size of aortic suture lines")
//                     ("bnd.dphi", po::value<double>(&dphi)->default_value(0.2), "Angle-size of sutured lines distance, radians")
//                     ("bnd.Hc", po::value<double>(&Hc)->default_value(1.5), "Vertical allowance of the line, mm")
//                     ("suture_relations", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{1, 1, 1}), "Relation of free lengths of leaflets")
//                     ("bnd.Hb", po::value<double>(&Hb)->default_value(12), "Deflection depth of the line, mm")
//                     ("add_dist", po::value<double>(&add_dist)->default_value(0), "Additional distance for template from line, mm")
//                     ;

//             po::options_description cmdline_options;
//             cmdline_options.add(generic).add(config);

//             po::variables_map vm;
//             po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
//             po::notify(vm);

//             if (vm.count("help")) {
//                 cout << cmdline_options << "\n";
//                 exit(0);
//             }
//             if (vm.count("version")){
//                 cout << "SewOnCilinder version 0.1.0\n";
//                 exit(0);
//             }
//             {
//                 auto& _phi = vm["bnd.phi"].as<std::vector<double>>();
//                 std::copy(_phi.data(), _phi.data() + 3, phi_all.data());
//                 auto& l = vm["suture_relations"].as<std::vector<double>>();
//                 std::copy(l.begin(), l.end(), leaf_lengths.data());
//             }
//             with_gui_window = vm.count("view");
//             save_steps = vm.count("save_steps");
//             mtd.with_bending = vm.count("with_bending");
//             is_bezier = !vm.count("with_planes");
//             efd.with_free_edge_penalty = vm.count("with_free_edge_penalty");
//             efd.with_mass_force = vm.count("with_mass_force");
//             efd.with_collision_walls = vm.count("with_collision_walls");
//             auto tg = vm["template_geometry"].as<std::vector<double>>();
//             auto& ot = mhd.templ;
//     #define IN_T(X, Y) if (tg.size() > X) ot.Y = tg[X]
//             IN_T(0, D); IN_T(1, h); IN_T(2, w); IN_T(3, alpha); IN_T(4, beta); IN_T(5, s); IN_T(6, m);
//     #undef IN_T
//             if (efd.with_mass_force) {
//                 auto mf = vm["with_mass_force"].as<std::vector<double>>();
//                 if (mf.size() < 3) throw std::runtime_error("Mass force should define direction array that has 3 values");
//                 std::copy(mf.data(), mf.data() + 3, efd.mass_force_dir.data());
//                 efd.mass_force_params.resize(mf.size() - 3);
//                 std::copy(mf.data() + 3, mf.data() + mf.size(), efd.mass_force_params.data());
//             }
//             if (mtd.with_bending){
//                 if (mtd.E_bend < 0) mtd.E_bend = mtd.E;
//                 if (mtd.nu_bend < 0) mtd.nu_bend = mtd.nu;
//                 if (mtd.thickness_bend < 0) mtd.thickness_bend = mtd.thickness;
//                 if (mtd.MatType_bend < 0) mtd.MatType_bend = mtd.MatType_bend;
//             }
//             if (efd.with_free_edge_penalty){
//                 if (efd.free_edge_penalty_params.empty())
//                     efd.free_edge_penalty_params.push_back(0.1);
//             }
//             if (save_steps || read_step >= 0){
//                 path p(root_directory);
//                 if (!exists(status(p)) || !is_directory(status(p)))
//                     throw std::runtime_error("root_directory = \"" + root_directory + "\" is not exists");
//                 if (read_step >= 0){
//                     p = root_directory + "/" + read_directory;
//                     if (!exists(status(p)) || !is_directory(status(p)))
//                         throw std::runtime_error("read_directory = \"" + p.string() + "\" is not exists");
//                     p = root_directory + "/" + read_directory + "/" + "sew_step_" + std::to_string(read_step) + "_x.txt";
//                     if (!exists(status(p)) || !is_regular_file(status(p)))
//                         throw std::runtime_error("step file = \"" + p.string() + "\" is not exists");
//                     p = root_directory + "/" + read_directory + "/" + "sew_step_" + std::to_string(read_step) + "_x0.tag";
//                     if (!exists(status(p)) || !is_regular_file(status(p)))
//                         throw std::runtime_error("step file = \"" + p.string() + "\" is not exists");
//                     p = root_directory + "/" + read_directory + "/" + "sew_step_" + std::to_string(read_step) + "_init_x.tag";
//                     if (!exists(status(p)) || !is_regular_file(status(p)))
//                         throw std::runtime_error("step file = \"" + p.string() + "\" is not exists");
//                 }
//                 if (save_steps){
//                     p = root_directory + "/" + result_directory;
//                     if (!exists(status(p))) create_directory(p);
//                     else if (!is_directory(status(p)))
//                         throw std::runtime_error("save_directory = \"" + p.string() + "\" but it's exists as nondirectory");
//                 }
//             }
//             read_directory = root_directory + "/" + read_directory;
//             result_directory = root_directory + "/" + result_directory;
//         }
//         catch(exception& e) {
//             cerr << "error: " << e.what() << "\n";
//             exit(-1);
//         }
//         catch(...) {
//             cerr << "Exception of unknown type!\n";
//             exit(-2);
//         }
        
//         {
//             auto& phi_relations = phi_all;
//             DReal fphi = 2*M_PI - 3*dphi;
//             DReal c = phi_relations[0] + phi_relations[1] + phi_relations[2];
//             for (int i = 0; i < 3; ++i) phi_relations[i] *= fphi / c;
//             phi[0] = -phi_relations[0]/2; phi[1] = phi_relations[0]/2;
//             phi[2] = phi[1] + dphi; phi[3] = phi[2] + phi_relations[1];
//             phi[4] = phi[3] + dphi; phi[5] = phi[4] + phi_relations[2];
//         }

//         {
//             DReal L0 = mhd.templ.D + 2*mhd.templ.alpha;
//             for (int i = 2; i >= 0; --i) leaf_lengths[i] = L0 * leaf_lengths[i] / leaf_lengths[0];
//             for (int i = 0; i < 3; ++i){
//                 DReal phi0 = phi[2*((i+1)%3)+1];
//                 if (phi0 < 0) phi0 += 2*M_PI;
//                 DReal phi1 = phi[2*((i+2)%3)];
//                 if (phi1 < 0) phi1 += 2*M_PI;
//                 DReal lphi = (phi0 + phi1)/2;
//                 m_commissure.m_A[i] = Point(R*cos(lphi), R*sin(lphi), Hc);
//             }
//             m_commissure.setup();
//         }
//         printInput(std::cout);
//     }
//     static SewingParam makeSewingParam(int argc, char* argv[]){
//         SewingParam res;
//         res.parseArgs(argc, argv);
//         return res;
//     }
//     Mesh prepareMesh(){
//         Mesh m;
//         if (read_step < 0) {
//             AniMesh am = generate_ozaki(mhd.templ, mhd.mesh_size);
//             m = convert_to_Mesh(am, "v:boundary_lbl");
//         } else {
//             Object3D _obj;
//             _obj.read(read_directory + "/sew_step_" + std::to_string(read_step) + "_x.txt");
//             m = _obj.m_mesh;
//             ObReadTag(m, read_directory + "/sew_step_" + std::to_string(read_step) + "_x0.tag", "v:point");
//             ObReadTag(m, read_directory + "/sew_step_" + std::to_string(read_step) + "_init_x.tag", "init_x");
//         }
//         return m;
//     }
//     static void align_to_vectors(Mesh& m, Vector right, Vector up, std::string bnd_tag_name = "v:boundary_lbl"){
//         Vector f0 = right / sqrt(right.squared_length()),
//                f1 = up / sqrt(up.squared_length());
//         auto& m_boundary = m.property_map<V_ind, int>(bnd_tag_name).first;
//         auto& m_x = m.property_map<V_ind, Point>("v:point").first;
//         Point p[3];
//         for (auto v: m.vertices()){
//             if (m_boundary[v] == (2 | 1) )
//                 p[0] = m_x[v];
//             else if (m_boundary[v] == (2 | 4) )
//                 p[1] = m_x[v];
//         }
//         Vector e0 = p[1] - p[0];
//         e0 /= sqrt(e0.squared_length());
//         //here we invert mesh orientation
//         Vector e1(e0[1], -e0[0], 0);
//         for (auto v: m.vertices()){
//             m_x[v] = Point((m_x[v]-CGAL::ORIGIN)*e0, (m_x[v]-CGAL::ORIGIN)*e1, 0);
//         }

//         p[2] = *m_x.begin();
//         for (auto v: m.vertices()){
//             if (m_boundary[v] == (2 | 1) )
//                 p[0] = m_x[v];
//             else if (m_boundary[v] == (2 | 4) )
//                 p[1] = m_x[v];
//             if (p[2].y() > m_x[v].y())
//                 p[2] = m_x[v];
//         }

//         Point c = p[0] + (p[1] - p[0])/2;
//         for (auto v: m.vertices()){
//             m_x[v] = CGAL::ORIGIN + (m_x[v][0] - c[0])* f0 + (m_x[v][1] - c[1]) * f1;;
//         }
//     }
//     static void align_to_vectors(Object3D& obj, Vector right, Vector up){
//         align_to_vectors(obj.m_mesh, right, up);
//         for (auto v: obj.m_mesh.vertices()){ obj.m_x[v] = obj.m_x0[v]; }
//     }
//     static void cilindrize(double Kcil, Object3D& obj, Vector n, Vector r, Point new_origin){
//         if (Kcil == 0) return;
//         auto& m = obj.m_mesh;
//         auto f2 = -r;
//         auto f0 = n;
//         auto f1 = CGAL::cross_product(f2, f0);
//         auto    x0O = (new_origin - CGAL::ORIGIN) * f0,
//                 x1O = (new_origin - CGAL::ORIGIN) * f1,
//                 x2O = (new_origin - CGAL::ORIGIN) * f2;
//         for (auto v: m.vertices()) {
//             auto r = obj.m_x[v] - CGAL::ORIGIN;
//             auto x0 = r * f0, x1 = r * f1;
//             auto dx1 = x1 - x1O;
//             auto kdx1 = Kcil * dx1;
//             auto kdx1p2 = kdx1 * kdx1;
//             auto x0_ = x0 - x0O;
//             auto x1_ = x1 + dx1*kdx1p2*(kdx1p2 - 20)/120;
//             auto x2_ = x2O + dx1*kdx1/40320*(((kdx1p2-56)*kdx1p2+1680)*kdx1p2-20160);
//             obj.m_x[v] = CGAL::ORIGIN + x1_ * f1 + x0_ * f0 + x2_ * f2;
//         }
//     }
//     bool is_obj_intersect_suture_line(Object3D& obj, int line_num) { return false; }
//     void placeMeshNearCilinderAorta(Object3D& obj, int nplace){
//         if (nplace < 0 || nplace >= 3) throw std::runtime_error("nplace should be 0, 1 or 2");
//         auto i = nplace;
//         Vector f0 = (m_commissure.m_A[(i+2)%3] - m_commissure.m_A[(i+1)%3])/m_commissure.m_a[i], f1 = m_commissure.m_n;
//         Vector f2 = CGAL::cross_product(f1, f0);
//         align_to_vectors(obj, f0, f1);
//         Point po = m_commissure.m_A[(i+2)%3] + (m_commissure.m_A[(i+1)%3] - m_commissure.m_A[(i+2)%3])/2 
//                     + f2 * init_shift * m_commissure.m_R;
//         auto init_x = obj.m_mesh.add_property_map<V_ind, Point>("v:init_x").first;
//         for (auto v: obj.m_mesh.vertices()){
//             init_x[v] = obj.m_x0[v];
//             obj.m_x0[v] = obj.m_x[v] = po + (obj.m_x[v] - CGAL::ORIGIN);
//         }
//         cilindrize(Kcil, obj, f1, -f2, po);
//         while (is_obj_intersect_suture_line(obj, i)){
//             for (auto v: obj.m_mesh.vertices()){
//                 init_x[v] = obj.m_x0[v];
//                 obj.m_x0[v] = obj.m_x[v] += f2 * shift_coef * m_commissure.m_R;
//             }
//         }
//         obj.apply_updaters();
//     }
//     void map_cusp_to_suture_line(Object3D& obj, int line_num){
//         phi0 = phi[2*line_num+1] - phi[2*line_num];
//         std::function<std::array<double, 3> (double t)> _line;
//         if (bnd_type == 1)
//             _line = get_line1(R, phi0, Hc, Hb);
//         else if (bnd_type == 2)
//             _line = get_line2(R, phi0, Hc, Hb);
//         else
//             throw std::runtime_error("Faced unknown sew line type");
//         DReal dphi = (phi[2*line_num+1] + phi[2*line_num])/2;

//         line = [_line, dphi](double t)->std::array<double, 3> {
//             auto x = _line(t);
//             auto phi_x = atan2(x[1], x[0]);
//             auto R = hypot(x[0], x[1]);
//             x[0] = R*cos(dphi + phi_x);
//             x[1] = R*sin(dphi + phi_x);
//             return {x[0], x[1], x[2]};
//         };

//     }
// };

int make_sewing(int argc, char* argv[]){
    bool is_bezier = true;
    MeshD mhd;
    MaterialD mtd;
    RelaxSolverD rsd;
    ExternalForcesD efd;
    bool with_gui_window = false;
    bool save_steps  = false;
    int read_step = -1;
    std::string root_directory;
    std::string result_directory;
    std::string read_directory;
    mhd.templ.s = 0;
    mhd.templ.m = 0;
    mhd.templ.w = 0;
    mhd.templ.h = 11.3;
    mhd.templ.D = 9.5;
    mhd.templ.alpha = 12.85 - 9.5;
    mhd.templ.beta = 2.1;

    double Kcil = 0; //init cilidric curvature Rcil = 1/Kcil
    double R = 12;
    //double phi0 = 2*M_PI/3 - 0.4;
    std::array<DReal, 3> phi_all, leaf_lengths;
    std::array<DReal, 6> phi;
    DReal dphi = 0.2;
    double Hc = 1.5, Hb = R;
    double add_dist = 0;
    int bnd_type = 1;
    auto print_input = [&](){
        std::cout << mhd << mtd << rsd << efd;
        std::cout << "with_gui_window = " << (with_gui_window ? "true" : "false") << "\n";
        std::cout << "with_save_steps = " << (save_steps ? "true" : "false") << "\n";
        std::cout << "read_step = " << (read_step < 0 ? "false" : to_string(read_step)) << "\n";
        std::cout << "root_directory = " << root_directory << "\n"
                  << "result_directory = " << result_directory << "\n"
                  << "read_directory = " << read_directory << "\n";
        std::cout << "{\n"
                     "\tR = " << R << "\n"
                     "\tphi_sizes = {" << phi_all[0] << ", " << phi_all[1] << ", " << phi_all[2] << "}\n"
                     "\tphi_pos = {[" << phi[0] << ", " << phi[1] << "], [" << phi[2] << ", " << phi[3] << "], [" << phi[4] << ", " << phi[5] << "]}\n" 
                     "\tHc = " << Hc << "\n"
                     "\tHb = " << Hb << "\n"
                     "}\n";
        std::cout << "add_dist = " << add_dist << "\n";
    };

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
                ("thickness,t", po::value<double>(&mtd.thickness)->default_value(0.5), "Leafs thickness, mm")
                ("material_model,m", po::value<int>(&mtd.MatType)->default_value(0), "Choosing of the material model: 0 - SVKNoPoisson, 1 - SVK, 2 - NHK, 3 - LIN")
                ("young_modulus,E", po::value<double>(&mtd.E)->default_value(1000), "Young modulus, kPa")
                ("poisson_modulus,N", po::value<double>(&mtd.nu)->default_value(0), "Poisson modulus")
                ("init_cilindric", po::value<double>(&Kcil)->default_value(0), "Initial cilidric curvature Rcil = 1/Kcil, mm^{-1}")
                ("with_bending,b", "Take into account bending force")
                ("with_planes", "Use planes insted bezier curves")
                ("tb", po::value<double>(&mtd.thickness_bend)->default_value(-1), "Leafs thickness for bending force, mm")
                ("mb", po::value<int>(&mtd.MatType_bend)->default_value(-1), "Choosing of the material model for bending force: 0 - SVKNoPoisson, 1 - SVK, 2 - NHK")
                ("Eb", po::value<double>(&mtd.E_bend)->default_value(-1), "Young modulus for bending force, kPa")
                ("Nb", po::value<double>(&mtd.nu_bend)->default_value(-1), "Poisson modulus for bending force")
                ("mesh_step", po::value<double>(&mhd.mesh_size)->default_value(1), "Step of generated mesh, mm")
                ("template_geometry,g", po::value<std::vector<double>>()->multitoken()->default_value({19, 11.3, 0, 3.35, 2.1, 0, 0}), "Parameters of template geometry {D, h, w, alpha, beta, s, m}")
                ("delta,d", po::value<double>(&rsd.delta)->default_value(5e-8), "Set relaxation solver iteration parameter")
                ("err,e", po::value<double>(&rsd.err)->default_value(1e-7), "Set stop criterion based on relative error decrease")
                ("maxit,i", po::value<int>(&rsd.maxits)->default_value(1e5), "Set stop criterion based on maximum number of iterations")
                ("pfreq,q", po::value<int>(&rsd.pfreq)->default_value(1e3), "Frequency of printing of state")
                ("pressure,p", po::value<double>(&efd.pressure)->default_value(80), "Diastolic pressure, Hg mm")
                ("with_free_edge_penalty", po::value<std::vector<double>>(&efd.free_edge_penalty_params), "Turn on edge penalty")
                ("with_mass_force", po::value<std::vector<double>>()->multitoken(), "Turn on mass force arg{direction[3], params[*]]}")
                ("with_collision_walls", "Turn on collision walls")
                ("root_directory", po::value<std::string>(&root_directory)->default_value("../../../result/Sewing"), "Root directory for saving and reading of results")
                ("save_directory", po::value<std::string>(&result_directory)->default_value("default"), "Path relative root directory to save steps")
                ("read_directory", po::value<std::string>(&read_directory)->default_value("default"), "Path relative root directory to read start step")
                ("view", "Create ImGui debug view window")
                ("save_steps", "Save sewing steps")
                ("read_step", po::value<int>(&read_step)->default_value(-1), "Start from specific previously saved step")
                ("bnd.type", po::value<int>(&bnd_type)->default_value(1), "Type of sew line (1 - simple, 2 - complex)")
                ("bnd.R", po::value<double>(&R)->default_value(12), "Radius of cilinder, mm")
                ("bnd.phi", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{1, 1, 1}), "Relations[3] of angle size of aortic suture lines")
                ("bnd.dphi", po::value<double>(&dphi)->default_value(0.2), "Angle-size of sutured lines distance, radians")
                ("bnd.Hc", po::value<double>(&Hc)->default_value(1.5), "Vertical allowance of the line, mm")
                ("suture_relations", po::value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{1, 1, 1}), "Relation of free lengths of leaflets")
                ("bnd.Hb", po::value<double>(&Hb)->default_value(12), "Deflection depth of the line, mm")
                ("add_dist", po::value<double>(&add_dist)->default_value(0), "Additional distance for template from line, mm")
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
            cout << "SewOnCilinder version 0.1.0\n";
            return 0;
        }
        {
            auto& _phi = vm["bnd.phi"].as<std::vector<double>>();
            std::copy(_phi.data(), _phi.data() + 3, phi_all.data());
            auto& l = vm["suture_relations"].as<std::vector<double>>();
            std::copy(l.begin(), l.end(), leaf_lengths.data());
        }
        with_gui_window = vm.count("view");
        save_steps = vm.count("save_steps");
        mtd.with_bending = vm.count("with_bending");
        is_bezier = !vm.count("with_planes");
        efd.with_free_edge_penalty = vm.count("with_free_edge_penalty");
        efd.with_mass_force = vm.count("with_mass_force");
        efd.with_collision_walls = vm.count("with_collision_walls");
        auto tg = vm["template_geometry"].as<std::vector<double>>();
        auto& ot = mhd.templ;
#define IN_T(X, Y) if (tg.size() > X) ot.Y = tg[X]
        IN_T(0, D); IN_T(1, h); IN_T(2, w); IN_T(3, alpha); IN_T(4, beta); IN_T(5, s); IN_T(6, m);
#undef IN_T
        if (efd.with_mass_force) {
            auto mf = vm["with_mass_force"].as<std::vector<double>>();
            if (mf.size() < 3) throw std::runtime_error("Mass force should define direction array that has 3 values");
            std::copy(mf.data(), mf.data() + 3, efd.mass_force_dir.data());
            efd.mass_force_params.resize(mf.size() - 3);
            std::copy(mf.data() + 3, mf.data() + mf.size(), efd.mass_force_params.data());
        }
        if (mtd.with_bending){
            if (mtd.E_bend < 0) mtd.E_bend = mtd.E;
            if (mtd.nu_bend < 0) mtd.nu_bend = mtd.nu;
            if (mtd.thickness_bend < 0) mtd.thickness_bend = mtd.thickness;
            if (mtd.MatType_bend < 0) mtd.MatType_bend = mtd.MatType_bend;
        }
        if (efd.with_free_edge_penalty){
            if (efd.free_edge_penalty_params.empty())
                efd.free_edge_penalty_params.push_back(0.1);
        }
        if (save_steps || read_step >= 0){
            path p(root_directory);
            if (!exists(status(p)) || !is_directory(status(p)))
                throw std::runtime_error("root_directory = \"" + root_directory + "\" is not exists");
            if (read_step >= 0){
                p = root_directory + "/" + read_directory;
                if (!exists(status(p)) || !is_directory(status(p)))
                    throw std::runtime_error("read_directory = \"" + p.string() + "\" is not exists");
                p = root_directory + "/" + read_directory + "/" + "sew_step_" + std::to_string(read_step) + "_x.txt";
                if (!exists(status(p)) || !is_regular_file(status(p)))
                    throw std::runtime_error("step file = \"" + p.string() + "\" is not exists");
                p = root_directory + "/" + read_directory + "/" + "sew_step_" + std::to_string(read_step) + "_x0.tag";
                if (!exists(status(p)) || !is_regular_file(status(p)))
                    throw std::runtime_error("step file = \"" + p.string() + "\" is not exists");
                p = root_directory + "/" + read_directory + "/" + "sew_step_" + std::to_string(read_step) + "_init_x.tag";
                if (!exists(status(p)) || !is_regular_file(status(p)))
                    throw std::runtime_error("step file = \"" + p.string() + "\" is not exists");
            }
            if (save_steps){
                p = root_directory + "/" + result_directory;
                if (!exists(status(p))) create_directory(p);
                else if (!is_directory(status(p)))
                    throw std::runtime_error("save_directory = \"" + p.string() + "\" but it's exists as nondirectory");
            }
        }
        read_directory = root_directory + "/" + read_directory;
        result_directory = root_directory + "/" + result_directory;
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return -1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        return -2;
    }
    
    {
        auto& phi_relations = phi_all;
        DReal fphi = 2*M_PI - 3*dphi;
        DReal c = phi_relations[0] + phi_relations[1] + phi_relations[2];
        for (int i = 0; i < 3; ++i) phi_relations[i] *= fphi / c;
        phi[0] = -phi_relations[0]/2; phi[1] = phi_relations[0]/2;
        phi[2] = phi[1] + dphi; phi[3] = phi[2] + phi_relations[1];
        phi[4] = phi[3] + dphi; phi[5] = phi[4] + phi_relations[2];
    }

    std::array<Point, 3> m_commissure_points;
    {
        DReal L0 = mhd.templ.D + 2*mhd.templ.alpha;
        for (int i = 2; i >= 0; --i) leaf_lengths[i] = L0 * leaf_lengths[i] / leaf_lengths[0];
        for (int i = 0; i < 3; ++i){
            DReal phi0 = phi[2*((i+1)%3)+1];
            if (phi0 < 0) phi0 += 2*M_PI;
            DReal phi1 = phi[2*((i+2)%3)];
            if (phi1 < 0) phi1 += 2*M_PI;
            DReal lphi = (phi0 + phi1)/2;
            m_commissure_points[i] = Point(R*cos(lphi), R*sin(lphi), Hc);
        }
    }

    DReal phi0 = phi[1] - phi[0];
    print_input();
    std::function<std::array<double, 3>(double t)> line;
    if (bnd_type == 1)
        line = get_line1(R, phi0, Hc, Hb);
    else if (bnd_type == 2)
        line = get_line2(R, phi0, Hc, Hb);
    else
        throw std::runtime_error("Faced unknown sew line type");

    Object3D _obj;
    if (read_step < 0) {
        AniMesh am = generate_ozaki(mhd.templ, mhd.mesh_size);
        _obj = Object3D(convert_to_Mesh(am, "v:boundary_lbl"));
    } else {
        _obj.read(read_directory + "/sew_step_" + std::to_string(read_step) + "_x.txt");
        ObReadTag(_obj, read_directory + "/sew_step_" + std::to_string(read_step) + "_x0.tag", "v:point");
        ObReadTag(_obj, read_directory + "/sew_step_" + std::to_string(read_step) + "_init_x.tag", "init_x");
    }

    AVSimulator model;
    std::thread* front_end;
    if (with_gui_window) {
        front_end = new std::thread(Listener, argc, argv);
        auto renderer = std::make_unique<World3d::DefaultRenderer>();
        renderer->set_interconnector(&g_interconnection);
        model.m_w.setRenderer(std::move(renderer));
    }
    QuickLineIntegral qli;
    qli.m_curve = [line](double t) { auto p = line(t); return Point{p[0], p[1], p[2]}; };
    qli.setup();
    model.sld.commissure_points = std::array<Point, 3>{qli.m_curve(0), qli.m_curve(1), Point{-R, 0, 0}};
    model.sld.lines.resize(1);
    model.sld.lines[0].resize(qli.quick_index.size());
    for (int i = 0; i < model.sld.lines[0].size(); ++i) model.sld.lines[0][i] = qli.quick_index[i].p - CGAL::ORIGIN;
    auto &com_p = model.sld.commissure_points;
    model.sld.orientation = CGAL::cross_product(com_p[1] - com_p[0], com_p[2] - com_p[0]);
    model.sld.orientation /= sqrt(model.sld.orientation.squared_length());
    model.sld.cdfs.emplace_back(onethirdCDF);//gaussableCDF);
    auto lid = model.m_w.addObject3D(std::move(_obj), 1);
    auto& obj = model.m_w.obj(lid);
    auto save_new_step = [f = result_directory, &obj, save_steps](int nstep){
        if (!save_steps) return;
        obj.save(f + "/sew_step_" + to_string(nstep) + "_x.txt");
        obj.save(f + "/sew_step_" + to_string(nstep) + ".vtk");
        ObSaveTag(obj, f + "/sew_step_" + to_string(nstep) + "_x0.tag", "v:point");
        ObSaveTag(obj, f + "/sew_step_" + to_string(nstep) + "_init_x.tag", "init_x");
    };

    if (read_step < 0) {
        auto f0 = model.sld.orientation;
        auto f1 = model.sld.commissure_points[1] - model.sld.commissure_points[0];
        f1 /= sqrt(f1.squared_length());
        auto midcenter = ((com_p[2] - CGAL::ORIGIN) + (com_p[1] - CGAL::ORIGIN) + (com_p[0] - CGAL::ORIGIN)) / 3;
        auto dif = midcenter - (com_p[0] - CGAL::ORIGIN);
        auto f2 = dif - (dif * f1) * f1;
        double pdist = sqrt((dif - (dif * f1) * f1).squared_length());
        f2 /= pdist;
        double min_dist_expected = pdist / 2;
        using Plane = CGAL::Plane_3<CGAL::Simple_cartesian<double>>;
        auto n = -cross_product(f0, f1);
        n /= sqrt(n.squared_length());
        Plane p(com_p[0], -cross_product(f0, f1));
        double min_dist_real = 0;
        for (auto dp: model.sld.lines[0])
            if ((CGAL::ORIGIN + dp - com_p[0]) * n < min_dist_real) min_dist_real = (CGAL::ORIGIN + dp - com_p[0]) * n;
        rotateToNewBasis(&obj, f1, f0);
        Point new_origin = com_p[0] + (com_p[1] - com_p[0]) / 2 + f2 * (min_dist_expected - min_dist_real + add_dist);
        auto init_x = obj.m_mesh.add_property_map<V_ind, Point>("init_x").first;
        for (auto v: obj.m_mesh.vertices()) {
            obj.m_x[v] = new_origin + (obj.m_x[v] - CGAL::ORIGIN);
            init_x[v] = obj.m_x0[v];
            obj.m_x0[v] = obj.m_x[v];
        }
        bool to_cilindric = (Kcil > 0);
        if (to_cilindric) {
            if (Kcil > 2/R) throw std::runtime_error("Too big initial curvature of cilindric surface");
            auto f2 = cross_product(f0, f1);
            auto    x0O = (new_origin - CGAL::ORIGIN) * f0,
                    x1O = (new_origin - CGAL::ORIGIN) * f1,
                    x2O = (new_origin - CGAL::ORIGIN) * f2;
            for (auto v: obj.m_mesh.vertices()) {
                auto r = obj.m_x[v] - CGAL::ORIGIN;
                auto x0 = r * f0, x1 = r * f1;
                auto dx1 = x1 - x1O;
                auto kdx1 = Kcil * dx1;
                auto kdx1p2 = kdx1 * kdx1;
                auto x0_ = x0 - x0O;
                auto x1_ = x1 + dx1*kdx1p2*(kdx1p2 - 20)/120;
                auto x2_ = x2O + dx1*kdx1/40320*(((kdx1p2-56)*kdx1p2+1680)*kdx1p2-20160);
                obj.m_x[v] = CGAL::ORIGIN + x1_ * f1 + x0_ * f0 + x2_ * f2;
                obj.m_x0[v] = obj.m_x[v];
            }

            // double RO = 1 / Kcil;
            // auto rO = (new_origin - CGAL::ORIGIN) - RO * f3;
            // auto x1O = rO * f1, x2O = rO * f0, x3O = rO * f3;
            // for (auto v: obj.m_mesh.vertices()) {
            //     auto r = obj.m_x[v] - CGAL::ORIGIN;
            //     auto x1 = r * f1, x2 = r * f0, x3 = r * f3;
            //     auto phi = (x1 - x1O) / RO;
            //     auto x1_ = x1O + RO * sin(phi), x2_ = x2, x3_ = x3O + RO * cos(phi);
            //     obj.m_x[v] = CGAL::ORIGIN + x1_ * f1 + x2_ * f0 + x3_ * f3;
            //     obj.m_x0[v] = obj.m_x[v]; //TODO: remove this
            // }
        }

        model.sewLeafBndToSewBnd(std::vector<Object3D*>{&obj});
        save_new_step(0);
    }

    static const int FIXED = 2, FREE = 1 | 4 | 8;
    DirichletCond dc = [](Object3D &obj, const V_ind &v) -> unsigned char {
        if (obj.m_boundary[v] & FIXED) return 0;
        return 7;
    };
    obj.setDirichletCond(dc);
    auto bdata = computeBndDataOnCilAligned(obj,
                [](Object3D& leaf, E_ind e, V_ind* v, V_ind vop){
                    return (leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED) && !(leaf.m_boundary[vop] & FREE);
                });
    auto bdit = addBdataTag(obj, bdata);
    //obj.save(result_directory + "/start_configuration.vtk");
    auto _mass = obj.m_mesh.add_property_map<V_ind, DReal>("v:mass").first;
    auto _Ap = set_Ap(obj);
    for (auto v: obj.m_mesh.vertices()){
        auto ff = face_around(obj.m_mesh, v);
        DReal s = 0;
        for (auto f: ff) s += _Ap[f];
        if (!ff.empty()){
            _mass[v] = s * mtd.density * mtd.thickness / ff.size();
        } else 
            _mass[v] = 1e20;
        if (!obj.is_movable(v)) _mass[v] = 1e20;
    }

    Force_ID pf_id;
    Force_ID ef_id, bf_id;
    auto set_force = [&](int MatType, double H, double E, double nu,
                        bool with_bending = false, int MatType_bend = -1, double H_bend = -1, double E_bend = -1, double nu_bend = -1){
        if (with_bending){
            if (MatType_bend < 0) MatType_bend = MatType;
            if (H_bend < 0) H_bend = H;
            if (E_bend < 0) E_bend = E;
            if (nu_bend < 0) nu_bend = nu;
        }
        {
            Force elastic;
            double lambda = E * nu / (1 - nu * nu), mu = E / (2 * (1 + nu));
            if (MatType == 0) {
                elastic = SVKirchhoffNoPoissonModel(H, lambda, mu, "../../../generated", false);
                elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
            } else if (MatType == 1) {
                elastic = SVKirchhoffModel(H, lambda, mu, "../../../generated", false);
                elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
            } else if (MatType == 2) {
                elastic = NeoGookModelOpt(mu, H);
            } else if (MatType == 3){
                elastic = LinearElasticModel1(H, lambda, mu, "../../../generated", false);
                elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
            } else {
                throw std::runtime_error("Unknown elastic model type = " + to_string(mtd.MatType));
            }
            ef_id = model.m_w.addForce(elastic, lid);
        }
        if (with_bending){
            Force bending;
            double lambda = E_bend * nu_bend / (1 - nu_bend * nu_bend), mu = E_bend / (2 * (1 + nu_bend));
            if (MatType_bend == 0) {
                Force proxy = SVKirchhoffNoPoissonModel(H_bend, lambda, mu, "../../../generated", false);
                proxy.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
                bending = BendingForce(proxy.target<HyperElasticForceBase>()->f, false);
            } else if (MatType_bend == 1) {
                Force proxy = SVKirchhoffModel(H_bend, lambda, mu, "../../../generated", false);
                proxy.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
                bending = BendingForce(proxy.target<HyperElasticForceBase>()->f, false);
            } else if (MatType_bend == 2) {
                Force proxy = NeoGookModel(mu, H_bend, "../../../generated", false);
                proxy.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
                bending = BendingForce(proxy.target<HyperElasticForceBase>()->f, false);
            } else if (MatType_bend == 3){
                Force proxy = LinearElasticModel1(H_bend, lambda, mu, "../../../generated", false);
                proxy.target<HyperElasticForceBase>()->prepareJacobianFunction(false);
                bending = BendingForce(proxy.target<HyperElasticForceBase>()->f, false);
            } else {
                throw std::runtime_error("Unknown elastic model type = " + to_string(mtd.MatType));
            }
            bending.target<BendingForce>()->set_boundary_data(obj, bdata);
            bending.target<BendingForce>()->prepareJacobianFunction(false);
            bf_id = model.m_w.addForce(std::move(bending), lid);
        }
    };
    set_force(mtd.MatType, mtd.thickness, mtd.E, mtd.nu, false, mtd.MatType_bend, mtd.thickness_bend, mtd.E_bend, mtd.nu_bend);
    if (efd.with_free_edge_penalty) {
        double S = 0;
        for (auto f: obj.m_mesh.faces()){
            auto v = vert_around(obj.m_mesh, f);
            S += sqrt((CGAL::cross_product(obj.m_x0[v[1]] - obj.m_x0[v[0]], obj.m_x0[v[2]] - obj.m_x0[v[0]]) / 2).squared_length());
        }
        Vector spel_orient = model.sld.orientation;
        double width = 0.1;
        if (efd.free_edge_penalty_params.size() > 1) width = efd.free_edge_penalty_params[1];
        auto plane_force = [f = model.sld.orientation, p0 = model.sld.commissure_points[0], width](double* x, double* F, double& d){
            double dist = (Point(x[0], x[1], x[2]) - p0) * f - d;
            double w_f = (2*dist < width) ? (1 - 2*dist / width) : 0;
            for (int k = 0; k < 3; ++k) F[k] = w_f * f[k];
            return true;
        };
        auto plane_force_d = [f = model.sld.orientation, p0 = model.sld.commissure_points[0], width](double* x, double* J, double& d){
            std::fill(J, J+9, 0.0);
            double dist = (Point(x[0], x[1], x[2]) - p0) * f - d;
            if (2*dist < width){
                double w_j = -2 / width;
                for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    J[3*i + j] = w_j * f[i] * f[j];
            }
            return true;
        };
        std::vector<std::pair<V_ind, VertexLoad<double>::VertParam>> free_edge;
        for (auto v: obj.m_mesh.vertices()){
            auto& p0 = model.sld.commissure_points[0];
            auto& f = model.sld.orientation;
            if ((obj.m_boundary[v] & (4|1)) && !(obj.m_boundary[v] & 2))
                free_edge.emplace_back(v, VertexLoad<double>::VertParam{1, (obj.m_x[v] - p0)*f});
        }
        double P = operator""_mmHg(efd.pressure) / 1000;
        double W = P*S*efd.free_edge_penalty_params[0] / free_edge.size();
        for (auto& d: free_edge) d.second.weight *= W;
        Force Pre = VertexLoad<double>().insertVertices(free_edge).setForceField(plane_force, plane_force_d);

        model.m_w.addForce(Pre);
    }
    {
        double P = operator""_mmHg(efd.pressure) / 1000; //kPa
        Force Pr = SimplePressureLoad(P);
        pf_id = model.m_w.addForce(Pr, lid);
    }
    // WorldForce_ID ctt_id;
    // {
    //     ContactForce cf;
    //     cf.setContactType(ContactForce::VertFaceSelf | ContactForce::VertFace
    //                         | ContactForce::VertVertSelf | ContactForce::VertVert);
    //     double sigma = +0.5*0.5*2*(operator""_mmHg(efd.pressure) / 1000)*20*50;
    //     cf.setDistanceFunc([sigma](int contact_type, const ContactForce::Contact& c)->double{
    //         double thickness = (c.elems[0].obw->shape_thickness + c.elems[1].obw->shape_thickness)/2;
    //         return (thickness > c.d) ? (sigma * (thickness - c.d) / thickness)  : 0;
    //     }, 
    //     [sigma](int contact_type, const ContactForce::Contact& c)->double{
    //         double thickness = (c.elems[0].obw->shape_thickness + c.elems[1].obw->shape_thickness)/2;
    //         return (thickness > c.d) ? (-sigma/thickness) : 0;
    //     });
    //     ctt_id = model.m_w.addForce(std::move(cf));
    //     model.m_w.force(ctt_id).target<ContactForce>()->registerWorld(model.m_w);
    //     model.m_w.force(ctt_id).target<ContactForce>()->setShapeThickness(lid, 0.1);
    // }

    WorldForce_ID dmf_id;
    double damp_dt = 1;
    {
        auto damp_func = [&damp_dt, &_mass](Object3D* obj, V_ind v){
            return _mass[v] / (damp_dt * damp_dt);
        };
        WVertDampForce Df;
        auto dmpf = std::make_shared<WVertDampForce::DampFuncEvalDist>(2.0, 0.1, 50);
        dmpf->setDropSigmaCoef(0).setMaxSigmaCoef(1e4);
        Df.setDampFunc(std::move(dmpf));
        dmf_id = model.m_w.addForce(std::move(Df));
        model.m_w.force(dmf_id).target<WVertDampForce>()->registerWorld(model.m_w);
    }

    World3d::Timer sym_time;
    double abs_err_init = model.m_w.getWorldNorm();
    if (false) {
        double delta = rsd.delta;
        int maxits = rsd.maxits;
        double err = rsd.err;
        double diverge_lim = 1e15;
        int pfreq = rsd.pfreq;
        auto stopCond = World3d::StaticStopCondition(&err, &diverge_lim, &maxits, &pfreq, &sym_time, nullptr, true);
        StaticForceApplierP sfa(delta);
        model.m_w.setForceApplier(sfa);

        model.m_w.Simulation(stopCond);
        save_new_step(1);
    }

    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BridsonCollisionManager>();
    model.m_w.setCollider(std::move(colMan));
    BridsonCollisionManager* pcol = reinterpret_cast<BridsonCollisionManager*>(model.m_w.getCollider());
    DReal Ht = 1.0;
    CollisionThicknessCompressed<CollisionThicknessConst> cthickness{CollisionThicknessConst(Ht)};
    cthickness.setMaxCompression(4);
    pcol->set_thickness(lid, CollisionThickness(cthickness));
    pcol->set_Kspring(operator""_mmHg(efd.pressure) / 1000 * 0.5 * 0.5 / 0.05 * 50);
    pcol->set_Dt(0.1);
    // pcol->set_PenetrateDepthCoef(0.1);
    // pcol->set_RepulsSprCoef(0.5);
    pcol->setAsInitialNonCollidingState();
//    model.m_w.removeForce(lid, ef_id);
//    if (mtd.with_bending) model.m_w.removeForce(lid, bf_id);
//    set_force(3, mtd.thickness, mtd.E, mtd.nu, mtd.with_bending, 3, mtd.thickness_bend, mtd.E_bend, mtd.nu_bend);

    NSWorldWrapper nsww(model.m_w);
    NLProblem nlp = NSWorldWrapper::makeProblem(nsww);
    NonLinearSolverKinsol nlsp(nlp);
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
    nlsp.setInterIterationMonitorFunc([&nsww, &time = sym_time, pcol, &save_new_step, &read_step](const double * x){
        nsww.setX(x);
        // pcol->findCollisions();
        // pcol->solveCollisions();
        // save_new_step(read_step+1);
        std::cout << "time = " << time.elapsed() << std::endl;
        nsww.RenderIteration();
        return 0;
    });
    nlsp.SetFuncNormTol(1e-7 * abs_err_init);
    nlsp.SetNumMaxIters(50);
    nlsp.SetMaxNewtonStep(R / 10 * sqrt(nlp.m_dofs));
    bool slvFlag = false;
    nlsp.SetParameterInt("MaxSetupCalls", 1);
    nlsp.SetParameterInt("MaxSubSetupCalls", 1);
    nlsp.SetScaledStepTol(1e-5);
    nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
    nlsp.SetMaxBetaFails(80);

    pcol->set_MaxRelaxIts(0);
    pcol->set_RelaxEps(0.5);
    damp_dt = 0.05;
    model.m_w.force(dmf_id).target<WVertDampForce>()->resetCurrentState();
    auto solve_problem = [&](bool resolve_collisions = true){
        bool repeat = false;
        slvFlag = nlsp.Solve();
        bool make_col = resolve_collisions || (nlsp.GetReasonFlag() == KIN_LINESEARCH_NONCONV && nlsp.GetNumNolinSolvIters() <= 5);
        if (make_col){
            pcol->findCollisions();
            pcol->solveCollisions();
            if (!resolve_collisions && (nlsp.GetReasonFlag() == KIN_LINESEARCH_NONCONV && nlsp.GetNumNolinSolvIters() <= 5))
                repeat = true;
        }
        nlsp.SetInitialGuess(nsww.getCurrentX());
        nsww.RenderIteration();
        model.m_w.force(dmf_id).target<WVertDampForce>()->resetCurrentState();
        if (!make_col){
            pcol->findCollisions();
            pcol->solveCollisions();
        }
        return repeat;
    };
    for (int i = 0; i < 3; ++i){
        // if (i == 0){
        //     damp_dt = 0.01;
        //     auto dfc = std::make_shared<WVertDampForce::DampFuncConst>([&damp_dt](Object3D& obj, const ObjectID id, V_ind v){ return 1.0/(damp_dt*damp_dt); });
        //     model.m_w.force(dmf_id).target<WVertDampForce>()->setDampFunc(std::move(dfc));
        // }
        solve_problem();
        // slvFlag = nlsp.Solve();
        // pcol->findCollisions();
        // pcol->solveCollisions();
        // nlsp.SetInitialGuess(nsww.getCurrentX());
        // nsww.RenderIteration();
        if (i == 1){
            // damp_dt = 1.0;
            auto dfc = std::make_shared<WVertDampForce::DampFuncConst>([&damp_dt](Object3D& obj, const ObjectID id, V_ind v){ return 1.0/(damp_dt*damp_dt); });
            model.m_w.force(dmf_id).target<WVertDampForce>()->setDampFunc(std::move(dfc));
        }
        // model.m_w.force(dmf_id).target<WVertDampForce>()->resetCurrentState();
        //sleep(10);
    }
    damp_dt = 10;
    pcol->set_MaxRelaxIts(5);
    save_new_step(2);
//    model.m_w.Simulation(stopCond);
    nlsp.SetNumMaxIters(50);
//    double Pri = efd.pressure * 0.9;
//    for (int i = 0; i < 10; ++i)
    {
        {
            model.m_w.removeForce(lid, ef_id);
            if (mtd.with_bending) model.m_w.removeForce(lid, bf_id);
            set_force(mtd.MatType, mtd.thickness, mtd.E, mtd.nu, mtd.with_bending, mtd.MatType_bend, mtd.thickness_bend,
                      mtd.E_bend, mtd.nu_bend);
            model.m_w.removeForce(lid, pf_id);
            double Pri = 100;
            double P = operator ""_mmHg(Pri) / 1000; //kPa
            Force Pr = SimplePressureLoad(P);
            pf_id = model.m_w.addForce(Pr, lid);
            //slvFlag = nlsp.Solve();
            solve_problem();
            save_new_step(3);
//            model.m_w.Simulation(stopCond);
//            Pri *= 0.95;
        }

        {
            model.m_w.removeForce(lid, ef_id);
            if (mtd.with_bending) model.m_w.removeForce(lid, bf_id);
            set_force(mtd.MatType, mtd.thickness, mtd.E, mtd.nu, mtd.with_bending, mtd.MatType_bend, mtd.thickness_bend,
                      mtd.E_bend, mtd.nu_bend);
            model.m_w.removeForce(lid, pf_id);
            double Pri = 60;
            double P = operator ""_mmHg(Pri) / 1000; //kPa
            Force Pr = SimplePressureLoad(P);
            pf_id = model.m_w.addForce(Pr, lid);
            //slvFlag = nlsp.Solve();
            solve_problem();
            save_new_step(4);
//            model.m_w.Simulation(stopCond);
//            Pri *= 0.95;
        }
    }

    auto aorta_sdf = std::make_shared<CilDist>(Vector{0, 0, 0}, Vector{0, 0, 1}, R+0.5);
    Force aorta_f = SDFForce(aorta_sdf, 1 * operator ""_mmHg(efd.pressure) / 1000, 0.5);
    auto art_id = model.m_w.addForce(std::move(aorta_f), lid);
    if (!is_bezier){
        Force_ID lc_id = -1, rc_id = -1;
        double d_lr = -4.8;
        nlsp.SetNumMaxIters(10);
        for (int i = 0, cnt = 20, cnt_i = 25; i < cnt_i; ++i){
            std::cout << "\nSpurnPlane loop it = " << i << std::endl;
            if (lc_id) model.m_w.removeForce(lid, lc_id);
            if (rc_id) model.m_w.removeForce(lid, rc_id);
            Force left_constraint = SpurnPlane(1 * (operator""_mmHg(efd.pressure) / 1000) / cnt * i, 0.02 * R, SpurnPlane::Plane_3(Point(d_lr + 0.25*std::min(i, 20), 0, 0),
                                                                                Vector(-sin(-M_PI / 3), cos(-M_PI / 3),0.0)));
            Force right_constraint = SpurnPlane(1 * (operator""_mmHg(efd.pressure) / 1000) / cnt * i, 0.02 * R, SpurnPlane::Plane_3(Point(d_lr + 0.25*std::min(i, 20), 0, 0),
                                                                                    -Vector(-sin(M_PI / 3), cos(M_PI / 3), 0.0)));
            lc_id = model.m_w.addForce(left_constraint, lid);
            rc_id = model.m_w.addForce(right_constraint, lid);
            //slvFlag = nlsp.Solve();
            solve_problem();
            if (i % 5 == 4) save_new_step(5 + i / 5);
        }
    } else {
        damp_dt = 10;
        Force_ID sd_id = -1;
        std::array<Point, 3> A = m_commissure_points;
        std::array<DReal, 3> l = leaf_lengths;
        CilindricContactSurface ccs(A, l);
        auto curve = ccs.getLeafContactSdf(0);
        std::cout << curve->m_main_lines[0] << "\n" << curve->m_main_lines[1] << "\n";
        Force sdf_f = SDFForce(curve, 1 * operator ""_mmHg(efd.pressure) / 1000, 1);
        sd_id = model.m_w.addForce(std::move(sdf_f), lid);
        nlsp.SetNumMaxIters(20);
        // model.m_w.obj(lid).remove_force(bf_id);
        for (int st = 61, i = st, cnt = 60; i <= cnt; ++i){
            std::cout << "\nBezier loop it = " << i << std::endl;
            model.m_w.obj(lid).m_forces[sd_id].target<SDFForce>()->field<CilindricContactSurface::LeafContactSDF>()->set_lambda(i*1.0/cnt);
            solve_problem();
            if ((i-st) % 5 == 0 && i > st) save_new_step(5 + (i-st) / 5);
        }
        auto lc_sdf = model.m_w.obj(lid).m_forces[sd_id].target<SDFForce>()->field<CilindricContactSurface::LeafContactSDF>();
        lc_sdf->set_lambda(1);
        struct ContSurfIter{
            int maxit = 50;
            DReal cP = 1.0;
            DReal Ht = 0.5;
        };
        std::vector<ContSurfIter> csi{
                                    // ContSurfIter{10, 1, 1.0}, 
                                    // ContSurfIter{10, 0.6, 1.0}, 
                                    // ContSurfIter{10, 0.6, 0.5}, 
                                    // ContSurfIter{10, 0.6, 0.3},
                                    // ContSurfIter{10, 0.6, 0.2},
                                    // ContSurfIter{10, 0.6, 0.15}, 
                                    // ContSurfIter{10, 0.6, 0.1},
                                    // ContSurfIter{10, 0.6, 0.1},
                                    // ContSurfIter{10, 0.6, 0.1},
                                    // ContSurfIter{10, 0.6, 0.1},
                                    // ContSurfIter{10, 0.6, 0.1},
                                    // ContSurfIter{10, 0.6, 0.1},
                                    // ContSurfIter{10, 0.6, 0.1},
                                    
                                    ContSurfIter{100, 1, 1.0}, 
                                    ContSurfIter{50, 0.6, 1.0}, 
                                    ContSurfIter{50, 0.6, 0.5}, 
                                    ContSurfIter{50, 0.6, 0.3},
                                    ContSurfIter{50, 0.6, 0.2},
                                    ContSurfIter{50, 0.6, 0.15}, 
                                    ContSurfIter{50, 0.6, 0.1},
                                    ContSurfIter{50, 0.6, 0.1},
                                    ContSurfIter{50, 0.6, 0.1},

                                    // ContSurfIter{100, 0.5, 0.5}, 
                                    // ContSurfIter{100, 0.04, 0.5},
                                    // ContSurfIter{50, 0.08, 0.25},
                                    // //ContSurfIter{50, 0.1, 0.2},
                                    // //ContSurfIter{50, 0.16, 0.2},
                                    // //ContSurfIter{50, 0.2, 0.2},
                                    // /*ContSurfIter{50, 0.3, 0.2},*/ContSurfIter{50, 0.6, 0.5},
                                    // /*ContSurfIter{50, 0.4, 0.2},*/ContSurfIter{50, 0.6, 0.25}, 
                                    // /*ContSurfIter{50, 0.6, 0.2},*/ContSurfIter{50, 0.6, 0.17}, 
                                    // ContSurfIter{50, 0.6, 0.1},  
                                    // ContSurfIter{50, 0.6, 0.1}
                                    };
        for (int k = 0; k < csi.size(); ++k){
            auto shift = std::max(csi[k].Ht - dphi*R/2, 0.0);
            lc_sdf->make_shift(Vector{-shift, 0, 0});
            nlsp.SetNumMaxIters(csi[k].maxit);
            if (sd_id >= 0) model.m_w.obj(lid).remove_force(sd_id);
            sdf_f = SDFForce(curve, csi[k].cP * operator ""_mmHg(efd.pressure) / 1000, csi[k].Ht);
            sd_id = model.m_w.addForce(std::move(sdf_f), lid);
            bool repeat = solve_problem(k == csi.size()-1);
            // bool repeat = solve_problem();
            lc_sdf->make_shift(Vector{+shift, 0, 0});
            if (repeat) { std::cout << "Try repeat" << std::endl; --k; continue; }
            save_new_step(30 + k);
        }
        // nlsp.SetNumMaxIters(50);
        // for (int k = 1; k <= 10; ++k){
        //     // model.m_w.obj(lid).remove_force(sd_id);
        //     // Force sdf_f = SDFForce(curve, (3+k) * operator ""_mmHg(efd.pressure) / 1000, 0.02*R);
        //     // sd_id = model.m_w.addForce(std::move(sdf_f), lid);
        //     solve_problem();
        //     save_new_step(30 + k - 1);
        // }
        // model.m_w.obj(lid).remove_force(sd_id);
        // sdf_f = SDFForce(curve, 1 * operator ""_mmHg(efd.pressure) / 1000, 0.5);
        // sd_id = model.m_w.addForce(std::move(sdf_f), lid);
        // for (int k = 1; k <= 10; ++k){
        //     // model.m_w.obj(lid).remove_force(sd_id);
        //     // Force sdf_f = SDFForce(curve, (3+k) * operator ""_mmHg(efd.pressure) / 1000, 0.02*R);
        //     // sd_id = model.m_w.addForce(std::move(sdf_f), lid);
        //     solve_problem();
        //     save_new_step(40 + k - 1);
        // }
        // // nlsp.SetNumMaxIters(50);
        // // solve_problem();
    }

    if (with_gui_window){
        front_end->join();
        delete front_end;
    }

    return 0;
}

// struct HeightLeafParams{
//     double Hf;
//     double Hi;
//     double Rs;
// };

// int create_quasi_optimal_valve(int argc, char* argv[], 
//         std::array<HeightLeafParams, 3> templs, double R, std::array<double, 3> suture_relations){
//     std::vector<std::string> args(argv, argv+argc);
//     args.push_back("bnd.R=12");

// }

int test_CCD(){
    std::string pnts_str = "E0: v0: 0x1.5a14c6979fd7fp+1 -0x1.1656338d3c01p+4 0x1.33373fd660328p+1 -> -0x1.b3940f88a8cp-6 0x1.2f09a2c01fp-7 0x1.7f70519d4a2p-5 v1: 0x1.58998076ade86p+1 -0x1.0b2b1e12a909dp+4 0x1.3e02b1a161179p+1 -> -0x1.9a966f5bd9p-7 0x1.23ba555552p-9 0x1.8f477ff935ccp-5 E1: v0: 0x1.5b8c4c6e7fba2p+1 -0x1.21650212c9a4dp+4 0x1.288722d9c8c78p+1 -> 0x1.1714c6629f32p-2 0x1.a4bdcef51bep-2 0x1.67c8dfca3ca8p-4 v1: 0x1.670788ba119a2p+1 -0x1.252f87a5c5a32p+4 0x1.a9bee2e2e4a05p+0 -> 0x1.4c77768028dc8p-1 0x1.64b69196af93p+0 0x1.ddea189bf48p-7";
    int type = 0;
    if (pnts_str[0] == 'E') type = 1;
    const int Npnt = 3*4*2;
    DReal pd[Npnt];
    {
        std::regex rex(R"(([+-]?(\d+\.\d*))|([-]?0x\d(\.[0-9a-f]*)?p[+-]\d+))");
        auto nums_begin = 
            std::sregex_iterator(pnts_str.begin(), pnts_str.end(), rex);
        auto nums_end = std::sregex_iterator();

        if (Npnt != std::distance(nums_begin, nums_end)) throw std::runtime_error("Recognised wrong number of values");
        int n = 0;
        for (std::sregex_iterator i = nums_begin; i != nums_end; ++i) {
            std::smatch match = *i;
            std::string match_str = match.str();
            pd[n++] = stod(match_str);
        }
    }
    std::array<Point, Npnt/6> p;
    std::array<Vector, Npnt/6> dp;
    for (int i = 0; i < Npnt/3; ++i){
        if (i%2 == 0) p[i/2] = Point(pd[3*i + 0], pd[3*i + 1], pd[3*i + 2]);
        else dp[i/2] = Vector(pd[3*i + 0], pd[3*i + 1], pd[3*i + 2]);
    }
    DReal t = -1;
    bool res = false;
    if (type == 0){
        Point o = p[0], z = Point(0, 0, 0);
        res = GeomProjs::PerformVFCCD(Vector(CGAL::NULL_VECTOR), p[1] - o, p[2] - o, p[3] - o, dp[0], p[1] - o + dp[1], p[2] - o + dp[2], p[3] - o + dp[3], &t);
        if (res){
            auto    ac = (p[0] - z) + t*dp[0], 
                    bc = (p[1] - z) + t*dp[1],
                    cc = (p[2] - z) + t*dp[2],
                    pc = (p[3] - z) + t*dp[3];
            auto sqdt = (GeomProjs::BaryEval(ac, bc, cc, GeomProjs::BaryCoord(ac, bc, cc, pc)) - pc).squared_length();
            std::cout << "sqdt = " << sqdt << std::endl;
        }
    } else {
        Point o = p[0], z = Point(0, 0, 0);
        res = World3d::GeomProjs::PerformEECCD(Vector(CGAL::NULL_VECTOR), p[1] - o, p[2] - o, p[3] - o, dp[0], p[1] - o + dp[1], p[2] - o + dp[2], p[3] - o + dp[3], &t);
        if (res){
            Vector prj0, prj1; double w[2]; double sqdt; Point o = CGAL::ORIGIN;
            auto    ac = (p[0] - o) + t*dp[0], 
                    bc = (p[1] - o) + t*dp[1],
                    cc = (p[2] - o) + t*dp[2],
                    dc = (p[3] - o) + t*dp[3];
            GeomProjs::SegmentsProjects(ac, bc, cc, dc, prj0, prj1, w[0], w[1], sqdt);
            std::cout << "sqdt = " << sqdt << std::endl;
        }
    }
    std::cout << "t = " << t << " res = " << res << std::endl;
    return 0;
}

int test_DCD(){
    using namespace GeomProjs;
    Vector V[] = {Vector{0, 0, 0  }, {1, 0, 0  }, {0, 1, 0  }, {1, 1, 0  },
                        {0, 0, 1  }, {1, 0, 1  }, {0, 1, 1  }, {1, 1, 1  },
                        {0, 0, 0.5}, {1, 0, 0.5}, {0, 1, 0.5}, {1, 1, 0.5},
                        {0.25,0.25,0}, {0.25,0.25,0.5}, {0.25, 0.25,1}, {0.25, 0.25,2},
                        {3, 4, 5  }, {9, 11, 7 }, {13, 2, 4 }, {-1,-1, -1}, 
                        {0.25,0.25,-0.01}};
    Vector prj[2];
    DReal wf[3] = {0}, wq[3] = {0}, we[2] = {0}, sqd = -1;
    auto print_state = [&](int* inds, int n, int stat = 1 | 2 | 4){
        std::cout << "V[" << n << "] = {";
        for (int i = 0; i < n-1; ++i)
            std::cout << V[inds[i]] << ", ";
        if (n > 0) std::cout << V[inds[n-1]];
        std::cout << "}\n";
        std::cout << "dist = " << sqrt(sqd) << " ( sqrt(" << sqd << "))" << "\n";
        std::cout   << "prj[0] = " << prj[0] << "\n" 
                    << "prj[1] = " << prj[1] << "\n";
        if (stat & (1 << 0)) 
            std::cout << "wf = {" << wf[0] << ", " << wf[1] << ", " << wf[2] << "}\n";
        if (stat & (1 << 1)) 
            std::cout << "wq = {" << wq[0] << ", " << wq[1] << ", " << wq[2] << "}\n";
        if (stat & (1 << 2))
            std::cout << "we = {" << we[0] << ", " << we[1] << "}\n";
    };
    auto perform_edge_test = [&](std::array<int, 5> inds){
        GeomProjs::FaceEdgeProjects(V[inds[0]], V[inds[1]], V[inds[2]], V[inds[3]], V[inds[4]], prj[0], prj[1], wf, we, sqd);
        print_state(inds.data(), 5, 1|4);
    };
    auto perform_face_test = [&](std::array<int, 6> inds){
        GeomProjs::FaceFaceProjects(V[inds[0]], V[inds[1]], V[inds[2]], V[inds[3]], V[inds[4]], V[inds[5]], prj[0], prj[1], wf, wq, sqd);
        print_state(inds.data(), 6, 1|2);
    };
    perform_edge_test({0, 1, 2, 4, 8});
    perform_edge_test({8, 9, 10, 12, 14});
    perform_edge_test({8, 9, 10, 14, 19});
    perform_edge_test({0, 1, 2, 2, 1});
    perform_edge_test({0, 1, 2, 2, 2});
    perform_edge_test({12, 13, 14, 2, 2});
    std::cout << "###################################\n";
    perform_face_test({0, 1, 2, 4, 8, 18});
    perform_face_test({8, 9, 10, 12, 14, 18});
    perform_face_test({8, 9, 10, 14, 19, 18});
    perform_face_test({0, 1, 2, 2, 1, 18});
    perform_face_test({0, 1, 2, 2, 1, 0});
    perform_face_test({0, 1, 2, 12, 16, 17});
    perform_face_test({0, 1, 2, 20, 16, 17});
    perform_face_test({0, 1, 2, 12, 13, 14});
    auto res = GeomProjs::Proximity::Query(GeomProjs::Proximity::Primitive(1, {V[0], V[1], V[2]}), GeomProjs::Proximity::Primitive(2, {V[12], V[13], V[14]}));
    std::cout << "res.sqd = " << res.sqd << std::endl; 
    std::cout << res << std::endl;

    return 0;
}

int test_ContactSurface(){
    double alpha = 1.4;
    std::array<double, 3> a{6, 9, 13};
    std::array<double, 3> l;
    for (int i = 0; i < 3; ++i) l[i] = alpha*a[i];
    double cos_x = (a[0]*a[0] + a[2]*a[2] - a[1]*a[1])/(2*a[0]*a[2]);
    std::array<Vector, 3> A{Vector{a[2]*cos_x, a[2]*sqrt(1 - cos_x*cos_x), 0}, {0, 0, 0}, {a[0], 0, 0}};
    CilindricContactSurface ccs(A, l);
    auto curve = ccs.getLeafContactSdf(0);
    std::cout << ccs.m_C << std::endl;
    std::cout << curve->m_main_lines[0] << "\n" << curve->m_main_lines[1] << "\n";
    std::cout << ArcLength(curve->m_main_lines[0]) << " " << ArcLength(curve->m_main_lines[1]) << "\n";
    std::cout << curve->m_main_derivatives[0] << "\n" << curve->m_main_derivatives[1] << "\n";
    std::cout << curve->m_current_lines[0] << "\n" << curve->m_current_lines[1] << "\n";
    DReal lmb[] = {0, 0.001, 0.1, 0.5, 0.999, 1.0};
    for (int i = 0; i < sizeof(lmb)/sizeof(DReal); ++i){
        curve->set_lambda(lmb[i]);
        std::cout << "lambda = " << lmb[i] << std::endl;
        std::cout << (*curve)(A[1]).to_string() << "\n"
                << (*curve)(A[0]).to_string() << "\n"  
                << (*curve)(ccs.m_C).to_string() << "\n"
                << (*curve)(ccs.m_C - Vector(0, 0.3, 0)).to_string() << "\n"
                << (*curve)(ccs.m_C + Vector(0, 0.3, 0)).to_string() << "\n"
                << (*curve)(ccs.m_C + Vector(1, 0.3, 0)).to_string() << "\n"
                << (*curve)(ccs.m_C - Vector(1, 0.3, 0)).to_string() << "\n" << std::endl;
    }
    auto line = [C = ccs.m_C, x = Vector(0, -10, 10), &curve](double t){
        return (*curve)(C + x*t);
    };
    std::cout << line(0).to_string() << "\n " << line(0.5).to_string() << "\n " << line(1.0).to_string() << "\n" <<  line(2.0).to_string() << "\n";

    return 0;
}

int main(int argc, char* argv[]){
    //return test_ContactSurface();
    //return test_DCD();
    //return test_CCD();
    return make_sewing(argc, argv);
}
