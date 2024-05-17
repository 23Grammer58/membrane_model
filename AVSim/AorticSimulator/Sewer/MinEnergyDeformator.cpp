#include "MinEnergyDeformator.h"

class LengthFunctor: public MinEnergyDeformator::EnergyFunctor{
public:
    double weight;
    LengthFunctor (double weigth = 0) : weight{weigth} {}
    virtual double operator()(const EnergyDeformatorMesh& mesh, int edge_id, double df[6]) override {
        double alpha = weight;
        auto l = mesh.l(edge_id), l_0 = mesh.l0(edge_id);
        double eps = (l - l_0) / l_0;
        double f = alpha * eps * eps;
        if (df) {
            auto ve = mesh.edge(edge_id);
            auto x = std::array<Vector,2>{mesh.x(ve[0]), mesh.x(ve[1])};
            double dcoef = 2 * alpha * eps / l / l_0;
            for (int i = 0; i < 3; ++i){
                df[i] = (x[0][i] - x[1][i]) * dcoef;
                df[i+3] = -df[i];
            }
        }
        return f;
    }
};

class DigedralFunctor: public MinEnergyDeformator::EnergyFunctor {
public:
    double weight;
    double scale;
    bool convexity;

    DigedralFunctor (double weight = 0, double scale = 1, bool convexity = true) :
    weight{weight}, scale{scale}, convexity{convexity} {}
    virtual double operator()(const EnergyDeformatorMesh& mesh, int edge_id, double df[12]) override {
        auto [ve, is_digedral] = mesh.digedral(edge_id);
        Vector p[4] = {mesh.x(ve[2]), mesh.x(ve[3]), mesh.x(ve[0]), mesh.x(ve[1])};
        double x[12] = {}, k;
        if (df) {
            double phi = deriv_f3(p, x, 1.0, scale, &k);
            double weightphik = weight * phi;
            double _edge_coef = -2 * weightphik;
            for (int i = 0; i < 12; ++i) {
                x[i] *= _edge_coef;
                df[i] += x[i];
            }
            return weightphik * phi * k;
        }
        else{
            return weight * get_phi3(p, x, 1.0, scale);
        }
    }

private:
    inline double get_phi(const Vector p[3], Vector res[3], double sc[2]){
        double a = p[0].stableNorm();
        Vector axb = p[0].cross(p[1]), cxa = p[2].cross(p[0]);
        double abc = axb.dot(p[2]), axb_axc = -axb.dot(cxa);
        double phi = abc * a;
        double psi = axb_axc;
        double f = phi / psi;
        sc[0] = phi, sc[1] = -psi;

        return f;
    }
   inline double get_phi2(const Vector p[3], Vector res[3], double beta, double scal){
        double sc[2];
        double x = get_phi(p, res, sc);
        double f = atan2(sc[0], sc[1]) / M_PI;
        double k = ((convexity) ? (f > 0) : (f < 0)) ? 0.1 : (scal * 0.1);

        return k * f * f;
    }

    inline double get_phi3(const Vector p[4], double x[12], double beta, double scal){
        int mask[3] = {3, 0, 1};
        Vector pp[3];
        for (int i = 0; i < 3; ++i)
            pp[i] = p[mask[i]] - p[2];
        double f = get_phi2(pp, NULL, beta, scal);
        return f;
    }

    double deriv_f1(const Vector p[3], Vector res[3], double sc[2]){
        //f = (a, b, c) * |a| / ([a, b], [a, c])
        //p[0] <-> a
        //p[1] <-> b
        //p[2] <-> c
        double a2 = p[0].squaredNorm();
        double a = p[0].stableNorm();
        double ab = p[0].dot(p[1]), ac = p[0].dot(p[2]), bc = p[1].dot(p[2]);
        Vector axb = p[0].cross(p[1]), cxa = p[2].cross(p[0]), bxc = p[1].cross(p[2]);
        double abc = axb.dot(p[2]), axb_axc = -axb.dot(cxa);
        //std::array<point_t&, 3> abc_p = {bxc, cxa, axb};
        double phi = abc * a;
        //point_t abc_adna = SCAL(abc/a, p[0]);
        //for (int i = 1; i < 3; ++i)
        //    SCAL_S(a, &abc_p);
        //point_t phi_da = SCAL_SUM(a, bxc, abc/a, p[0]);
        std::array<Vector, 3> phi_p = {a*bxc + abc/a *p[0], a*cxa, a*axb};
        double psi = axb_axc;
        Vector psi_da = 2 * bc * p[0] - ac * p[1] - ab * p[2];
        Vector psi_p[3] = {psi_da, a2 * p[2] - ac * p[0], a2 * p[1] - ab * p[0]};
        double f = phi / psi;
        sc[0] = phi, sc[1] = -psi;
        double d1_psi = 1.0/psi, phid2psi = -phi / psi / psi;
        for (int i = 0; i < 3; ++i){
            res[i] = d1_psi * phi_p[i] + phid2psi * psi_p[i];
        }
        return f;
    }

    double deriv_f2(const Vector p[3], Vector res[3], double beta, double scal, double* kk){
        double sc[2];
        double x = deriv_f1(p, res, sc);
        double f = atan2(sc[0], sc[1]) / M_PI;
        double coef_arctg = 1.0 / (1.0 + x * x);
        double k = ((convexity) ? (f > 0) : (f < 0)) ? 0.1 : (scal * 0.1);
        double coef = coef_arctg * k * beta / M_PI;
        for (int i = 0; i < 3; ++i)
            res[i] *= coef;
        *kk = k;
        return f;
    }

    double deriv_f3(const Vector p[4], double x[12], double beta, double scal, double* kk){
        int mask[3] = {3, 0, 1};
        Vector pp[3];
        for (int i = 0; i < 3; ++i)
            pp[i] = p[mask[i]] - p[2];
        Vector grad[3];
        double f = deriv_f2(pp, grad, beta, scal, kk);

        if (x != nullptr)
            for (int i = 0; i < 3; ++i) {
                x[3 * 0 + i] = grad[1][i];
                x[3 * 1 + i] = grad[2][i];
                x[3 * 2 + i] = -(grad[0][i] + grad[1][i] + grad[2][i]);
                x[3 * 3 + i] = grad[0][i];
            }

        return f;
    }

};

class PlaneConstrFunctor: public MinEnergyDeformator::EnergyFunctor {
public:
    Vector normal;
    double m_b;   // (n, x) = b
    double weight;

    PlaneConstrFunctor(double weight, const Vector& normal, double b) :
    weight{weight}, normal{normal}, m_b{b} {}
    virtual double operator()(const EnergyDeformatorMesh& mesh, int node_id, double* df) override {
        double pn = normal.dot(mesh.x(node_id));
        double pnmb = pn - m_b;
        double mpnmb2 = pnmb * 2 * weight;
        bool intersect = (pnmb > 0);
        if (df && intersect){
            for (int i = 0; i < 3; ++i)
                df[i] += mpnmb2 * normal[i];
        }
        return (intersect) ? weight * pnmb * pnmb : 0;
    }
};

class ForceFunctor: public MinEnergyDeformator::EnergyFunctor {
public:
    Vector force;
    double weight;
private:
    Vector m_f;
public:

    ForceFunctor(double weight, const Vector& force) :
            weight{weight}, force{force}, m_f{-weight * force} {}
    virtual double operator()(const EnergyDeformatorMesh& mesh, int node_id, double* df) override {
        if (df)
            for (int i = 0; i < 3; ++i)
                df[i] += m_f[i];
        double f = m_f.dot(mesh.x(node_id));
        return f;
    }
};

void set_default_length_constr(MinEnergyDeformator& m, double weight){
    if (weight == 0) return;
    auto en = std::make_shared<LengthFunctor>(weight);
    m.addEdgeEnergyComponent(en);
//        m.addEdgeEnergyComponent(new LengthFunctor(weight));
}

void set_default_digedral_angle_constr(MinEnergyDeformator& m, double weight, double scale, bool convexity){
    if (weight == 0) return;
    auto en = std::make_shared<DigedralFunctor>(weight, scale, convexity);
    m.addDigedralEnergyComponent(en);
//        m.addDigedralEnergyComponent(new DigedralFunctor(weight, scale, convexity));
}

void set_plane_constr(MinEnergyDeformator& m, double weight, const Eigen::Vector3d& normal, double b){
    if (weight == 0) return;
    auto en = std::make_shared<PlaneConstrFunctor>(weight, normal, b);
    m.addNodeEnergyComponent(en);
//        m.addNodeEnergyComponent(new PlaneConstrFunctor(weight, normal, b));
}

void set_isotrop_force(MinEnergyDeformator& m, double weight, const Eigen::Vector3d& force){
    if (weight == 0) return;
    auto en = std::make_shared<ForceFunctor>(weight, force);
    m.addNodeEnergyComponent(en);
}
