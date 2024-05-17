#include "ContactSurface.h"
#include "find_polynomial_roots_jenkins_traub.h"
#include <gsl/gsl_poly.h>

using namespace World3d;
void CilindricContactSurface::setup(){
    auto pn = computeNormalAndArea(m_A);
    std::array<DReal, 3> c;
    {
        DReal p = (m_l[0] + m_l[1] + m_l[2])/2;
        for (int i = 0; i < 3; ++i) {
            c[i] = p - m_l[i];
            if (c[i] <= 0) throw std::runtime_error("Leaflets are incompatible sizes");
        }
    }
    
    auto pd = computeDilationPoint(m_A, c);
    m_n = pn.first;
    m_C = pd.C;
    m_B = compute_B_points(m_A, m_C);
    for (int i = 0; i < 3; ++i){
        auto pd = compute_line(m_A[i], m_B[i], m_C, c[i]);
        m_lines[i] = pd.first; m_sigma[i] = pd.second;
    }
}
void CilindricContactSurface::setPointsAndWeights(const std::array<Vector, 3>& commissure_points, const std::array<DReal, 3>& leaf_width){
    m_A = commissure_points, m_l = leaf_width;
    setup();
}
std::shared_ptr<CilindricContactSurface::LeafPlaneContactSDF> CilindricContactSurface::getLeafPlaneContactSdf(int ileaf){
    LeafPlaneContactSDF l;
    l.m_n = m_n;
    l.m_lambda = 0.0;
    l.m_Cs = m_A[ileaf];
    l.m_A = m_A[(ileaf+1)%3], l.m_B = m_A[(ileaf+2)%3];
    l.m_C = m_C;
    return std::make_shared<LeafPlaneContactSDF>(std::move(l));
}

std::shared_ptr<CilindricContactSurface::LeafContactSDF> CilindricContactSurface::getLeafContactSdf(int ileaf){
    LeafContactSDF l;
    l.m_n = m_n;
    l.m_sigma[0] = m_sigma[(ileaf+1)%3], l.m_sigma[1] = m_sigma[(ileaf+2)%3];
    l.m_main_lines[0] = m_lines[(ileaf+1)%3], l.m_main_lines[1] = m_lines[(ileaf+2)%3];
    l.m_main_derivatives[0] = l.m_main_lines[0].derivative(), l.m_main_derivatives[1] = l.m_main_lines[1].derivative();
    l.opposit_point = m_A[ileaf]; 
    l.set_lambda(0.0);
    return std::make_shared<LeafContactSDF>(std::move(l));
}
std::pair<Vector,DReal> CilindricContactSurface::computeNormalAndArea(const std::array<Vector, 3>& A){
    Vector n = CGAL::cross_product(A[1]-A[0], A[2]-A[0]);
    DReal S = sqrt(n.squared_length()); 
    n /= S;
    S /= 2;
    return {n,S};
}
CilindricContactSurface::DilationPointRes CilindricContactSurface::computeDilationPointFromLengths(const std::array<Point, 3>& A, const std::array<DReal, 3>& l){
    std::array<DReal, 3> c;
    DReal p = (l[0] + l[1] + l[2])/2;
    for (int i = 0; i < 3; ++i) {
        c[i] = p - l[i];
        if (c[i] <= 0) throw std::runtime_error("Leaflets are incompatible sizes");
    }
    return computeDilationPoint(std::array<Vector, 3>{A[0] - CGAL::ORIGIN, A[1] - CGAL::ORIGIN, A[2] - CGAL::ORIGIN}, c);
}
CilindricContactSurface::DilationPointRes CilindricContactSurface::computeDilationPoint(const std::array<Vector, 3>& A, const std::array<DReal, 3>& c){
    for (int i = 0; i < 3; ++i){
        Vector K = (c[(i+2)%3]*A[(i+1)%3] + c[(i+1)%3]*A[(i+2)%3])/(c[(i+1)%3] + c[(i+2)%3]);
        DReal v_2 = (A[(i+1)%3] - A[(i+2)%3]).squared_length()/((c[(i+1)%3]+c[(i+2)%3])*(c[(i+1)%3]+c[(i+2)%3])), c_2 = c[i]*c[i];
        DReal AK_2 = (K - A[i]).squared_length();
        if (AK_2 <= v_2 * c_2) return DilationPointRes{K, sqrt(v_2), true};
    }
    DReal c_min = *min_element(c.begin(), c.end()), c_max = *max_element(c.begin(), c.end());
    if (c_max - c_min < 64 * (abs(c_max) + abs(c_min)) * std::numeric_limits<DReal>::epsilon()){
        std::array<DReal, 3> a;
        for (int i = 0; i < 3; ++i) a[i] = (A[(i+1)%3] - A[(i+2)%3]).squared_length();
        std::array<DReal, 3> l;
        DReal s = 0;
        for (int i = 0; i < 3; ++i){
            l[i] = a[i] * (a[(i+1)%3] + a[(i+2)%3] - a[i]);
            s += l[i];
        }
        for (int i = 0; i < 3; ++i) l[i] /= s;
        Vector C = l[0] * A[0] + l[1] * A[1] + l[2] * A[2];
        return DilationPointRes{C, sqrt((C - A[0]).squared_length())/c[0], true};
    }
    std::array<DReal, 3> a;
    for (int i = 0; i < 3; ++i) a[i] = sqrt((A[(i+1)%3] - A[(i+2)%3]).squared_length());
    DReal cos_x = (a[1]*a[1] + a[2]*a[2] - a[0]*a[0])/(2*a[1]*a[2]),
            cos_y = (a[0]*a[0] + a[2]*a[2] - a[1]*a[1])/(2*a[0]*a[2]);
    DReal sin_x = sqrt(1 - cos_x*cos_x);
    DReal sin_y = a[1]/a[0]*sin_x;
    DReal cos_x2 = sqrt((1+cos_x)/2), sin_x2 = sqrt((1-cos_x)/2);
    DReal r0 = (c[0]-c[2])*(c[0]+c[2])/a[1], r1 = (c[0]-c[1])*(c[0]+c[1])/a[2];
    DReal bx0 = (r0 + r1)/(4*cos_x2), bx1 = (a[1] + a[2])/(4*cos_x2),
            by0 = (r0 - r1)/(4*sin_x2), by1 = (a[1] - a[2])/(4*sin_x2);
    DReal ae = bx0 * bx0 + by0 * by0, 
            be = 2 * ( bx0 * bx1 + by0 * by1) - c[0]*c[0],
            ce = bx1*bx1 + by1*by1;
    DReal k = be / (2*ae), q = ce / ae;
    DReal D = sqrt(std::max(k*k - q, 0.0));
    std::array<DReal, 2> v_2;
    v_2[0] = -k + D, v_2[1] = -k - D;
    int nroot = 2;
    if (v_2[1] < 0) nroot = 1;
    if (v_2[0] <= 0) throw std::runtime_error("Equation doesn't have numerical solutions");
    std::array<DReal, 2> err = {0};
    for (int i = 0; i < nroot; ++i){
        DReal x = r1/2*v_2[i] + a[2]/2;
        DReal y = (r0/2*v_2[i] + a[1]/2 - cos_x*x)/sin_x;
        err[i] = std::max(0.0, -y) + std::max(0.0, -x*sin_x + y*cos_x) + std::max(0.0, y*cos_y - (a[2]-x)*sin_y);
    }
    DReal v = v_2[0];
    if (nroot == 2 && err[1] < err[0]) v = v_2[1];
    std::array<DReal, 3> l;
    v = sqrt(v);
    DReal s = 0;
    for (int i = 0; i < 3; ++i){
        DReal p = (v*c[(i+1)%3] + v*c[(i+2)%3] + a[i])/2;
        l[i] = sqrt(p*(p-a[i])*(p-v*c[(i+1)%3])*(p-v*c[(i+2)%3]));
        s += l[i];
    }
    for (int i = 0; i < 3; ++i) l[i] /= s;
    Vector C = l[0] * A[0] + l[1] * A[1] + l[2] * A[2];
    return DilationPointRes{C, v, true};
}
std::array<Vector, 3> CilindricContactSurface::compute_B_points(const std::array<Vector, 3>& A, const Vector& C){
    std::array<Vector, 3> B;
    Vector n = CGAL::cross_product(A[1]-A[0], A[2]-A[0]);
    DReal S = sqrt(n.squared_length()); 
    n /= S;
    S /= 2;
    DReal R = sqrt((A[1]-A[0]).squared_length() * (A[2]-A[0]).squared_length() * (A[2]-A[1]).squared_length())/(4*S);
    for (int i = 0; i < 3; ++i){
        if ((A[i] - C) * (A[(i+1)%3]-C) > 0){
            B[i] = C + (A[(i+1)%3]-C) * std::min(R / sqrt((C-A[(i+1)%3]).squared_length()), 1.0); 
        } else{
            B[i] = CGAL::cross_product(n, A[i] - C);
            B[i] /= sqrt(B[i].squared_length());
            B[i] = C + B[i]*std::min(R, sqrt((C-A[(i+1)%3]).squared_length()));
        }
    }
    return B;
}
std::pair<BezierCurve<3, Vector>, DReal> CilindricContactSurface::compute_line(const Vector& A, const Vector& B, const Vector& C, DReal c, DReal rel_err){
    if ((A - C).squared_length() >= c*c) return {BezierCurve<3, Vector>(std::array<Vector, 4>{C, C, C, A}), 0.0};
    const DReal D = 0.5;
    auto len_func = [&](DReal sigma)->DReal{
        std::array<Vector, 4> Y = {C, C + sigma*(B-C), C + sigma/(D + sigma)*(A-C), A};
        return ArcLength(BezierCurve<3, Vector>(Y), 0.0, 1.0, rel_err/8);
    };
    if ((B-C).squared_length() == 0) throw std::runtime_error("Wrong B point");
    DReal m = 0, M = 1;
    DReal fm = sqrt((A - C).squared_length()), fM = len_func(M);
    while(fM < c){
        m = M; fm = fM;
        M *= 2; fM = len_func(M);
    }
    while (fM - fm > rel_err * c){
        DReal a = (M + m)/2;
        DReal fa = len_func(a);
        if (fa > c) M = a, fM = fa;
        else m = a, fm = fa;
    }
    DReal sigma = (M + m)/2;
    return {BezierCurve<3, Vector>(std::array<Vector, 4>{C, C + sigma*(B-C), C + sigma/(D + sigma)*(A-C), A}), sigma};
} 

void CilindricContactSurface::LeafPlaneContactSDF::set_lambda(DReal lambda){
    m_lambda = lambda;
}

void CilindricContactSurface::LeafContactSDF::set_lambda(DReal lambda){ 
    m_lambda = lambda; 
    for (int i = 0; i < 2; ++i){
        m_current_lines[i] = lambda * m_main_lines[i] + (1-lambda)*BezierCurve<1, Vector>(std::array<Vector, 2>{opposit_point, m_main_lines[i].B[3]});
        m_current_derivatives[i] = lambda * m_main_derivatives[i] + (1 - lambda)*(m_main_lines[i].B[3] - opposit_point);
        m_current_lines[i] = m_current_lines[i] - (m_current_lines[i]*m_n)*m_n;
        m_current_derivatives[i] = m_current_derivatives[i] - (m_current_derivatives[i]*m_n)*m_n;
        m_current_dderivs[i] = m_current_derivatives[i].derivative();
        m_current_dist_func0[i] = m_current_lines[i] * m_current_derivatives[i];
    }
}
void CilindricContactSurface::LeafContactSDF::make_shift(const Vector& v){
    for (int i = 0; i < 2; ++i){
        m_main_lines[i] = m_main_lines[i] + v;
    }
    opposit_point += v;
    set_lambda(m_lambda);
}
SignedDistanceField::SDF CilindricContactSurface::LeafPlaneContactSDF::operator()(const Vector& P) const{
    Vector X = P - (P*m_n)*m_n;
    Vector Cl = (1 - m_lambda) * m_Cs + m_lambda * m_C;
    std::array<Vector, 2> Al{m_A, m_B};

    Vector p = Cl;
    DReal d2 = (X - p).squared_length();
    int ncurve = 0;
    bool is_apex = true;
    auto set_position = [&](Vector pp, int _ncurve, bool _is_apex = false){
        DReal dd2 = (X-pp).squared_length();
        if (dd2 < d2)
            p = pp, d2 = dd2, ncurve = _ncurve, is_apex = _is_apex;
    };
    for (int i = 0; i < 2; ++i){
        DReal t = (X - Cl)*(Al[i] - Cl)/(Al[i] - Cl).squared_length();
        if (t >= 0) set_position(Cl + (Al[i]-Cl)*t, i);
    }
    SDF res;
    if (is_apex){
        if (d2 > 0){
            res.grad_sdf = (X - p);
            bool is_positive =   (res.grad_sdf * CGAL::cross_product(m_n, Cl - Al[0]) >= 0)
                               && (res.grad_sdf * CGAL::cross_product(m_n, Cl - Al[1]) <= 0);
            res.sdf = (is_positive ? 1 : -1) * sqrt(d2);
            res.grad_sdf = (X - p) / res.sdf;
            Vector f = CGAL::cross_product(res.grad_sdf, m_n);
            DReal K = 1.0 / abs(res.sdf);
            if (std::isinf(K) || std::abs(K) > 1e7) K = ((K > 0) ? 1 : -1)*1e7;
            for (int i = 0; i < 3; ++i) for (int j = i; j < 3; ++j)
                res.grad_grad_sdf(i,j) = K*f[i]*f[j];
        } else {
            res.grad_sdf = 0.5 * (Cl - Al[0]) + 0.5 * (Cl - Al[1]);
            res.grad_sdf /= sqrt(res.grad_sdf.squared_length());
            res.sdf = 0;
        }
    } else {
        res.grad_sdf = CGAL::cross_product(m_n, Al[ncurve] - Cl);
        bool is_positive = (X - p) * res.grad_sdf >= 0;
        res.grad_sdf /= sqrt(res.grad_sdf.squared_length());
        if (ncurve == 1) {
            is_positive = !is_positive;
            res.grad_sdf *= -1;
        }
        res.sdf = (is_positive ? 1 : -1) * sqrt(d2);
    }
    return res;
}

SignedDistanceField::SDF CilindricContactSurface::LeafContactSDF::operator()(const Vector& P) const{
    Vector X = P - (P*m_n)*m_n;
    Vector Cl = m_current_lines[0].B[0];
    std::array<Vector, 2> Al{m_current_lines[0].B[3], m_current_lines[1].B[3]};

    Vector p = Cl;
    DReal d2 = (X - p).squared_length(), tt = -1;
    int ncurve = 0;
    bool is_apex = true;
    auto set_position = [&](Vector pp, int _ncurve, DReal t = -1, bool _is_apex = false){
        DReal dd2 = (X-pp).squared_length();
        if (dd2 < d2)
            p = pp, d2 = dd2, tt = t, ncurve = _ncurve, is_apex = _is_apex;
    };

    for (int i = 0; i < 2; ++i){
        bool is_linear = (m_lambda == 0.0 || m_sigma[i] == 0);
        DReal t = (X - Cl)*(Al[i] - Cl)/(Al[i] - Cl).squared_length();
        if (t >= 1) set_position(Cl + (Al[i]-Cl)*t, i, t);
        else if (t >= 0 && is_linear) set_position(Cl + (Al[i]-Cl)*t, i, t);
        else set_position(Al[i], i);
        if (!is_linear){
            BezierCurve<5, DReal> m_current_dist_func = m_current_dist_func0[i] - X*m_current_derivatives[i];
            auto coefs = m_current_dist_func.polynomial();
            Eigen::VectorXd _coefs, _re, _im; _coefs.resize(6);
            for (int k = 0; k < 6; ++k) _coefs[k] = coefs[5-k];
            bool isSolved = rpoly_plus_plus::FindPolynomialRootsJenkinsTraub(_coefs, &_re, &_im);
            if (!isSolved) {
                const int n = 6;
                gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(n);
                double z[2 * (n - 1)];
                if (gsl_poly_complex_solve (coefs.data(), n, w, z) == 0) {
                    isSolved = true;
                    gsl_poly_complex_workspace_free (w);
                    _re.resize(n-1), _im.resize(n-1);
                    for (int k = 0; k < n-1; ++k) _re[i] = z[2*i], _im[i] = z[2*i+1];
                    continue; 
                }
            } 
            if (isSolved){
                for (int k = 0; k < _re.size(); ++k) 
                    if (_im[k] == 0 && _re[k] > 0 && _re[k] < 1) 
                        set_position(m_current_lines[i].eval(_re[k]), i, _re[k]);
            } else
                std::cerr << "Warning: rpoly is failed" << std::endl;
        }
    }
    SDF res;
    if (is_apex){
        if (d2 > 0){
            res.grad_sdf = (X - p);
            bool is_positive =   (res.grad_sdf * CGAL::cross_product(m_n, m_current_lines[0].B[1] - m_current_lines[0].B[0]) >= 0)
                               && (res.grad_sdf * CGAL::cross_product(m_n, m_current_lines[1].B[1] - m_current_lines[0].B[1]) <= 0);
            res.sdf = (is_positive ? 1 : -1) * sqrt(d2);
            res.grad_sdf = (X - p) / res.sdf;
            Vector f = CGAL::cross_product(res.grad_sdf, m_n);
            DReal K = 1.0 / abs(res.sdf);
            if (std::isinf(K) || std::abs(K) > 1e7) K = ((K > 0) ? 1 : -1)*1e7;
            for (int i = 0; i < 3; ++i) for (int j = i; j < 3; ++j)
                res.grad_grad_sdf(i,j) = K*f[i]*f[j];
        } else {
            res.grad_sdf = 0.5 * (m_current_lines[0].B[1] - m_current_lines[0].B[0]) + 0.5 * (m_current_lines[1].B[1] - m_current_lines[1].B[0]);
            res.grad_sdf /= sqrt(res.grad_sdf.squared_length());
            res.sdf = 0;
        }
    } else {
        res.grad_sdf = (X - p);
        bool is_positive = true, is_curve = false;
        if (tt >= 1.0 || m_lambda == 0.0 || m_sigma[ncurve] == 0.0){
            res.grad_sdf = CGAL::cross_product(m_n, Al[ncurve] - Cl);
            is_positive = (X - p) * res.grad_sdf >= 0;
            res.grad_sdf /= sqrt(res.grad_sdf.squared_length());
            if (ncurve == 1) {
                is_positive = !is_positive;
                res.grad_sdf *= -1;
            }
            res.sdf = (is_positive ? 1 : -1) * sqrt(d2);
        }else {
            is_curve = true;
            auto dy = m_current_derivatives[ncurve].eval(tt);
            res.grad_sdf = CGAL::cross_product(m_n, dy);
            is_positive = (X - p) * res.grad_sdf >= 0;
            res.grad_sdf /= sqrt(res.grad_sdf.squared_length());
            if (ncurve == 1) {
                is_positive = !is_positive;
                res.grad_sdf *= -1;
            }
            res.sdf = (is_positive ? 1 : -1) * sqrt(d2);
            Vector f = CGAL::cross_product(res.grad_sdf, m_n);
            auto ddy = m_current_dderivs[ncurve].eval(tt);
            DReal r = (dy*dy - res.sdf * (res.grad_sdf*ddy));
            DReal sigma = (res.grad_sdf*ddy) / dy.squared_length();
            DReal K = 0.0;
            if (sigma != 0.0) {
                K = -(is_positive ? 1 : -1)*sigma / (1-res.sdf*sigma);
                //avoiding the singularity
                if (std::isnan(K)) K = 0.0;
                else if (std::isinf(K) || std::abs(K) > 1e7) K = ((K > 0) ? 1 : -1)*1e7;
            }
            for (int i = 0; i < 3; ++i) for (int j = i; j < 3; ++j)
                res.grad_grad_sdf(i,j) = K*f[i]*f[j];
        }
    }
    return res;
}

CilindricContactSurface::CilindricContactSurface(const std::array<Vector, 3>& commissure_points, const std::array<DReal, 3>& leaf_width): m_A{commissure_points}, m_l{leaf_width}{
    setup();
} 
CilindricContactSurface::CilindricContactSurface(std::array<Point, 3> commissure_points, const std::array<DReal, 3>& leaf_width): 
    CilindricContactSurface(std::array<Vector, 3>{commissure_points[0]-CGAL::ORIGIN, commissure_points[1]-CGAL::ORIGIN, commissure_points[2]-CGAL::ORIGIN}, leaf_width) { }