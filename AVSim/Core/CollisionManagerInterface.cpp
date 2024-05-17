//
// Created by Liogky Alexey on 01.06.2022.
//

#include "CollisionManagerInterface.h"

///\remark implementation of algorithm from Tang M. et al. Fast and exact continuous collision detection with bernstein sign classification

using namespace World3d;

template<int N>
static DReal DeCasteljau(DReal t0, const DReal* coefs){
    DReal tmp[N];
    std::copy(coefs, coefs + N, tmp);
    auto q = (1 - t0);
    for (int j = 1; j < N; ++j) {
        for (int k = 0; k < N - j; ++k) {
            tmp[k] = tmp[k] * q + tmp[k + 1] * t0;
        }
    }
    return tmp[0];
}

template<int N>
static DReal DeCasteljauSplit(DReal t0, const DReal* coefs, DReal* coefs1, DReal* coefs2){
    DReal tmp[N];
    std::copy(coefs, coefs + N, tmp);
    coefs1[0] = tmp[0];
    coefs2[N-1] = tmp[N-1];
    auto q = (1 - t0);
    for (int j = 1; j < N; ++j) {
        for (int k = 0; k < N - j; ++k) {
            tmp[k] = tmp[k] * q + tmp[k + 1] * t0;
        }
        coefs1[j] = tmp[0];
        coefs2[N-j-1] = tmp[N-j-1];
    }
    return tmp[0];
}

template<int N>
static bool isNegativeBezier(DReal* L){
    if (!((L[0] < 0) && (L[N-1] < 0))) return false;
    for (int i = 1; i < N - 1; ++i)
        if (!(L[i] <= 0)) return false;
    return true;
}

template<int N>
static bool isPositiveBezier(DReal* L){
    if (!((L[0] > 0) && (L[N-1] > 0))) return false;
    for (int i = 1; i < N - 1; ++i)
        if (!(L[i] >= 0)) return false;
    return true;
}

///Approximate solver for function f(x) on [a, b] such that
/// f(x) is strictly monotonic on [a, b], f(a)*f(b) < 0
/// f'(x) != 0 on (a, b)
/// f''(x) is strictly monotonic on [a, b] and f''(x) != 0 on (a, b)
template <class Ftype, class FDtype, class FDDtype, typename NT>
static NT quick_special_nonlin_solver(const Ftype& f, const FDtype& df, const FDDtype& ddf, NT a, NT b,
                                        NT ftol = (1e-7 + 10*std::numeric_limits<NT>::epsilon()),
                                        NT xtol = (1e-7 + 10*std::numeric_limits<NT>::epsilon())){
    if (a > b) std::swap(a, b);
    auto lerp = [a, b](double t){ return a*(1-t) + t*b; };
    int bin_its = 0, newton_its = 0;
    int f_eval = 0, df_eval = 0, ddf_eval = 0;
    NT f_left  = f(a), t_left = 0; ++f_eval;
    if (abs(f_left) < ftol) return f_left;
    NT f_right = f(b), t_right = 1; ++f_eval;
    if (abs(f_right) < ftol) return f_right;
    int sc = (f_left < f_right) ? 1 : -1;
    if (f_left * f_right >= 0) return (abs(f_left) < abs(f_right)) ? a : b;


    auto binary_search = [&](int maxits) -> bool {
        int i = 0;
        while (true){
            ++bin_its;
            NT t_try = (t_right + t_left) / 2;
            NT f_try = f(lerp(t_try)); ++f_eval;
            if (abs(f_try) < ftol) return t_right = t_left = t_try, true;
            if (sc * f_try < 0) f_left  = f_try, t_left  = t_try;
            else                f_right = f_try, t_right = t_try;

            if (t_right - t_left < xtol) return true;
            if (++i >= maxits) return false;
        }
    };

    NT init_border = 1e-5 + 10*std::numeric_limits<NT>::epsilon();
    NT f_try = f(lerp(init_border)); ++f_eval;
    if (f_try * sc > 0){
        f_right = f_try, t_right = init_border;
        binary_search(100);
        return lerp((t_left + t_right) / 2);
    }
    f_left = f_try, t_left = init_border;
    f_try = f(lerp(1-init_border)); ++f_eval;
    if (f_try * sc < 0){
        f_left = f_try, t_left = init_border;
        binary_search(100);
        return lerp((t_left + t_right) / 2);
    }

    NT ddf_left = ddf(lerp(t_left)), ddf_right = ddf(lerp(t_right)); ddf_eval += 2;
    bool choose_dd_side = (abs(ddf_left) > abs(ddf_right)); //choose left side?
    int dd_sg = (ddf_left + ddf_right > 0) ? 1 : -1;
    bool choose_d_side = (sc > 0) != (dd_sg > 0); //choose left side?
    auto max_abs_ddf = [&](){ ++ddf_eval; return abs(ddf(lerp(choose_dd_side ? t_left : t_right))); };
    auto min_abs_df = [&](){ ++df_eval; return abs(df(lerp((!choose_d_side) ? t_left : t_right))); };
    NT M = (choose_dd_side ? abs(ddf_left) : abs(ddf_right)) / min_abs_df() / (b - a);
    auto newton_search = [&](int maxits) -> bool{
        int i = 0;
        while (true){
            ++newton_its;
            NT t_try = (t_right + t_left) / 2;
            NT f_try = f(lerp(t_try)); ++f_eval;
            if (abs(f_try) < ftol) return t_right = t_left = t_try, true;
            NT df_left = df(lerp(t_left)), df_right = df(lerp(t_right)); df_eval += 2;
            NT  t1 = t_try - f_try / ((b - a) * df_left),
                t2 = t_try - f_try / ((b - a) * df_right);
            if (t1 > t2) std::swap(t1, t2);
            if (t1 > t_left) t_left = t1;
            if (t_left > t_right) t_left = t_right;
            if (t2 < t_right) t_right = t2;
            if (t_right < t_left) t_right = t_left;

            if (t_right - t_left < xtol) return true;
            if (++i >= maxits) return false;
        }
    };

    bool found = false;
    int it = 0, maxit = 15;
    while (it < maxit && !found){
        NT r = (t_right - t_left) * M / 4;
        if (r < 1) {
            found = newton_search(5);
        } else {
            bool make_binary = false;
            for (int k = 4; k >= 0; --k)
                if (r >= (1 << k)) {
                    found = binary_search(k+1);
                    make_binary = true;
                    break;
                }
            if (!make_binary) found = binary_search(5);
            M = (b - a) * max_abs_ddf() / min_abs_df();
        }
        ++it;
    }
    return lerp((t_left + t_right) / 2);
}
struct CCDQuery{
    DReal t_res = 0;
    bool has_collision = false;
    bool should_find_time = false;

    const Vector &m_a0, &m_a1, &m_a2, &m_a3;
    const Vector &m_at0, &m_at1, &m_at2, &m_at3;
    CCDQuery(const Vector& a0, const Vector& a1, const Vector& a2, const Vector& a3,
             const Vector& at0, const Vector& at1, const Vector& at2, const Vector& at3):
            m_a0{a0}, m_a1{a1}, m_a2{a2}, m_a3{a3}, m_at0{at0}, m_at1{at1}, m_at2{at2}, m_at3{at3}{}
    void ComputeCoplanarityCoefficients(DReal* coefs){
        auto abs_vec = [](const Vector& a) { return Vector(abs(a[0]), abs(a[1]), abs(a[2])); };
        auto mk_add_tol = [&](const Vector& a, const Vector& b, const Vector& r){ return std::numeric_limits<DReal>::epsilon() * (Vector(1, 1, 1) + abs_vec(a) + abs_vec(b) + abs_vec(r));};
        auto mk_dot_tol = [&](const Vector& a, const Vector& atol, const Vector& b, const Vector& btol) { return atol*abs_vec(b) + abs_vec(a)*btol + abs(a*b)*std::numeric_limits<DReal>::epsilon(); };

        auto v0 = m_at0 - m_a0, v1 = m_at1 - m_a1, v2 = m_at2 - m_a2;
        auto v0tol = mk_add_tol(m_at0, m_a0, v0), v1tol = mk_add_tol(m_at1, m_a1, v1), v2tol = mk_add_tol(m_at2, m_a2, v2);
        auto n0 = cross_product((m_a1 - m_a0),(m_a2 - m_a0)), n1 = cross_product((m_at1 - m_at0),(m_at2 - m_at0));
        auto nhat = (n0 + n1 - cross_product((v1 - v0),(v2 - v0))) / 2;
        coefs[0] = (m_a3 - m_a0) * n0;
        coefs[1] = (2*(m_a3 - m_a0)*nhat + (m_at3 - m_at0)*n0)/3;
        coefs[2] = (2*(m_at3 - m_at0)*nhat + (m_a3 - m_a0)*n1)/3;
        coefs[3] = (m_at3 - m_at0) * n1;

        auto n0tol = abs_vec(m_a1 - m_a0)*abs_vec(m_a2 - m_a0)*std::numeric_limits<DReal>::epsilon()*Vector(1, 1, 1);
        auto n1tol = abs_vec(m_at1 - m_at0)*abs_vec(m_at2 - m_at0)*std::numeric_limits<DReal>::epsilon()*Vector(1, 1, 1);
        auto nhattol = (n0tol + n1tol + abs_vec(v1 - v0)*abs_vec(v2 - v0)*std::numeric_limits<DReal>::epsilon()*Vector(1, 1, 1) + v0tol + v1tol) / 2;
        if (abs(coefs[0]) < mk_dot_tol(m_a3 - m_a0, mk_add_tol(m_a3, m_a0, m_a3 - m_a0), n0, n0tol))
            coefs[0] = 0;
        if (abs(coefs[1]) < (2*mk_dot_tol(m_a3 - m_a0, mk_add_tol(m_a3, m_a0, m_a3 - m_a0), nhat, nhattol) + mk_dot_tol(m_at3 - m_at0, mk_add_tol(m_at3, m_at0, m_at3 - m_at0), n0, n0tol))/3)
            coefs[1] = 0;
        if (abs(coefs[2]) < (2*mk_dot_tol(m_at3 - m_at0, mk_add_tol(m_at3, m_at0, m_a3 - m_a0), nhat, nhattol) + mk_dot_tol(m_a3 - m_a0, mk_add_tol(m_a3, m_a0, m_a3 - m_a0), n1, n1tol))/3)
            coefs[2] = 0;
        if (abs(coefs[3]) < mk_dot_tol(m_at3 - m_at0, mk_add_tol(m_at3, m_at0, m_at3 - m_at0), n1, n1tol))
            coefs[3] = 0;
    }
    bool ConservativeFilter(DReal* coefs){
        return (coefs[0] > 0 && coefs[1] >= 0 && coefs[2] >= 0 && coefs[3] > 0)
               || (coefs[0] < 0 && coefs[1] <= 0 && coefs[2] <= 0 && coefs[3] < 0);
    }
    void simplify_inequality(const DReal* Y/*[4]*/, const DReal* L/*[5]*/, DReal* P/*[3]*/){
        DReal K[5] = {Y[0], (Y[0] + 3*Y[1])/4, (Y[1] + Y[2])/2, (3*Y[2] + Y[3])/4, Y[3]};
        if (abs(Y[0]) >= abs(Y[3])) {
            DReal s[4] = {(L[1] * K[0] - L[0] * K[1]) * 4,
                          (L[2] * K[0] - L[0] * K[2]) * 2,
                          ((L[3] * K[0] - L[0] * K[3]) * 4) / 3,
                          (L[4] * K[0] - L[0] * K[4])};
            P[0] = (s[1] * Y[0] - s[0] * Y[1]) * 3;
            P[1] = ((s[2] * Y[0] - s[0] * Y[2]) * 3) / 2;
            P[2] = s[3] * Y[0] - s[0] * Y[3];
        } else {
            DReal s[4] = {(L[3] * K[4] - L[4] * K[3]) * 4,
                          (L[2] * K[4] - L[4] * K[2]) * 2,
                          ((L[1] * K[4] - L[4] * K[1]) * 4) / 3,
                          (L[0] * K[4] - L[4] * K[0])};
            P[2] = (s[1] * Y[3] - s[0] * Y[2]) * 3;
            P[1] = ((s[2] * Y[3] - s[0] * Y[1]) * 3) / 2;
            P[0] = s[3] * Y[3] - s[0] * Y[0];
        }
    }
    static void GetInequality(DReal* L/*[5]*/,
                              const Vector& a0, const Vector& a1, const Vector& a2, const Vector& a3,
                              const Vector& at0, const Vector& at1, const Vector& at2, const Vector& at3){
        auto abs_vec = [](const Vector& a) { return Vector(abs(a[0]), abs(a[1]), abs(a[2])); };
        auto nrm1 = [](const Vector& a) { return abs(a[0]) + abs(a[1]) + abs(a[2]); };
        auto v0 = at0 - a0, v1 = at1 - a1, v2 = at2 - a2, v3 = at3 - a3;
        auto    da10 = a1 - a0, da20 = a2 - a0, da13 = a1 - a3, da23 = a2 - a3,
                da10t = at1 - at0, da20t = at2 - at0, da13t = at1 - at3, da23t = at2 - at3;
        auto n0 = cross_product(da10,da20), n1 = cross_product(da10t,da20t), vn = cross_product((v1 - v0),(v2 - v0));
        if (nrm1(n0) < std::numeric_limits<DReal>::epsilon() * nrm1(da10) * nrm1(da20)) n0 = Vector(0, 0, 0);
        if (nrm1(n1) < std::numeric_limits<DReal>::epsilon() * nrm1(da10t) * nrm1(da20t)) n1 = Vector(0, 0, 0);
        if (nrm1(vn) < std::numeric_limits<DReal>::epsilon() * nrm1(v1 - v0) * nrm1(v2 - v0)) vn = Vector(0, 0, 0);
        auto m0 = cross_product(da13, da23), m1 = cross_product(da13t, da23t), vm = cross_product((v1 - v3),(v2 - v3));
        if (nrm1(m0) < std::numeric_limits<DReal>::epsilon() * nrm1(da13) * nrm1(da23)) m0 = Vector(0, 0, 0);
        if (nrm1(m1) < std::numeric_limits<DReal>::epsilon() * nrm1(da13t) * nrm1(da23t)) m1 = Vector(0, 0, 0);
        if (nrm1(vm) < std::numeric_limits<DReal>::epsilon() * nrm1(v1 - v3) * nrm1(v2 - v3)) vm = Vector(0, 0, 0);
        auto nhat = (n0 + n1 - vn) / 2;
        auto mhat = (m0 + m1 - vm) / 2;
        L[0] = m0*n0;
        L[1] = (m0*nhat + mhat*n0)/2;
        L[2] = (m0*n1 + 4 * mhat*nhat + m1*n0)/6;
        L[3] = (mhat*n1 + m1*nhat)/2;
        L[4] = m1*n1;
    }
    struct InPlaneData{
        Vector_2 a0s, at0s;
        Vector_2 a1s, at1s;
        Vector_2 a2s, at2s;
        Vector_2 a3s, at3s;
        InPlaneData(CCDQuery& c){
            auto abs_vec = [](const Vector& a) { return Vector(abs(a[0]), abs(a[1]), abs(a[2])); };
            auto mk_add_tol = [&](const Vector& a, const Vector& b, const Vector& r){ return std::numeric_limits<DReal>::epsilon() * (Vector(1, 1, 1) + abs_vec(a) + abs_vec(b) + abs_vec(r));};
            auto mk_dot_tol = [&](const Vector& a, const Vector& atol, const Vector& b, const Vector& btol) { return atol*abs_vec(b) + abs_vec(a)*btol + abs(a*b)*std::numeric_limits<DReal>::epsilon(); };
            auto compare_vectol = [&](const Vector& a, const Vector& atol){
                for (int k = 0; k < 3; ++k) if (abs(a[k]) > atol[k]) return true;
                return false;
            };
            auto n0 = cross_product((c.m_a1 - c.m_a0),(c.m_a2 - c.m_a0)), n1 = cross_product((c.m_at1 - c.m_at0),(c.m_at2 - c.m_at0));
            auto n0tol = std::numeric_limits<DReal>::epsilon()*(abs_vec(c.m_a1) + 2*abs_vec(c.m_a0) + abs_vec(c.m_a2) + abs_vec(c.m_a1 - c.m_a0)*abs_vec(c.m_a2 - c.m_a0)*Vector(1, 1, 1));
            auto n1tol = std::numeric_limits<DReal>::epsilon()*(abs_vec(c.m_at1) + 2*abs_vec(c.m_at0) + abs_vec(c.m_at2) + abs_vec(c.m_at1 - c.m_at0)*abs_vec(c.m_at2 - c.m_at0)*Vector(1, 1, 1));
            auto cmp0 = compare_vectol(n0, n0tol), cmp1 = compare_vectol(n1, n1tol);
            if (cmp0 && cmp1){
                n0 = n0 / sqrt(n0.squared_length()), n1 = n1 / sqrt(n1.squared_length());
            } else if (!cmp0 && !cmp1){
                auto e0 = c.m_a1 - c.m_a0, e1 = c.m_at1 - c.m_at0;
                e0 = e0 / sqrt(e0.squared_length()), e1 = e1 / sqrt(e1.squared_length());
                a0s = CGAL::NULL_VECTOR, at0s = CGAL::NULL_VECTOR;
                a1s = Vector_2((c.m_a1 - c.m_a0)*e0, 0), at1s = Vector_2((c.m_at1 - c.m_at0)*e1, 0);
                a2s = Vector_2((c.m_a2 - c.m_a0)*e0, 0), at2s = Vector_2((c.m_at2 - c.m_at0)*e1, 0);
                DReal d3s = (c.m_a3 - c.m_a0)*e0, dt3s = (c.m_at3 - c.m_at0)*e1;
                DReal r3s = sqrt(((c.m_a3 - c.m_a0) - d3s*e0).squared_length()), rt3s = sqrt(((c.m_at3 - c.m_at0) - dt3s*e1).squared_length());
                bool sgt = true;
                //if (rt3s > dt3s * std::numeric_limits<DReal>::epsilon()){
                    //TODO: check of sign changing (if (c.m_a1 - c.m_a0) x (c.m_a3 - c.m_a0)(t) = 0 has root)
                    //here we may loose root, but probability of such scenario in 3D is very low, so we neglect them
                //}
                a3s = Vector_2(d3s, r3s), at3s = Vector_2(dt3s, sgt ? rt3s : -rt3s);
                return;
            } else if (cmp0) {
                n0 = n0 / sqrt(n0.squared_length());
                n1 = 0.5 * ( cross_product((c.m_at1 - c.m_at0),(c.m_a2 - c.m_a0)) + cross_product((c.m_a1 - c.m_a0),(c.m_at2 - c.m_at0)));
                n1 = n1 / sqrt(n1.squared_length());
            } else if (cmp1) {
                n1 = n1 / sqrt(n0.squared_length());
                n0 = -0.5 * ( cross_product((c.m_at1 - c.m_at0),(c.m_a2 - c.m_a0)) + cross_product((c.m_a1 - c.m_a0),(c.m_at2 - c.m_at0)));
                n0 = n0 / sqrt(n0.squared_length());
            }
            auto e0 = cross_product(n0, n1);
            auto e0tol = n0tol + n1tol + std::numeric_limits<DReal>::epsilon()*(abs_vec(n0)*abs_vec(n1)*Vector(1, 1, 1));
            if (!compare_vectol(e0, e0tol)){
                int k = 0;
                DReal dk = abs(n0[0]);
                for (int i = 1; i < 3; ++i) if (abs(n0[k]) < dk){ dk = abs(n0[k]), k = i; }
                std::array<DReal, 3> vals = {0, 0, 0}; vals[k] = 1;
                e0 = Vector(vals[0], vals[1], vals[2]);
                e0 = e0 - (e0*n0)*n0;
            }
            e0 = e0 / sqrt(e0.squared_length());
            Vector e1 = e0;
            auto f0 = cross_product(n0,e0), f1 = cross_product(n1,e1);
            a0s = Vector_2((c.m_a0 - c.m_a3)*e0, (c.m_a0 - c.m_a3)*f0), at0s = Vector_2((c.m_at0 - c.m_at3)*e1, (c.m_at0 - c.m_at3)*f1);
            a1s = Vector_2((c.m_a1 - c.m_a3)*e0, (c.m_a1 - c.m_a3)*f0), at1s = Vector_2((c.m_at1 - c.m_at3)*e1, (c.m_at1 - c.m_at3)*f1);
            a2s = Vector_2((c.m_a2 - c.m_a3)*e0, (c.m_a2 - c.m_a3)*f0), at2s = Vector_2((c.m_at2 - c.m_at3)*e1, (c.m_at2 - c.m_at3)*f1);
            a3s = CGAL::NULL_VECTOR, at3s = CGAL::NULL_VECTOR;
        }
    };
    static void GetInequality2D(DReal* L/*[3]*/,
                              const Vector_2& a0, const Vector_2& a1, const Vector_2& a2, const Vector_2& a3,
                              const Vector_2& at0, const Vector_2& at1, const Vector_2& at2, const Vector_2& at3){
        auto cross2d = [](const Vector_2& a, const Vector_2& b){ return a[0]*b[1] - a[1]*b[0];};

        L[0] = cross2d(a1 - a3, a2 - a3);
        L[1] = (cross2d(at1 - at3, a2 - a3) + cross2d(a1 - a3, at2 - at3)) / 2;
        L[2] = cross2d(at1 - at3, at2 - at3);
    }
    struct VFInsideTest{
        static void exec(CCDQuery& c, DReal* L/*[3*5]*/){
            CCDQuery::GetInequality(L + 0*5, c.m_a0, c.m_a1, c.m_a2, c.m_a3, c.m_at0, c.m_at1, c.m_at2, c.m_at3);
            CCDQuery::GetInequality(L + 1*5, c.m_a1, c.m_a2, c.m_a0, c.m_a3, c.m_at1, c.m_at2, c.m_at0, c.m_at3);
            CCDQuery::GetInequality(L + 2*5, c.m_a2, c.m_a0, c.m_a1, c.m_a3, c.m_at2, c.m_at0, c.m_at1, c.m_at3);
        }
        static void exec2D(InPlaneData& c, DReal* L/*[3*3]*/){
            CCDQuery::GetInequality2D(L + 0*3, c.a0s, c.a1s, c.a2s, c.a3s, c.at0s, c.at1s, c.at2s, c.at3s);
            CCDQuery::GetInequality2D(L + 1*3, c.a1s, c.a2s, c.a0s, c.a3s, c.at1s, c.at2s, c.at0s, c.at3s);
            CCDQuery::GetInequality2D(L + 2*3, c.a2s, c.a0s, c.a1s, c.a3s, c.at2s, c.at0s, c.at1s, c.at3s);
        }
    };
    struct EEInsideTest{
        static void exec(CCDQuery& c, DReal* L/*[3*5]*/){
            CCDQuery::GetInequality(L + 0*5, c.m_a0, c.m_a1, c.m_a2, c.m_a3, c.m_at0, c.m_at1, c.m_at2, c.m_at3);
            CCDQuery::GetInequality(L + 1*5, c.m_a1, c.m_a2, c.m_a0, c.m_a3, c.m_at1, c.m_at2, c.m_at0, c.m_at3);
            CCDQuery::GetInequality(L + 2*5, c.m_a2, c.m_a0, c.m_a1, c.m_a3, c.m_at2, c.m_at0, c.m_at1, c.m_at3);
            for (int k = 0; k < 5; ++k) L[2*5 + k] *= -1;
        }
        static void exec2D(InPlaneData& c, DReal* L/*[3*3]*/){
            CCDQuery::GetInequality2D(L + 0*3, c.a0s, c.a1s, c.a2s, c.a3s, c.at0s, c.at1s, c.at2s, c.at3s);
            CCDQuery::GetInequality2D(L + 1*3, c.a1s, c.a2s, c.a0s, c.a3s, c.at1s, c.at2s, c.at0s, c.at3s);
            CCDQuery::GetInequality2D(L + 2*3, c.a2s, c.a0s, c.a1s, c.a3s, c.at2s, c.at0s, c.at1s, c.at3s);
            for (int k = 0; k < 3; ++k) L[2*3 + k] *= -1;
        }
    };
    static int getLinearSignAtEqTypeC(const DReal* Y, const DReal* L){
        if (L[0]*L[1] >= 0) return ((L[0] + L[1]) >= 0 ? 1 : -1);
        DReal ts = L[0] / (L[0] - L[1]);
        if (ts > 1) ts = 1;
        if (ts < 0) ts = 0;
        if (abs(Y[0]) >= abs(Y[3])){
            if (DeCasteljau<4>(ts, Y)*Y[0] >= 0) return (L[1] >= 0) ? 1 : -1;
            else return (L[0] >= 0) ? 1 : -1;
        }
        else {
            if (DeCasteljau<4>(ts, Y)*Y[3] >= 0) return (L[0] >= 0) ? 1 : -1;
            else return (L[1] >= 0) ? 1 : -1;
        }
    }
    static std::pair<int, int> getLinearSignAtEqTypeB(const DReal* Y, const DReal* L){
        if (L[0]*L[1] >= 0) return ((L[0] + L[1]) >= 0 ? std::pair<int, int>{1, 1} : std::pair<int, int>{-1, -1});
        DReal ts = L[0] / (L[0] - L[1]);
        if (ts > 1) ts = 1;
        if (ts < 0) ts = 0;
        auto yt = DeCasteljau<4>(ts, Y);
        auto sign = [](auto x) -> int { return (x >= 0) ? 1 : -1; };
        if (!(yt*(Y[0]+Y[3]) > 0)){
            return {sign(L[0]), sign(L[1])};
        } else {
            DReal dY[3] = {3 * (Y[1] - Y[0]), 3 * (Y[2] - Y[1]), 3 * (Y[3] - Y[2])};
            auto dyt = DeCasteljau<3>(ts, dY);
            if (sign(dyt) == sign(dY[0])) return {sign(L[1]), sign(L[1])};
            else return {sign(L[0]), sign(L[0])};
        }
    }
    //return number of proper root
    template<int nIneq>
    static int PerformInequalityTest(int nroots, const DReal* Y/*[4]*/, const DReal* P/*[nIneq*3]*/){
        if (nroots == 2) {
            int res = 3;
            for (int k = 0; k < nIneq && res > 0; ++k) {
                DReal L[2], K[2];
                if (P[3*k + 0] >= 0 && P[3*k + 1] >= 0 && P[3*k + 2] >= 0) continue;
                if (P[3*k + 0] < 0 && P[3*k + 1] <= 0 && P[3*k + 2] < 0) return 0;
                if (P[3*k + 0] == 0 && Y[0] == 0){
                    if ((P[3*k + 1] < 0) || (P[3*k + 1] == 0 && P[3*k + 2] < 0)) res &= ~(1);
                }
                if (!DecomposePoly(Y, P + k*3, L, K)){
                    DReal ts = P[3*k+0] / (P[3*k+2] - P[3*k+0]);
                    if (ts > 1) ts = 1;
                    else if (ts < 0) ts = 0;
                    auto yt = DeCasteljau<4>(ts, Y);
                    auto sign = [](auto x) -> int { return (x >= 0) ? 1 : -1; };
                    if (sign(yt) != sign(Y[0])){
                        if (sign(P[3*k+0]) < 0) res &= ~(1);
                        if (sign(P[3*k+2]) < 0) res &= ~(2);
                    } else {
                        DReal dY[3] = {3 * (Y[1] - Y[0]), 3 * (Y[2] - Y[1]), 3 * (Y[3] - Y[2])};
                        auto dyt = DeCasteljau<3>(ts, dY);
                        if (sign(dyt) == sign(dY[0])){
                            if (sign(P[3*k+2]) < 0) res = 0;
                        } else 
                            if (sign(P[3*k+0]) < 0) res = 0;
                    }
                } else { 
                    auto rL = getLinearSignAtEqTypeB(Y, L);
                    auto rK = getLinearSignAtEqTypeB(Y, K);
                    if (!(rL.first * rK.first <= 0)) res &= ~(1);
                    if (!(rL.second * rK.second <= 0)) res &= ~(2);
                }
            }
            if (res & 1) return 1;
            else if (res & 2) return 2;
            else return 0;
        } else if (nroots == 1){
            bool res = true;
            for (int k = 0; k < nIneq && res; ++k) {
                DReal L[2], K[2];
                if (P[3*k + 0] >= 0 && P[3*k + 1] >= 0 && P[3*k + 2] >= 0) continue;
                if (P[3*k + 0] < 0 && P[3*k + 1] <= 0 && P[3*k + 2] < 0) return 0;
                if (P[3*k + 0] == 0 && Y[0] == 0){
                    if ((P[3*k + 1] < 0) || (P[3*k + 1] == 0 && P[3*k + 2] < 0)) return 0;
                    else continue;
                }
                if (!DecomposePoly(Y, P + k*3, L, K)){
                    DReal ts = P[3*k+0] / (P[3*k+2] - P[3*k+0]);
                    if (ts > 1) ts = 1;
                    else if (ts < 0) ts = 0;
                    auto yt = DeCasteljau<4>(ts, Y);
                    auto sign = [](auto x) -> int { return (x >= 0) ? 1 : -1; };
                    if (sign(yt) == sign(Y[0])){
                        if (sign(P[3*k+2]) < 0) res = false; 
                    } else 
                        if (sign(P[3*k+0]) < 0) res = false; 
                } else {
                    auto rL = getLinearSignAtEqTypeC(Y, L);
                    auto rK = getLinearSignAtEqTypeC(Y, K);
                    res &= (rL * rK <= 0);
                }
            }
            return res ? 1 : 0;
        }
        return 0;
    }
    void solve_equation(DReal* Y/*[4]*/, DReal t0, DReal t1, char type, int nroot, int root_id){
        if (!should_find_time) return;
        DReal dY[3] = {3 * (Y[1] - Y[0]), 3 * (Y[2] - Y[1]), 3 * (Y[3] - Y[2])};
        DReal ddY[2] = {6*(Y[2] - 2*Y[1] + Y[0]), 6*(Y[3] - 2 * Y[2] + Y[1])};
        auto f_func = [Y](DReal t) {return DeCasteljau<4>(t, Y); };
        auto df_func = [dY](DReal t) {return DeCasteljau<3>(t, dY); };
        auto ddf_func = [ddY](DReal t) {return DeCasteljau<2>(t, ddY); };
        if (type == 1){
            DReal nrm = abs(Y[3] - Y[0]);
            t_res = quick_special_nonlin_solver(f_func, df_func, ddf_func, t0, t1, nrm*(1e-7 + 10*std::numeric_limits<DReal>::epsilon()));
        } else {
            auto dddf_func = [K = ddY[1] - ddY[0]](DReal t){ return K; };
            DReal d_nrm = abs(dY[2] - dY[0]);
            DReal q_res = quick_special_nonlin_solver(df_func, ddf_func, dddf_func, t0, t1, d_nrm * (1e-7 + 10*std::numeric_limits<DReal>::epsilon()));
            auto Yq = f_func(q_res);
            if (nroot == 2) {
                if (root_id == 1) {
                    DReal nrm = abs(Yq - Y[0]);
                    t_res = quick_special_nonlin_solver(f_func, df_func, ddf_func, t0, q_res,
                                                        nrm * (1e-7 + 10 * std::numeric_limits<DReal>::epsilon()));
                } else {
                    DReal nrm = abs(Y[1] - Yq);
                    t_res = quick_special_nonlin_solver(f_func, df_func, ddf_func, q_res, t1,
                                                        nrm * (1e-7 + 10 * std::numeric_limits<DReal>::epsilon()));
                }
            } else { //nroot == 1
                bool sgddY = (ddY[0] + ddY[1]) > 0;
                bool sgY0 = (Y[0] > 0);
                if (sgddY == sgY0){
                    DReal nrm = abs(Yq - Y[0]);
                    t_res = quick_special_nonlin_solver(f_func, df_func, ddf_func, t0, q_res,
                                                        nrm * (1e-7 + 10 * std::numeric_limits<DReal>::epsilon()));
                } else {
                    DReal nrm = abs(Y[1] - Yq);
                    t_res = quick_special_nonlin_solver(f_func, df_func, ddf_func, q_res, t1,
                                                        nrm * (1e-7 + 10 * std::numeric_limits<DReal>::epsilon()));
                }
            }
        }
    }
    //return true on antiinterval
    static bool FindNonNegativeSqrFuncPart(DReal* P/*[3]*/, DReal& t0, DReal& t1){
        if (P[0] >= 0 && P[1] >= 0 && P[2] >= 0) return true;
        if (P[0] < 0 && P[1] <= 0 && P[2] < 0) { t1 = t0 - 1; return true; }
        DReal a = P[0] - 2*P[1] + P[2], k = P[1] - P[0], c = P[0];
        if (a != 0){
            DReal D = k*k - a*c;
            if (D <= 0) {
                if (a < 0) t1 = t0 - 1;
                return true;
            }
            DReal sqrtD = sqrt(D);
            DReal x0 = (-k - sqrtD) / a, x1 = (-k + sqrtD) / a;
            if (x0 > x1) std::swap(x0, x1);
            if (a < 0){
                if (x0 > t0) t0 = x0;
                if (x1 < t1) t1 = x1;
                return true;
            } else {
                if (x0 >= t1 || x1 <= t0) return true;
                if (x1 > t1){
                    if (x0 < t0) t1 = t0 - 1;
                    else t1 = x0;
                    return true;
                } else if (x0 < t0){
                    if (x1 > t1) t1 = t0 - 1;
                    else t0 = x1;
                    return true;
                } else {
                    t0 = x0;
                    t1 = x1;
                    return false;
                }
            }
        } else {
            if (P[0] >= 0 && 2*P[1] - P[0] >= 0) return true;
            if (P[0] < 0 && 2*P[1] - P[0] <= 0) { t1 = t0 - 1; return true; };
            if (P[1] - P[0] != 0){
                DReal ts = P[0] / (P[0] - P[1]) / 2;
                if (ts < t0 || ts > t1) {
                    if (P[0] + P[1] < 0) t1 = t0 - 1;
                    return true;
                } else {
                    if (P[0] < 0) t0 = ts;
                    else if (P[0] > 0) t1 = ts;
                    else if (P[1] < 0) t1 = t0 - 1;
                    return true;
                }
            } else{
                if (P[0] < 0) t1 = t0 - 1;
                return true;
            }
        }
    }
    template<int NIneq, typename GetIneqs>
    void Perform2DTest(){
        InPlaneData ipd(*this);
        DReal L[NIneq][3];
        GetIneqs::exec2D(ipd, reinterpret_cast<DReal *>(L));
        auto check_k = [&](int j){ bool is_ok = true; for (int k = 0; k < NIneq; ++k) is_ok &= (L[k][j] > 0); return is_ok; };
        if (check_k(0)){ has_collision = true, t_res = 0; return; }
        DReal t0 = 0, t1 = 1;
        std::array<std::pair<DReal, DReal>, NIneq> anti_intervals;
        int j = 0;
        for (int i = 0; i < NIneq; ++i){
            DReal tt0 = t0, tt1 = t1;
            bool is_interval = FindNonNegativeSqrFuncPart(L[i], tt0, tt1);
            if (is_interval){
                t0 = tt0, t1 = tt1;
                if (t0 > t1) return;
            } else {
                anti_intervals[j++] = std::pair<DReal, DReal>{tt0, tt1};
            }
        }
        std::array<DReal, 2*NIneq+2> brds; std::fill(brds.begin(), brds.end(), 2);
        brds[0] = t0, brds[1] = t1;
        for (int l = 0; l < j; ++l) {
            if (anti_intervals[l].first > t0) brds[2 + 2*l] = anti_intervals[l].first;
            if (anti_intervals[l].second < t1)brds[2 + 2*l + 1] = anti_intervals[l].second;
        }
        std::sort(brds.data(), brds.data() + 2 + 2 * j);
        for (int r = 0; r < 2 + 2 * j; ++r){
            bool is_part = true;
            if (brds[r] < t0 || brds[r] > t1) is_part = false;
            for (int l = 0; l < j && is_part; ++l)
                is_part &= !((brds[r] > anti_intervals[l].first) && (brds[r] < anti_intervals[l].second));
            if (is_part){  has_collision = true, t_res = brds[r]; return; }
        }

        return;
    }

    template<int NIneq, typename GetIneqs>
    void PerformTest(){
        has_collision = false; t_res = 0;
        DReal coefs[4];
        ComputeCoplanarityCoefficients(coefs);
        if (ConservativeFilter(coefs)) return;
        if (coefs[0] == 0 && coefs[3] == 0){
            if (coefs[1] == 0 && coefs[2] == 0) {
                Perform2DTest<NIneq, GetIneqs>();
                return;
            }
            DReal L[NIneq][5];
            GetIneqs::exec(*this, reinterpret_cast<DReal *>(L));
            auto check_k = [&](int j){ bool is_ok = true; for (int k = 0; k < NIneq; ++k) is_ok &= (L[k][j] > 0); return is_ok; };
            if (check_k(0)) { has_collision = true, t_res = 0; return; }
            if (coefs[1] * coefs[2] < 0){
                DReal ts = coefs[1] / (coefs[1] - coefs[2]);
                bool is_ok = true;
                for (int k = 0; k < NIneq; ++k)
                    is_ok &= (DeCasteljau<5>(ts, L[k]) >= 0);
                if (is_ok) { has_collision = true, t_res = ts; return; }
            }
            if (check_k(4)) { has_collision = true, t_res = 1; return; }

            return;
        }
        DReal r0 = coefs[2] - 2*coefs[1]+coefs[0], r1 = coefs[3] - 2*coefs[2] + coefs[1];
#define EPS std::numeric_limits<DReal>::epsilon()
        if (abs(r0) < (abs(coefs[2]) + 2*abs(coefs[1]) + abs(coefs[0]))*EPS && abs(r1) < (abs(coefs[3]) + 2*abs(coefs[2]) + abs(coefs[1]))*EPS){
            //ddY is zero
            if (coefs[0]*coefs[3] > 0) return;
            DReal ts = coefs[0] / (coefs[0] - coefs[3]);
            if (ts < 0) ts = 0;
            else if (ts > 1) ts = 1;
            DReal L[NIneq][5];
            GetIneqs::exec(*this, reinterpret_cast<DReal *>(L));
            for (int k = 0; k < NIneq; ++k){
                if (DeCasteljau<5>(ts, L[k]) < 0) return;
            }
            has_collision = true; t_res = ts;
        }
#undef EPS
        if (r0 * r1 < 0){
            DReal t_s = r0 / (r0 - r1);
            DReal coefs1[4], coefs2[4];
            DReal L1[NIneq][5], L2[NIneq][5], P[NIneq][3];
            DeCasteljauSplit<4>(t_s, coefs, coefs1, coefs2);
            auto res1 = PerformCaseBCCoplanarityTest(coefs1);
            if (res1.ncols > 0){
                GetIneqs::exec(*this, reinterpret_cast<DReal *>(L1));
                for (int k = 0; k < NIneq; ++k)
                    DeCasteljauSplit<5>(t_s, L1[k], L1[k], L2[k]);
            }
            if (res1.ncols > 0){
                for (int k = 0; k < NIneq; ++k){ 
                    if (isPositiveBezier<5>(L1[k])) P[k][0] = P[k][1] = P[k][2] = 1;
                    else if (isNegativeBezier<5>(L1[k])) P[k][0] = P[k][1] = P[k][2] = -1;
                    else simplify_inequality(coefs1, L1[k], P[k]); 
                }
                int nroot = PerformInequalityTest<NIneq>(res1.ncols, coefs1, reinterpret_cast<const DReal *>(P));
                if (nroot > 0 && nroot <= res1.ncols){
                    has_collision = true;
                    solve_equation(coefs1, 0, 1, res1.cid, res1.ncols, nroot);
                    t_res *= t_s;
#ifndef NO_COLLISION_ZERO_CHECK
                    if (coefs1[0] == 0){
                        has_collision = false;
                    } else 
                        return;
#else
                    return;
#endif
                }
            }
            auto res2 = PerformCaseBCCoplanarityTest(coefs2);
            if (res1.ncols <= 0 && res2.ncols > 0){
                GetIneqs::exec(*this, reinterpret_cast<DReal *>(L1));
                for (int k = 0; k < NIneq; ++k)
                    DeCasteljauSplit<5>(t_s, L1[k], L1[k], L2[k]);
            }
            if (res2.ncols > 0){
                for (int k = 0; k < NIneq; ++k) {
                    if (isPositiveBezier<5>(L2[k])) P[k][0] = P[k][1] = P[k][2] = 1;
                    else if (isNegativeBezier<5>(L2[k])) P[k][0] = P[k][1] = P[k][2] = -1;
                    else simplify_inequality(coefs2, L2[k], P[k]);
                }
                int nroot = PerformInequalityTest<NIneq>(res2.ncols, coefs2, reinterpret_cast<const DReal *>(P));
                if (nroot > 0 && nroot <= res2.ncols){
                    has_collision = true;
                    solve_equation(coefs2, 0, 1, res2.cid, res2.ncols, nroot);
                    t_res = t_s + t_res * (1 - t_s);

#ifndef NO_COLLISION_ZERO_CHECK
                    if (coefs2[0] == 0 && t_s == 0){
                        has_collision = false;
                    } else 
                        return;
#else
                    return;
#endif
                }
            }
        } else {
            DReal L[NIneq][5], P[NIneq][3];
            auto res = PerformCaseBCCoplanarityTest(coefs);
            if (res.ncols == 0) return;
            GetIneqs::exec(*this, reinterpret_cast<DReal *>(L));
            for (int k = 0; k < NIneq; ++k){ 
                if (isPositiveBezier<5>(L[k])) P[k][0] = P[k][1] = P[k][2] = 1;
                else if (isNegativeBezier<5>(L[k])) P[k][0] = P[k][1] = P[k][2] = -1;
                else simplify_inequality(coefs, L[k], P[k]); 
            }
            int nroot = PerformInequalityTest<NIneq>(res.ncols, coefs, reinterpret_cast<const DReal *>(P));
            if (nroot > 0 && nroot <= res.ncols){
                has_collision = true;
                solve_equation(coefs, 0, 1, res.cid, res.ncols, nroot);

#ifndef NO_COLLISION_ZERO_CHECK
                    if (coefs[0] == 0){
                        has_collision = false;
                    } else 
                        return;
#else
                    return;
#endif
            }
        }
    }

    void PerformVFTest(){
#ifndef NO_COLLISION_ZERO_CHECK
        { //special check for initial flat configuration
            auto sqdt = (GeomProjs::BaryEval(m_a0, m_a1, m_a2, GeomProjs::BaryCoord(m_a0, m_a1, m_a2, m_a3)) - m_a3).squared_length();
            auto sqd_cmp = (m_a0).squared_length() + (m_a1).squared_length() + (m_a2).squared_length() + (m_a3).squared_length();
            if (sqdt < sqd_cmp * std::numeric_limits<DReal>::epsilon()) {
                has_collision = true;
                t_res = 0;
                return;
            }
        }
#endif
        PerformTest<3, VFInsideTest>();
    }
    void PerformEETest(){
#ifndef NO_COLLISION_ZERO_CHECK
        { //special check for initial flat configuration
            Vector prj0, prj1; double w[2]; double sqdt;
            GeomProjs::SegmentsProjects(m_a0, m_a1, m_a2, m_a3, prj0, prj1, w[0], w[1], sqdt);
            auto sqd_cmp = (m_a0).squared_length() + (m_a1).squared_length() + (m_a2).squared_length() + (m_a3).squared_length();
            if (sqdt < sqd_cmp * std::numeric_limits<DReal>::epsilon()) {
                has_collision = true;
                t_res = 0;
                return;
            }
        }
#endif
        PerformTest<3, EEInsideTest>();
    }

    struct ComplTestRes{
        char ncols = 0;
        char cid = -1;
    };
    ComplTestRes PerformCaseBCCoplanarityTest(DReal* coefs){
        DReal r0 = coefs[1] - coefs[0], r1 = coefs[3] - coefs[2];
        if (r0 * r1 < 0){ //there is extreme point
            if (coefs[0] * coefs[3] < 0) { //there is single root
                return ComplTestRes{1, 0};
            } 
            if (coefs[0] == 0){
                if (coefs[1] == 0){
                    if (coefs[2]*coefs[3] >= 0)
                        return ComplTestRes{1, 0};
                    else 
                        return ComplTestRes{2, 0};
                } else {
                    if (coefs[1] * coefs[3] < 0) 
                        return ComplTestRes{2, 0};
                    else 
                        return ComplTestRes{1, 0};
                }
            }
            if (coefs[3] == 0){
                if (coefs[2] == 0){
                    if (coefs[0]*coefs[1] >= 0)
                        return ComplTestRes{1, 0};
                    else 
                        return ComplTestRes{2, 0};
                } else {
                   if (coefs[0] * coefs[2] < 0) 
                        return ComplTestRes{2, 0}; 
                    else 
                        return ComplTestRes{1, 0}; 
                }
            }
            {
                DReal dcoef[3] = {3 * (coefs[1] - coefs[0]), 3 * (coefs[2] - coefs[1]), 3 * (coefs[3] - coefs[2])};
                DReal Sc[2], Tc[2];
                DecomposePoly(coefs, dcoef, Sc, Tc);
                if (Tc[1] * Tc[0] > 0){ //Tc has no root on [0, 1]
                    if (coefs[0]*Tc[0] > 0) return ComplTestRes{0, 0};
                    else { //has 2 roots
                        return ComplTestRes{2, 0};
                    }
                } else {
                    DReal ts = Tc[0] / (Tc[0] - Tc[1]);
                    if (ts < 0) ts = 0;
                    if (ts > 1) ts = 1;
                    if (DeCasteljau<3>(ts,dcoef)*dcoef[0] > 0){
                        if (coefs[0]*Tc[1] > 0) return ComplTestRes{0, 0};
                        else { //has 2 root
                            return ComplTestRes{2, 0};
                        }
                    } else {
                        if (coefs[0]*Tc[0] > 0) return ComplTestRes{0, 0};
                        else {
                            return ComplTestRes{2, 0};
                        }
                    }
                }
            }
        } else {  //there is no extreme point
            if (coefs[0] * coefs[3] > 0){ //there is no root
                return ComplTestRes{0, 1};
            } else {
                return ComplTestRes{1, 1};
            }
        }
    }
    static void DecomposePoly10(const DReal* i/*[1]*/, DReal* u/*[2]*/, DReal* v/*[2]*/){
        u[0] = u[1] = 0;
        v[0] = v[1] = i[0];
    }
    static void DecomposePoly21(const DReal* i/*[2]*/, const DReal* j/*[1]*/, DReal* u/*[2]*/, DReal* v/*[2]*/){
        u[0] = u[1] = i[0] / j[0];
        v[0] = 0;
        v[1] = i[1] - i[0]; 
    }
    static void DecomposePoly32(const DReal* i/*[3]*/, const DReal* j/*[2]*/, DReal* u/*[2]*/, DReal* v/*[2]*/){
        auto a = (i[2] - 2*i[1] + i[0]) / (j[1] - j[0]);
        u[0] = a;
        u[1] = 2*a;
        v[0] = i[0] - u[0] * j[0];
        v[1] = i[2] - j[1] * u[1];
    }
    static void DecomposePoly43(const DReal* i/*[4]*/, const DReal* j/*[3]*/, DReal* u/*[2]*/, DReal* v/*[2]*/){
        DReal a = 2 * (j[1] - j[0]), b = j[0] - j[2], c = 2 * (j[1] - j[2]);
        DReal p = 3 * i[1] - 2 * i[0] - i[3], q = 3*i[2] - 2*i[3] - i[0];
        DReal kj = j[0] - 2*j[1] + j[2];
        DReal det = kj*kj;//c*a + b*b;
        u[0] = (c * p - b * q) / det;
        u[1] = (a * q + b * p) / det;
        v[0] = i[0] - u[0] * j[0];
        v[1] = i[3] - u[1] * j[2];
    }
    static bool DecomposePoly(const DReal* i/*[4]*/, const DReal* j/*[3]*/, DReal* u/*[2]*/, DReal* v/*[2]*/){
        DReal kj2 = j[0] - 2*j[1] + j[2], kj2_erc = 64*(abs(j[0]) + 2*abs(j[1]) + abs(j[1]));
        if (abs(kj2) > kj2_erc*std::numeric_limits<DReal>::epsilon()){
            DecomposePoly43(i, j, u, v);
            return true;
        } else {// j is linear polynomial
            DReal ki3 = i[3] - 3*i[2] + 3*i[1] - i[0], ki3_erc = 64*(abs(i[3]) + 3*abs(i[2]) + 3*abs(i[1]) + abs(i[0]));
            if (abs(ki3) < ki3_erc*std::numeric_limits<DReal>::epsilon()){ //i is quadratic polynomial
                DReal i1[3], j1[2];
                DReal kj1 = j[2] - j[0], kj1_erc = 64*(abs(j[2]) + abs(j[0]));
                if (abs(kj1) > kj1_erc*std::numeric_limits<DReal>::epsilon()){
                    i1[0] = i[0], i1[1] = (3*i[1] - i[0])/2, i1[2] = i[3];
                    j1[0] = j[0], j1[1] = j[2];
                    DecomposePoly32(i1, j1, u, v);
                    return true;
                } else { //j is constant
                    DReal ki2 = i[3] - 2*i[2] + i[1], ki2_erc = 64*(abs(i[3]) + 2*abs(i[2]) + abs(i[1]));
                    if (abs(ki2) < ki2_erc*std::numeric_limits<DReal>::epsilon()){ //i is linear
                        if (j[0] == 0) {
                            DReal ki1 = i[3] - i[0], ki1_erc = 64*(abs(i[3]) + abs(i[0]));
                            if (abs(ki1) < ki1_erc*std::numeric_limits<DReal>::epsilon()){//i is constant
                                i1[0] = i[0];
                                DecomposePoly10(i1, u, v);
                                return true;
                            } else return false;
                        } else {
                            i1[0] = i[0], i1[1] = i[3];
                            j1[0] = j[0];
                            DecomposePoly21(i1, j1, u, v);
                            return true;
                        }
                    } else 
                        return false;
                }
            } else //i is true cubic polynomial
                return false;
        }
    }
};

///macro NO_COLLISION_ZERO_CHECK turn off special check of case when t = 0 for initially flat configuration
bool GeomProjs::PerformVFCCD(const Vector& a0, const Vector& b0, const Vector& c0, const Vector& p0,
                             const Vector& a1, const Vector& b1, const Vector& c1, const Vector& p1, DReal* t){
    CCDQuery ccd(a0, b0, c0, p0, a1, b1, c1, p1);
    ccd.should_find_time = (t != nullptr);
    ccd.PerformVFTest();
    if (t != nullptr && ccd.has_collision) *t = ccd.t_res;
    return ccd.has_collision;
}

///macro NO_COLLISION_ZERO_CHECK turn off special check of case when t = 0 for initially flat configuration
bool GeomProjs::PerformEECCD(const Vector& a0, const Vector& b0, const Vector& c0, const Vector& d0,
                             const Vector& a1, const Vector& b1, const Vector& c1, const Vector& d1, DReal* t){
    CCDQuery ccd(a0, b0, c0, d0, a1, b1, c1, d1);
    ccd.should_find_time = (t != nullptr);
    ccd.PerformEETest();
    if (t != nullptr && ccd.has_collision) *t = ccd.t_res;
    return ccd.has_collision;
}







