//
// Created by alex on 01.04.2021.
//

#include "Elliptic.h"
#include <stdexcept>

gsl_root_fdfsolver * gS;
gsl_function_fdf gFDF;
double gE_k;
Ellint_E_params* gP;

static double ellint_E_f(double x, void* params)
{
    return gsl_sf_ellint_E(x, ((Ellint_E_params*)params)->k, GSL_PREC_DOUBLE) - ((Ellint_E_params*)params)->u;
}
static double ellint_E_df(double x, void* params)
{
    double k = ((Ellint_E_params*)params)->k;
    double sinx = sin(x);
    return sqrt(1 - k * k * sinx * sinx);
}

static void ellint_E_fdf(double x, void* params, double * f, double * df)
{
    *f = ellint_E_f(x, params);
    *df = ellint_E_df(x, params);
}

static double inv_ellint_E_guess(double u, double k)
{
    return M_PI_2 * u / gE_k;
}

void inv_ellint_E_init(double k)
{
    Ellint_E_params* p = (Ellint_E_params*)malloc(sizeof(Ellint_E_params)*1);
    p->k = k;
    p->u = 0;
    gP = p;

    const gsl_root_fdfsolver_type * T = gsl_root_fdfsolver_newton;
    gS = gsl_root_fdfsolver_alloc (T);
    gFDF.f = ellint_E_f;
    gFDF.df = ellint_E_df;
    gFDF.fdf = ellint_E_fdf;
    gFDF.params = p;
    gE_k = gsl_sf_ellint_Ecomp(p->k, GSL_PREC_DOUBLE);
}

void inv_ellint_E_free()
{
    gsl_root_fdfsolver_free (gS);
    free(gP);
}

double inv_ellint_E(double u, double k)
{
    gP->u = u;
    gP->k = k;
    double r_guess = inv_ellint_E_guess(u, k);
    double x0 = 0, x = r_guess;
    double err = 1e-6;
    gsl_root_fdfsolver_set(gS, &gFDF, x);

    int status;
    int iter = 0, max_iter = 120;
    do
    {
        iter++;
        status = gsl_root_fdfsolver_iterate(gS);
        x0 = x;
        x = gsl_root_fdfsolver_root (gS);
        status = gsl_root_test_delta (x, x0, 0, err);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    return x;
}

/*
    solve following system:
        nx * x + ny * y = d
        (x / a)^2 + (y / b)^2 = 1
        condition
    sign should be in {-2, -1, 1, 2}
    if (sign % 2) == 0 -> crd = x
    else crd = y
    condition: choose higher ([sign / 2] > 0) or lower ([sign / 2] < 0) root according crd coordinate
    return  (-1) - there is no roots
            (-2) - wrong ellips parameters
            (-3) - wrong sign
            (-4) - there is infty roots (may mean error)
            (0) - success
*/
int solve_2d_line_to_ellips(double nx, double ny, double d, double a, double b, int sign, double res[2])
{
    if (fabs(nx) < DBL_EPSILON  && fabs(ny) < DBL_EPSILON)
    {
        if (fabs(d) < DBL_EPSILON)
        {
            res[0] = a;
            res[1] = 0;
            return -4;
        }
        else return -1;
    }
    if (a < DBL_EPSILON || b < DBL_EPSILON) return -2;
    if (sign == 0 || abs(sign) > 2) return -3;
    if (fabs(nx) < DBL_EPSILON)
    {
        int map[5] = {-1, -2, 0, 2, 1};
        int newsign = map[(sign + 2)];
        int err = solve_2d_line_to_ellips(ny, nx, d, b, a, newsign, res);
        double tmp = res[0];
        res[0] = res[1]; res[1] = tmp;
        return err;
    }

    double al = d / nx / a, bt = -ny / nx / a;
    double a1 = bt * bt + 1.0 / (b * b), k1 = al * bt, c1 = al * al - 1;
    double p = k1 / a1, q = c1 / a1;
    double D = p * p - q;
    if (D < 0) return -1;
    double y = -p, dy = sqrt(D);
    int sg = 0;
    if (!(sign % 2)) sg = sign / 2;
    else sg = (((- ny * dy / nx) > 0) ? 1 : -1) * sign;
    y = y + sg * dy;
    double x = d / nx - ny * y / nx;
    res[0] = x, res[1] = y;

    return 0;
}

static void _swap(double* x, double* y)
{
    double t = *x;
    *x = *y;
    *y = t;
}

// return eigen vectors h1 and h2: ||h1^T A h1|| = ||h2^T A h2|| = 1
static int _get_eigen_params_2d_sym_matrix(/*matrix*/ double diag[2], double nodiag, /*out*/ double h[2][2])
{
    double p = (diag[0] + diag[1])/2, q = diag[0] * diag[1] - nodiag * nodiag;
    double D = sqrt(p * p - q);
    double l[2] = {p + D, p - D};
    int event = 0;
    if (fabs(nodiag) > DBL_EPSILON)
        for (int i = 0; i < 2; ++i)
        {
            double r = sqrt(fabs(l[i]) * (nodiag * nodiag + (l[i] - diag[0]) * (l[i] - diag[0])) );
            h[0][i] = nodiag / r;
            h[1][i] = (l[i] - diag[0])/ r;
        }
    else
    {
        l[0] = diag[0], l[1] = diag[1];
        if (fabs(l[0]) > DBL_EPSILON)
            h[0][0] = 1 / sqrt(fabs(l[0])), h[1][0] = 0, event = 2;
        else
            h[0][0] = 1, h[1][0] = 0;
        if (fabs(l[1]) > DBL_EPSILON)
            h[0][1] = 0, h[1][1] = 1 / sqrt(fabs(l[1])), event = 4;
        else
            h[0][0] = 0, h[1][0] = 1;

        if (fabs(l[0]) < fabs(l[1]))
            _swap(&h[0][0], &h[0][1]), _swap(&h[1][0], &h[1][1]);
    }


    return (l[1] < 0 || l[0] < 0 ) | event;
}

static int _cross_product(double p1[3], double p2[3], double res[3])
{
    for (int i = 0; i < 3; ++i)
        res[i] = p1[(i + 1) % 3] * p2[(i + 2) % 3] - p2[(i + 1) % 3] * p1[(i + 2) % 3];
    return 0;
}

static int _norm(double p[3], double res[3])
{
    double len = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    for (int i = 0; i < 3; ++i)
        res[i] = p[i]/len;
    return 0;
}

static int _scal(double p[3], double coef, double res[3])
{
    for (int i = 0; i < 3; ++i)
        res[i] = p[i] * coef;
    return 0;
}

static int _get_max_index(double p[3])
{
    double s[3] = {fabs(p[0]), fabs(p[1]), fabs(p[2])};
    if (s[0] > s[1])
    {
        if (s[0] > s[2]) return 0;
        else return 2;
    }
    else
    {
        if (s[1] > s[2]) return 1;
        else return 2;
    }
}

static double* _get_nonparallel_vector(double n[3], double res[3])
{
    int id = _get_max_index(n);
    for (int i = 0; i < 3; ++i)
        res[i] = 0;
    res[(id + 1)%3] = 1;
    return res;
}

/* find intersection between cylinder (1) and plane (2)
    (1) (x/a)^2 + (y/b)^2 = 1
    (2) (n, r) = 0
    return points of long and small axis intersection of plane and cylinder with positive z component
*/
int solve_plane_to_cylinder_ellips(double a, double b, double n[3], double pa[3], double pb[3])
{
    double f[2][3], mem[3];
    _cross_product(n, _get_nonparallel_vector(n, mem), &(f[0][0]));
    _cross_product(n, &(f[0][0]), &(f[1][0]));
    _norm(&(f[0][0]), &(f[0][0]));
    _norm(&(f[1][0]), &(f[1][0]));
    double diag[2] = {(f[0][0]/a) * (f[0][0]/a) + (f[0][1]/b) * (f[0][1]/b), (f[1][0]/a) * (f[1][0]/a) + (f[1][1]/b) * (f[1][1]/b)};
    double nodiag = (f[0][0]/a) * (f[1][0]/a) + (f[0][1]/b) * (f[1][1]/b);
    double c[2][2];
    _get_eigen_params_2d_sym_matrix(diag, nodiag, c);
    for (int i = 0; i < 3; ++i)
    {
        pa[i] = c[0][1] * f[0][i] + c[1][1] * f[1][i];
        pb[i] = c[0][0] * f[0][i] + c[1][0] * f[1][i];
    }
    if (pa[2] < 0) _scal(pa, -1, pa);
    if (pb[2] < 0) _scal(pb, -1, pb);

    return 0;
}

InvEllInt::InvEllInt(double k) {
    m_S = std::make_shared<M_gsl_root_fdfsolver>();
    m_FDF.f = ellint_E_f;
    m_FDF.df = ellint_E_df;
    m_FDF.fdf = ellint_E_fdf;
    m_k = k;
    if (m_k >= 0) m_E_k = gsl_sf_ellint_Ecomp(m_k, GSL_PREC_DOUBLE);
}

void InvEllInt::setK(double k) {
    if (k < 0) throw std::runtime_error("k value of inverse elliptic integral must be nonnegative");
    m_k = k;
    m_E_k = gsl_sf_ellint_Ecomp(m_k, GSL_PREC_DOUBLE);
}

double InvEllInt::compute(double u, double k) const {
    Ellint_E_params prms; prms.u = u; prms.k = k;
    m_FDF.params = &prms;
    double E_k = (m_k == k) ? m_E_k : gsl_sf_ellint_Ecomp(prms.k, GSL_PREC_DOUBLE);
    if (prms.k < 0) throw std::runtime_error("k value of inverse elliptic integral must be nonnegative");

    double r_guess = M_PI_2 * u / E_k;
    double x0 = 0, x = r_guess;
    double err = 1e-6;
    gsl_root_fdfsolver_set(m_S->s, &m_FDF, x);

    int status;
    int iter = 0, max_iter = 120;
    do
    {
        iter++;
        status = gsl_root_fdfsolver_iterate(m_S->s);
        x0 = x;
        x = gsl_root_fdfsolver_root (m_S->s);
        status = gsl_root_test_delta (x, x0, 0, err);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    return x;
}

