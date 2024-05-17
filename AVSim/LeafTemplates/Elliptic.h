//
// Created by alex on 01.04.2021.
//

#ifndef AORTIC_VALVE_ELLIPTIC_H
#define AORTIC_VALVE_ELLIPTIC_H

#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <memory>

typedef struct { double k; double u; } Ellint_E_params;

void inv_ellint_E_init(double k);

void inv_ellint_E_free();

double inv_ellint_E(double u, double k);

struct InvEllInt{
    struct M_gsl_root_fdfsolver{
        gsl_root_fdfsolver * s = nullptr;
        M_gsl_root_fdfsolver() {
            const gsl_root_fdfsolver_type * T = gsl_root_fdfsolver_newton;
            s = gsl_root_fdfsolver_alloc (T);
        }
        ~M_gsl_root_fdfsolver(){ if (s) gsl_root_fdfsolver_free (s); }
    };
    std::shared_ptr<M_gsl_root_fdfsolver> m_S;
    mutable gsl_function_fdf m_FDF;
    double m_E_k, m_k;

    explicit InvEllInt(double k = 0);
    void setK(double k);
    double compute(double u) const { return compute(u, m_k); }
    double compute(double u, double k) const;
};

typedef struct
{
    double a;
    double b;
} elliptic_prms;

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
int solve_2d_line_to_ellips(double nx, double ny, double d, double a, double b, int sign, double res[2]);

/* find intersection between elleptic cylinder (1) and plane (2)
    (1) (x/a)^2 + (y/b)^2 = 1
    (2) (n, r) = 0
    return points of long and small axis intersection of plane and cylinder with nonegative z component
*/
int solve_plane_to_cylinder_ellips(double a, double b, double n[3], double pa[3], double pb[3]);


#endif //AORTIC_VALVE_ELLIPTIC_H
