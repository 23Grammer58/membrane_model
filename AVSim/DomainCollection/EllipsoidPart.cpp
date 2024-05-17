//
// Created by alex on 18.11.2021.
//

#include "AVSim/Core/MeshGen/Helper.h"

static int create_EllipsoidPart0(AniMesh& am, double a, double b, double c, double theta, double size){
    double R = max({a, b, c});
    double sclu = 1.0/a, sclv = 1.0/c;
    auto line_func = [theta, sclu, sclv](int i, double t, double *pu, double *pv){
        switch(i){
            case 1: {
                if (t <= 0 || t >= 2*M_PI) {
                    *pu = sin(theta), *pv = cos(theta);
                } else if (t == M_PI){
                    *pu = -sin(theta), *pv = cos(theta);
                } else {
                    *pu = sin(theta)*cos(t); *pv = cos(theta);
                }
                break;
            }
            case 2: {
                if (t <= theta){
                    *pu = sin(theta), *pv = cos(theta);
                } else if (t >= 2*M_PI - theta){
                    *pu = -sin(theta), *pv = cos(theta);
                } else if (t == M_PI){
                    *pu = 0, *pv = -1;
                } else if (t < M_PI){
                    *pu = sin(t), *pv = cos(t);
                } else if (t > M_PI){
                    *pu = sin(t), *pv = cos(2*M_PI - t);
                }

                break;
            }
            default: {
                std::cout << "Something wrong" <<std::endl;
            }
        }
        *pu *= sclu, *pv *= sclv;

        return 0;
    };
    auto surf_func = [a,b,c, sclu, sclv](int i, double u, double v, double *px, double *py, double *pz) {
        double x = u / sclu;
        double z = v / sclv;
        double y = 1 - (x*x + z*z);
        y  = (y <= 0) ? 0 : sqrt(y);
        switch (i) {
            case 1: {
                y *= -1;
                break;
            }
            case 2: {
                break;
            }
            default: {
                std::cout << "Something wrong" <<std::endl;
            }
        }
        *px = a * x;
        *py = b * y;
        *pz = c * z;

        return 0;
    };
    int nVVert = 3,  nLine = 4,  nSurface = 2;
    double VVert[] = {-2*R, -2*R, -2*R, 2*R, 2*R, 2*R, 0, 0, -c};  /* these will be recomputed, just define the bbox */
    int LineD[] = {1, 2, 1,  1, 3, 2, 3, 2, 2,  1, 2, 1}; /*edges (p_st, p_end, count common domains)*/
    int LineP[] = {1, 1,     //1: (number of subdomain, label)
                   1, 2,  2, 2, //2: (number of subdomain, label)
                   1, 2,  2, 2, //3: (number of subdomain, label)
                   2, 1}; //4: (number of subdomain, label)
    int exportCurves[] = {2, 1, 1, 2}; /*curve colors*/
    double LineT[] = {0, M_PI,
                      theta, M_PI, theta, M_PI, M_PI,
                      2*M_PI - theta, M_PI, 2*M_PI - theta,
                      0, M_PI}; /*parameter pairs*/
    int SurfL[] = {3/*nlines*/, 1/*subdom number*/, 1/*count volume subdom*/, 0, 0,
                   3, 2, 1, 0, 0}; /*main surface*/
    int SurfI[] = {1, 0, 3, 1, 2, 1,
                   4, 0, 3, 1, 2, 1};
    double SurfT[] = {sclu*(-1-DBL_EPSILON), sclu*(1+DBL_EPSILON), sclv*(-1-DBL_EPSILON), sclv*(cos(theta)+DBL_EPSILON),
                      sclu*(-1-DBL_EPSILON), sclu*(1+DBL_EPSILON), sclv*(-1-DBL_EPSILON), sclv*(cos(theta)+DBL_EPSILON),}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = line_func;
    asd.bounsurf = surf_func;
    int stat = makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));

    return stat;
}

static int create_EllipsoidPart1(AniMesh& am, double a, double b, double c, double theta, double size){
    double R = max({a, b, c});
    auto line_func = [theta](int i, double t, double *pu, double *pv){
        switch(i){
            case 1: { *pu = t; *pv = cos(theta); break; }
            case 2: { *pu = sin(t), *pv = (t <= M_PI) ? cos(t) : cos(2*M_PI - t); break; }
            case 3: { *pu = t; *pv = 0; break; }
            default: {
                std::cout << "Something wrong" <<std::endl;
            }
        }

        return 0;
    };
    auto surf_func = [a,b,c](int i, double u, double v, double *px, double *py, double *pz) {
        double x = u;
        double z = v;
        double y = 1 - (x*x + z*z);
        y  = (y <= 0) ? 0 : sqrt(y);
        switch (i) {
            case 1: {
                y *= -1;
                break;
            }
            case 2: {
                break;
            }
            default: {
                std::cout << "Something wrong" <<std::endl;
            }
        }
        *px = a * x;
        *py = b * y;
        *pz = c * z;

        return 0;
    };
    int nVVert = 5,  nLine = 8,  nSurface = 4;
    double VVert[] = {-2*R, -2*R, -2*R, 2*R, 2*R, 2*R, 0, 0, -c, a, 0, 0, -a, 0, 0};  /* these will be recomputed, just define the bbox */
    int LineD[] = {1, 2, 1,  1, 4, 2,  4, 3, 2,  3, 5, 2,  5, 2, 2,  1, 2, 1,  4, 5, 1, 4, 5, 1}; /*edges (p_st, p_end, count common domains)*/
    int LineP[] = {1, 1,     //1: (number of subdomain, label)
                   1, 2,  2, 2, //2: (number of subdomain, label)
                   1, 2,  2, 2, //3: (number of subdomain, label)
                   1, 2,  2, 2, //4: (number of subdomain, label)
                   1, 2,  2, 2, //5: (number of subdomain, label)
                   2, 1, //6: (number of subdomain, label)
                   1, 3, //7: (number of subdomain, label)
                   2, 3}; //8: (number of subdomain, label)
    int exportCurves[] = {2, 1, 1, 1, 1, 2, 1, 1}; /*curve colors*/
    double LineT[] = {sin(theta), -sin(theta),
                      theta, M_PI_2, theta, M_PI_2,
                      M_PI_2, M_PI, M_PI_2, M_PI,
                      M_PI, 3*M_PI_2, M_PI, 3*M_PI_2,
                      3*M_PI_2, 2*M_PI - theta, 3*M_PI_2, 2*M_PI - theta,
                      sin(theta), -sin(theta),
                      1, -1,
                      1, -1}; /*parameter pairs*/
    int SurfL[] = {4/*nlines*/, 1/*subdom number*/, 1/*count volume subdom*/, 0, 0,
                   3, 1, 1, 0, 0,
                   4, 2, 1, 0, 0,
                   3, 2, 1, 0, 0}; /*main surface*/
    int SurfI[] = {1, 0, 5, 1, 7, 1, 2, 1,
                   7, 0, 4, 1, 3, 1,
                   6, 0, 5, 1, 8, 1, 2, 1,
                   8, 0, 4, 1, 3, 1};
    double SurfT[] = {-1, 1, (-1-DBL_EPSILON), (cos(theta)+DBL_EPSILON),
                      -1, 1, (-1-DBL_EPSILON), (cos(theta)+DBL_EPSILON),
                      -1, 1, (-1-DBL_EPSILON), (cos(theta)+DBL_EPSILON),
                      -1, 1, (-1-DBL_EPSILON), (cos(theta)+DBL_EPSILON),}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = line_func;
    asd.bounsurf = surf_func;
    int stat = makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));

    return stat;
}

static int create_EllipsoidPart2(AniMesh& am, double a, double b, double c, double theta, double size){
    double R = max({a, b, c});
    auto line_func = [theta, sin_theta = sin(theta)](int i, double t, double *pu, double *pv){
        switch(i){
            case 1: { *pu = sin_theta*cos(t); *pv = sin_theta*sin(t); break; }
            case 2: { *pu = t, *pv = 0; break; }
            case 3: { *pu = cos(t); *pv = sin(t); break; }
            default: {
                std::cout << "Something wrong" <<std::endl;
            }
        }

        return 0;
    };
    auto surf_func = [a,b,c](int i, double u, double v, double *px, double *py, double *pz) {
        double x = u;
        double y = v;
        double z = 1 - (x*x + y*y);
        z  = (z <= 0) ? 0 : sqrt(z);
        switch (i) {
            case 1: {
                z *= -1;
                break;
            }
            case 2: {
                break;
            }
            default: {
                std::cout << "Something wrong" <<std::endl;
            }
        }
        *px = a * x;
        *py = b * y;
        *pz = c * z;

        return 0;
    };
    int nVVert = 5,  nLine = 8,  nSurface = 4;
    double VVert[] = {-2*R, -2*R, -2*R, 2*R, 2*R, 2*R, 0, 0, -c, a, 0, 0, -a, 0, 0};  /* these will be recomputed, just define the bbox */
    int LineD[] = {1, 2, 1,  1, 4, 1,  4, 3, 1,  3, 5, 1,  5, 2, 1,  1, 2, 1,  4, 5, 2, 4, 5, 2}; /*edges (p_st, p_end, count common domains)*/
    int LineP[] = {1, 1,        //1: (number of subdomain, label)
                   1, 2,        //2: (number of subdomain, label)
                   2, 2,        //3: (number of subdomain, label)
                   2, 2,        //4: (number of subdomain, label)
                   1, 2,        //5: (number of subdomain, label)
                   1, 1,        //6: (number of subdomain, label)
                   1, 3, 2, 3,  //7: (number of subdomain, label)
                   1, 3, 2, 3}; //8: (number of subdomain, label)
    int exportCurves[] = {2, 1, 1, 1, 1, 2, 1, 1}; /*curve colors*/
    double LineT[] = {0, M_PI,
                      sin(theta), 1,
                      1, 0,
                      0, -1,
                      -1, -sin(theta),
                      0, -M_PI,
                      0, M_PI, 0, M_PI,
                      0, -M_PI, 0, -M_PI,}; /*parameter pairs*/
    int SurfL[] = {4/*nlines*/, 1/*subdom number*/, 1/*count volume subdom*/, 0, 0,
                   3, 2, 1, 0, 0,
                   4, 1, 1, 0, 0,
                   3, 2, 1, 0, 0}; /*main surface*/
    int SurfI[] = {2, 0, 7, 0, 5, 0, 1, 1,
                   7, 0, 4, 1, 3, 1,
                   6, 0, 5, 1, 8, 1, 2, 1,
                   3, 0, 4, 0, 8, 1};
    double SurfT[] = {-1, 1, -1, 1,
                      -1, 1, -1, 1,
                      -1, 1, -1, 1,
                      -1, 1, -1, 1,}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = line_func;
    asd.bounsurf = surf_func;
    int stat = makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
    invertFaceOrientation(am);

    return stat;
}

AniMesh generate_EllipsoidPart(double a, double b, double d, double size){
    if (a <= 0 || b <= 0 || d <= 0 || d > a) throw std::runtime_error("Wrong input parameters");
    AniMesh am;
    if (a > b)
        create_EllipsoidPart2(am, a, a, b, asin(d / a), size);
    else
        create_EllipsoidPart1(am, a, a, b, asin(d / a), size);
    mark_vertices(am);
    return am;
}

