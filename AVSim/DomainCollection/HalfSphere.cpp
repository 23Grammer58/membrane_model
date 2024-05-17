//
// Created by alex on 01.04.2021.
//

#include "AVSim/Core/MeshGen/Helper.h"

static int create_half_sphere(AniMesh& am, double R, double size){
    int nVVert = 2,  nLine = 2,  nSurface = 1;
    double VVert[] = {-2*R, -2*R, -2*R, 2*R, 2*R, 2*R};  /* these will be recomputed, just define the bbox */
    int LineD[] = {1, 2, 1,  1, 2, 1,  1, 2, 1}; /*edges*/
    int LineP[] = {1, 3,  1, 4,  1, 3}; /*param functions*/
    int exportCurves[] = {1, 2, 4}; /*curve colors*/
    double LineT[] = {0, M_PI,  1, -1, 0, -M_PI }; /*parameter pairs*/
    int SurfL[] = {2, 1, 1, 0, 0,
                   2, 1, 1, 0, 0}; /*main surface*/
    int SurfI[] = {1, 0,  2, 1,  2, 0, 3, 1};
    double SurfT[4] = {-1.1, 1.1, -1.1, 1.1}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [](int i, double t, double *pu, double *pv){
        switch(i){
            case 1: { *pu = 0.0   , *pv = sin(t); break; }
            case 2: { *pu = cos(t), *pv = 0     ; break; }
            case 3: { *pu = cos(t), *pv = sin(t); break; }
            case 4: { *pu = t,      *pv = 0;      break; }
        }
        return 1;
    };
    asd.bounsurf = [R](int i, double u, double v, double *px, double *py, double *pz) {
        double z = 1.0 - u*u - v*v;
        if (fabs(u) < DBL_EPSILON) u = 0;
        if (fabs(v) < DBL_EPSILON) v = 0;
        z = (z < 0.0 || fabs(z) < 10*DBL_EPSILON) ?  0.0 : sqrt(z);
        z *= ((i % 2) ? 1 : -1);

        *px = R * u,  *py = R * v, *pz = R * z;
        return 1;
    };
    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_half_sphere(double R, double size){
    AniMesh am;
    create_half_sphere(am, R, size);
    mark_vertices(am);

    return am;
}

