//
// Created by alex on 01.04.2021.
//

#include "AVSim/Core/MeshGen/Helper.h"

static int create_sphere_eigther(AniMesh& am, double radius, double size) {
    int nVVert = 3,  nLine = 3,  nSurface = 1;
    double VVert[] = {-2*radius, -2*radius, -2*radius,  2*radius, 2*radius, 2*radius,  0, 0, 2*radius};  /* these will be recomputed, just define the bbox */
    int LineD[] = {2, 3, 1,  3, 1, 1,  1, 2, 1}; /*edges*/
    int LineP[] = {1, 1,  1, 2,  1, 3}; /*param functions*/
    int exportCurves[] = {1, 2, 3}; /*curve colors*/
    double LineT[] = {M_PI/2, 0,  M_PI/2, 0,  0, M_PI/2 }; /*parameter pairs*/
    int SurfL[] = {3 /*edges*/, 1 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/,  3, 0 /*direction*/};
    double SurfT[4] = {-1.0,  1.0, -1.0, 1.0}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [](int i, double t, double *pu, double *pv){
        switch(i){
            case 1: { *pu = 0.0   , *pv = sin(t); break; }
            case 2: { *pu = cos(t), *pv = 0     ; break; }
            case 3: { *pu = cos(t), *pv = sin(t); break; }
        }
        return 1;
    };
    asd.bounsurf = [R = radius](int i, double u, double v, double *px, double *py, double *pz) {
        double z = 1.0 - u*u - v*v;
        if (fabs(u) < DBL_EPSILON) u = 0;
        if (fabs(v) < DBL_EPSILON) v = 0;
        z = (z < 0.0 || fabs(z) < 10*DBL_EPSILON) ?  0.0 : R * sqrt(z);

        *px = R * u,  *py = R * v, *pz = ((i % 2) ? 1 : -1) * z;
        return 1;
    };
    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_eigth_sphere_part(double R0, double size){
    AniMesh am;
    create_sphere_eigther(am,  R0, size);
    mark_vertices(am);
    auto& vm = am.vertex_label; auto vertex = am.vertices.data();
    for (int i = 0; i < 3 && i < vm.size(); ++i){
        double err = size * 2e-7;
        int mask[4] = {0, 1, 2, 4};
        vm[i] = 0;
        for (int j = 0; j < 3; ++j)
            if (fabs(vertex[3*i + j]) < err) vm[i] |= mask[j+1];
    }
    for (int i = 3; i < vm.size(); ++i){
        int mask[4] = {0, 1, 2, 4};
        vm[i] = mask[vm[i]];
    }

    return am;
}

