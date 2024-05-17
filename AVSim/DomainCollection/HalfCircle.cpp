//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"

static int create_half_circle(AniMesh& am, double R, double size){
    int nVVert = 3;
    int nLine = 3;
    int nSurface = 1;
    double VVert[] = {-R/2, 0, 0,  R/2, 0, 0,  0, 0, 0};
    int LineD[] = {1, 2, 1,  2, 3, 1,  3, 1, 1};
    int LineP[] = {1, 1,  0, 0,  0, 0};
    int exportCurves[] = {4, 1, 1};
    double LineT[] = {M_PI, 2 * M_PI,  0, 0, 0, 0};
    int SurfL[] = {3 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/, 3, 0};
    double SurfT[] = {-R, 0, -R, R};
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [R](int i, double t, double *pu, double *pv){ return *pu = R * cos(t), *pv = R * sin(t), 1;};
    asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1;};

    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_half_circle(double R, double size){
    AniMesh am;
    create_half_circle(am, R, size);
    mark_vertices(am);

    return am;
}
