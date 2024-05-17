//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"
static int create_cutted_annulus(AniMesh& am, double R1, double R2, double size){
    int nVVert = 6;
    int nLine = 7;
    int nSurface = 2;
    double VVert[] = {R1, 0, 0,  R2, 0, 0,  -R1, 0, 0,  -R2, 0, 0,   R1, 0, 0,  R2, 0, 0};
    int LineD[] = {1, 2, 1,  2, 4, 1,  4, 3, 1,  3, 1, 1,  4, 6, 1,  6, 5, 1,  5, 3, 1};
    int LineP[] = {0, 0,     1, 3,     0, 0,     1, 1,     1, 3,     0, 0,     1, 1};
    int exportCurves[] = {4, 1, 1, 1, 1, 2, 1};
    double LineT[] = {0, 0,  2*M_PI, M_PI,   0, 0,  M_PI, 2*M_PI,  M_PI, 0,  0, 0,  0, M_PI};
    int SurfL[] = {4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0, /*direction*/
                   4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0,  2, 0,  3, 0,  4, 0,//{1, 0,  2, 0,  7, 0,  6, 0,
                   3, 1, 5, 0, 6, 0, 7, 0};
    double SurfT[] = {-R2, -R2, R2, R2}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [R1, R2](int i, double t, double *pu, double *pv){
        switch (i) {
            case 1: { *pu = R1 * cos(t),        *pv = R1 * sin(t); break; }
            case 2: { *pu = R1 + (R2 - R1) * t, *pv = 0;           break; }
            case 3: { *pu = R2 * cos(t),        *pv = R2 * sin(t); break; }
            default:
                cout << "Something wrong in annulus generation, i = " << i << " t = " << t << endl;
        }
        return i;
    };
    asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1;};

    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_cutted_annulus(double R1, double R2, double size){
    if (R1 > R2) std::swap(R1, R2);
    AniMesh am;
    create_cutted_annulus(am, R1, R2, size);
    mark_vertices(am);

    return am;
}

