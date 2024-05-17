//
// Created by alex on 01.04.2021.
//

#include "AVSim/Core/MeshGen/Helper.h"

static int create_dual_aligned_rectangle(AniMesh& am, double a, double b, double size){
    const int nVVert = 9;
    const int nLine = 12;
    const int nSurface = 4;
    double VVert[3*nVVert] = {-a/2,-b/2, 0,   0, -b/2, 0,   a/2, -b/2, 0,   a/2, 0, 0,  a/2, b/2, 0,  0, b/2, 0,  -a/2, b/2, 0,  -a/2, 0, 0,  0, 0, 0};
    int LineD[3 * nLine] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 5, 1,  5, 6, 1,  6, 7, 1,  7, 8, 1,  8, 1, 1,  2, 9, 1,  4, 9, 1,   6, 9, 1,  8, 9, 1}; /*edges*/
    int LineP[2 * nLine] = {0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0};
    double LineT[2 * nLine] = {0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0};
    int exportCurves[nLine] = {4, 4, 2, 2, 8, 8, 1, 1, 32, 16, 32, 16};
    int SurfL[nSurface * 5] =
            { 4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/,
              4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/,
              4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/,
              4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {8, 0,  1, 0,   9, 0,  12, 1,
                   2, 0,  3, 0,  10, 0,   9, 1,
                   4, 0,  5, 0,  11, 0,  10, 1,
                   6, 0,  7, 0,  12, 0,  11, 1 };
    double SurfT[4 * nSurface] =
            { 0, 0, 0, 0,
              0, 0, 0, 0,};
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_dual_aligned_rectangle(double a, double b, double size){
    if (a < 0.7*size || b < 0.7*size) {
        throw runtime_error("Impossible to generate rectangle with such params a = " + to_string(a) + ", b = " + to_string(b) + ", h = " + to_string(size) + "\n");
    }
    AniMesh am;
    create_dual_aligned_rectangle(am, a, b, size);
    mark_vertices(am);

    return am;
}

