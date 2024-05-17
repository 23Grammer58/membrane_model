//
// Created by alex on 01.04.2021.
//

#include "AVSim/Core/MeshGen/Helper.h"

static int create_half_splited_rectangle(AniMesh& am, double a, double b, double size){
    int nVVert = 6;
    int nLine = 7;
    int nSurface = 2;
    double VVert[] = {-a/2,-b/2, 0,   a/2, -b/2, 0,   a/2, 0, 0,  -a/2, 0, 0,  a/2, b/2, 0,   -a/2, b/2, 0};
    int LineD[] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 1, 1,  3, 5, 1,  5, 6, 1,  6, 4, 1}; /*edges*/
    int LineP[] = {0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0};
    double LineT[] = {0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0};
    int exportCurves[] = {1, 2, 4, 2, 2, 1, 2};
    int SurfL[] = { 4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/,
                    4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/,  3, 0 /*direction*/,  4, 0,
                   3, 1,  5, 0,  6, 0,  7, 0};
    double SurfT[] = {0, 0, 0, 0,
                      0, 0, 0, 0,}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_half_splited_rectangle(double a, double b, double size){
    if (a < 0.5*size || b < 0.5*size) {
        throw runtime_error("Impossible to generate rectangle with such params a = " + to_string(a) + ", b = " + to_string(b) + ", h = " + to_string(size) + "\n");
    }
    AniMesh am;
    create_half_splited_rectangle(am, a, b, size);
    mark_vertices(am);

    return am;
}
