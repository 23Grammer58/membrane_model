//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"

static int create_cilinder_part(AniMesh& am, double radius, double len, double phi, double size) {
    int nVVert = 4,  nLine = 4,  nSurface = 1;
    //double bb = len + radius;
    double r = radius, l = len;
    double VVert[] = {r, 0, -l/2,  r, 0, l/2,  r*cos(phi), r*sin(phi), l/2,  r*cos(phi), r*sin(phi), -l/2};
    //{-2*bb, -2*bb, -2*bb,  2*bb, 2*bb, 2*bb,  0, 0, 0,  bb, bb, -bb};  /* these will be recomputed, just define the bbox */

    int LineD[] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 1, 1}; /*edges*/
    int LineP[] = {1, 1,  1, 2,  1, 3,  1, 4}; /*param functions*/
    int exportCurves[] = {1, 2, 4, 2}; /*curve colors*/
    double LineT[] = {-l/2, l/2,  0, phi,  l/2, -l/2,  phi, 0}; /*parameter pairs*/

    int SurfL[] = {4 /*edges*/, 1 /*param*/, 1, 0 /*backcolor*/, 1 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/,  3, 0 /*direction*/,  4, 0};
    double SurfT[4] = {0,  phi, -l/2, l/2}; /*parametrization bbox*/

    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [phi, len](int i, double t, double *pu, double *pv) {
        switch(i) {
            case 1: { *pu = 0.0,  *pv = t;      break; }
            case 2: { *pu = t,    *pv = len/2;  break; }
            case 3: { *pu = phi,  *pv = t;      break; }
            case 4: { *pu = t,    *pv = -len/2; break; }
        }
        return 1;
    };
    asd.bounsurf = [R=radius](int i, double u, double v, double *px, double *py, double *pz) {
        return *px = R * cos(u), *py = R * sin(u), *pz = v, 1;
    };
    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_cilinder_part(double R0, double l, double phi, double size){
    AniMesh am;
    create_cilinder_part(am, R0, l, phi, size);
    mark_vertices(am);

    return am;
}
