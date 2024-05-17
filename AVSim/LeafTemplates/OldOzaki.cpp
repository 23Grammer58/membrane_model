//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"

static int create_old_ozaki(AniMesh& am, int ozaki_size, double size) {
    int nVVert = 5,  nLine = 5,  nSurface = 1;
    double VVert[5*3] = {-60, -60, -1,  60, 60, 1,  0,0,0, 0,0,0, 0,0,0};
    int LineD[5*3] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 5, 1,  5, 1, 1}; /*edges*/
    int LineP[5*2] = {1, 1,  1, 2,  1, 6,  1, 4,  1, 5}; /*param functions*/
    int exportCurves[5] = {1, 2, 2, 2, 4}; /*curve colors*/
    double LineT[5*2] = {0, 1,  0, 1,  0, 1, 0, 1,  0, 1}; /*parameter pairs*/
    int SurfL[5] = {5 /*edges*/, 1 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[5*2] = {1, 1 /*direction*/,  2, 1 /*direction*/,  3, 1,  4, 1,  5, 1};
    double SurfT[4] = {-60.0,  60.0, -60.0, 60.0}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};

    struct OldOzakiPrm{ double a, b, alpha, beta, radius, shift; };
    OldOzakiPrm ozaki;

    if (ozaki_size == 13) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 11.4 * M_PI / 180.0,  ozaki.a = 10.6,  ozaki.b = 11.0;
    } else if (ozaki_size == 15) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 10.9 * M_PI / 180.0,  ozaki.a = 11.5,  ozaki.b = 11.2;
//    } else if (ozaki_size == 101) { //my size
//        ozaki.alpha = 3.1 * M_PI / 180.0,  ozaki.beta = 10.9 * M_PI / 180.0,  ozaki.a = 13.1,  ozaki.b = 11;
    } else if (ozaki_size == 17) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  9.0 * M_PI / 180.0,  ozaki.a = 12.8,  ozaki.b = 13.3;
    } else if (ozaki_size == 19) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.8 * M_PI / 180.0,  ozaki.a = 13.8,  ozaki.b = 13.1;
    } else if (ozaki_size == 21) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  7.8 * M_PI / 180.0,  ozaki.a = 14.7,  ozaki.b = 13.6;
    } else if (ozaki_size == 23) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  7.7 * M_PI / 180.0,  ozaki.a = 15.8,  ozaki.b = 13.5;
    } else if (ozaki_size == 25) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 10.0 * M_PI / 180.0,  ozaki.a = 16.9,  ozaki.b = 13.5;
    } else if (ozaki_size == 27) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta = 10.0 * M_PI / 180.0,  ozaki.a = 17.9,  ozaki.b = 13.6;
    } else if (ozaki_size == 29) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  9.3 * M_PI / 180.0,  ozaki.a = 18.8,  ozaki.b = 13.7;
    } else if (ozaki_size == 31) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.3 * M_PI / 180.0,  ozaki.a = 19.9,  ozaki.b = 14.0;
    } else if (ozaki_size == 33) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.2 * M_PI / 180.0,  ozaki.a = 20.8,  ozaki.b = 14.1;
    } else if (ozaki_size == 35) {
        ozaki.alpha = 5.0 * M_PI / 180.0,  ozaki.beta =  8.0 * M_PI / 180.0,  ozaki.a = 21.8,  ozaki.b = 14.1;
    } else {
        throw std::runtime_error("unknown ozaki size");
    }
    ozaki.radius = (ozaki.a * cos(ozaki.beta) - ozaki.b*sin(ozaki.alpha)) / cos(ozaki.alpha);
    ozaki.shift = -ozaki.b*cos(ozaki.alpha) + ozaki.radius * sin(ozaki.alpha);

    asd.bounline = [ozaki](int i, double t, double *pu, double *pv) {
        double ax = 0, ay = ozaki.a * sin(ozaki.beta);
        double bx = ozaki.a * cos(ozaki.beta), by = 0;
        double cx = ozaki.a * cos(ozaki.beta) - ozaki.b*sin(ozaki.alpha), cy = -ozaki.b*cos(ozaki.alpha);
        double dx = -ozaki.a * cos(ozaki.beta) + ozaki.b*sin(ozaki.alpha), dy = -ozaki.b*cos(ozaki.alpha);
        double ex = -ozaki.a * cos(ozaki.beta), ey = 0;
        switch(i){
            case 1: { *pu = (1-t)*ax + t*bx,  *pv = (1-t)*ay + t*by; break; }
            case 2: { *pu = (1-t)*bx + t*cx,  *pv = (1-t)*by + t*cy; break; }
            case 3: { *pu = (1-t)*cx + t*dx,  *pv = (1-t)*cy + t*dy; break; }
            case 4: { *pu = (1-t)*dx + t*ex,  *pv = (1-t)*dy + t*ey; break; }
            case 5: { *pu = (1-t)*ex + t*ax,  *pv = (1-t)*ey + t*ay; break; }
            case 6: { *pu = ozaki.radius*cos(ozaki.alpha + t*(M_PI-2*ozaki.alpha));
                *pv = -ozaki.radius*sin(ozaki.alpha + t*(M_PI-2*ozaki.alpha)) + ozaki.shift;
                break; }
        }
        return 1;
    };
    asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1; };

    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_old_ozaki(int ozaki_size, double size){
    AniMesh am;
    create_old_ozaki(am, ozaki_size, size);
    mark_vertices(am);
    divideBorderElems(am);

    return am;
}

