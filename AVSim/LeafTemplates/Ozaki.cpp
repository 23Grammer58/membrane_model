//
// Created by alex on 31.08.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"
#include "TemplatesCollection.h"

static int create_ozaki(AniMesh& am, OzakiTempl p, double size) {
    int nVVert = 4,  nLine = 4,  nSurface = 1;
    double VVert[4*3] = {-60, -60, -1,  60, 60, 1,  0,0,0, 0,0,0};
    int LineD[4*3] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 1, 1}; /*edges*/
    int LineP[4*2] = {1, 1,  1, 2,  1, 2,  1, 5}; /*param functions*/
    int exportCurves[4] = {1, 2, 2, 4}; /*curve colors*/
    double LineT[4*2] = {0, 1,  0, 1.5,  1.5, 3,  0, 1}; /*parameter pairs*/
    int SurfL[5] = {4 /*edges*/, 1 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[4*2] = {1, 1 /*direction*/,  2, 1 /*direction*/,  3, 1,  4, 1};
    double SurfT[4] = {-60.0,  60.0, -60.0, 60.0}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};

    asd.bounline = [ozaki = p](int i, double t, double *pu, double *pv) {
    /* Scheme:
     *   A
     * E   B
     *  D C
     */
        double ax, ay, bx, by, cx, cy, dx, dy, ex, ey, r, x, y, gamma;
        r = ozaki.D/2.0 + ozaki.w;
        ax = 0.0, ay = ozaki.h + ozaki.beta;
        bx = r + ozaki.alpha, by = ozaki.h;
        if (ozaki.m > 0.0) {
            x = ozaki.m*(bx-ax)/sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
            y = ozaki.m*(by-ay)/sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
            bx += x,  by += y;
        }

        gamma = atan2(bx, -by) - acos(r/sqrt(bx*bx+by*by));

        if (ozaki.s > 0.0) {
            x = ozaki.s*(bx-ax)/sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
            y = ozaki.s*(by-ay)/sqrt((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
            bx -= x,  by -= y;
            r -= ozaki.s;
        }
        ex = -bx, ey = by;
        cx = r*sin(gamma),  cy = -r*cos(gamma);
        dx = r*sin(-gamma),  dy = -r*cos(-gamma);
        while (t > 1.0) {
            t -= 1.0;
            i = i % 5 + 1;
        }
        if (i == 1) {
            *pu = (1-t)*ax + t*bx,  *pv = (1-t)*ay + t*by;
        } else if (i == 2) {
            *pu = (1-t)*bx + t*cx,  *pv = (1-t)*by + t*cy;
        } else if (i == 3) {
            *pu = r*sin(gamma-2*t*gamma),  *pv = -r*cos(gamma-2*t*gamma);
        } else if (i == 4) {
            *pu = (1-t)*dx + t*ex,  *pv = (1-t)*dy + t*ey;
        } else if (i == 5) {
            *pu = (1-t)*ex + t*ax,  *pv = (1-t)*ey + t*ay;
        }
        return 1;
    };
    asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1; };

    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_ozaki(OzakiTempl p, double size){
    AniMesh am;
    create_ozaki(am, p, size);
    mark_vertices(am);
    divideBorderElems(am);

    return am;
}

