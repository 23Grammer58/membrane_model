//
// Created by alex on 27.04.2021.
//

#include "AVSim/Core/MeshGen/Helper.h"

static int create_HalfEllipseWithRectangle(AniMesh& am, double H1, double H2, double l, double size){
    if (H2 < size / 2) {
        int nVVert = 2;
        int nLine = 2;
        int nSurface = 1;
        double VVert[] = {-l / 2, 0, 0,  l / 2, 0, 0 };
        int LineD[] = {1, 2, 1,   2, 1, 1};
        int LineP[] = {1, 1, 0, 0 };
        int exportCurves[] = {4, 1};
        double LineT[] = {M_PI, 2 * M_PI, 0, 0};
        int SurfL[] = {2 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
        int SurfI[] = {1, 0 /*direction*/, 2, 0};
        double SurfT[] = {-2*l, 2*l, -2*H1, 2*H1};
        Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
        asd.bounline = [a=l/2, b = H1](int i, double t, double *pu, double *pv) { return *pu = a * cos(t), *pv = b * sin(t), 1; };
        asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz) { return *px = u, *py = v, *pz = 0, 1; };

        return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
    } else {
        int nVVert = 6;
        int nLine = 4;
        int nSurface = 1;
        double VVert[] = { -l, -2*H1, 0,  l, 2*H2, 0,  -l / 2, 0, 0,  l / 2, 0, 0,  l / 2, H2, 0,   -l / 2, H2, 0};
        int LineD[] = {3, 4, 1,  4, 5, 1,  5, 6, 1,  6, 3, 1};
        int LineP[] = { 1, 1,  0, 0,  0, 0,  0, 0 };
        int exportCurves[] = {2, 2, 1, 2};
        double LineT[] = {M_PI, 2 * M_PI, 0, 0,  0, 0,  0, 0};
        int SurfL[] = {4 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
        int SurfI[] = {1, 0,  2,0,  3,0,  4,0};
        double SurfT[] = {-2*l, 2*l, -2*(H1+H2), 2*(H1+H2)};
        Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
        asd.bounline = [a=l/2, b = H1](int i, double t, double *pu, double *pv) { return *pu = a * cos(t), *pv = b * sin(t), 1; };
        asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz) { return *px = u, *py = v, *pz = 0, 1; };

        int status = makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
        for (int fi = 0; fi < am.faces.size(); ++fi) am.faces[fi] -= 2;
        for (int ei = 0; ei < am.edges.size(); ++ei) am.edges[ei] -= 2;
        for (int vi = 6; vi < am.vertices.size(); ++vi) am.vertices[vi-6] = am.vertices[vi];
        am.vertices.resize(am.vertices.size()-6);

        return status;
    }
}

AniMesh generate_HalfEllipseWithRectangle(double H1, double H2, double l, double size){
    AniMesh am;
    create_HalfEllipseWithRectangle(am, H1, H2, l, size);
    mark_vertices(am);
    invertFaceOrientation(am);

    return am;
}

