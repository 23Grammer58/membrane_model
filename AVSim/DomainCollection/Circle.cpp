//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"

static int create_circle(AniMesh& am, double R, double size){
    int nVVert = 3;
    int nLine = 4;
    int nSurface = 2;
    double VVert[] = {-R/2, 0, 0,  R/2, 0, 0,  0, 0, 0};
    int LineD[] = {1, 2, 1,  2, 3, 1,  3, 1, 1,  2, 1, 1};
    int LineP[] = {1, 1,  0, 0,  0, 0,  1, 3};
    int exportCurves[] = {1, 2, 4, 1};
    double LineT[] = {M_PI, 2 * M_PI,  0, 0, 0, 0,  0, M_PI};
    int SurfL[] = {3 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/,
                   3 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/, 3, 0,
                   3, 1 /*direction*/,  2, 1,  4, 0};
    double SurfT[] = {-R, 0, R, R,
                      -R, -R, 0, R}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [R](int i, double t, double *pu, double *pv){ return *pu = R * cos(t), *pv = R * sin(t), 1;};
    asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1;};

    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_circle(double R, double size){
    AniMesh am;
    create_circle(am, R, size);
    mark_vertices(am);

    return am;
}

// theta - parameter of coaxial circle radius: R_c = theta * R
// K - dimensional curvature of spherical surface, K_real = K / R
AniMesh generate_circle_with_coaxial(double R, double theta, double size, double K){
    assert(theta*R >= 2*size && theta*R <= R - 0.85*size && "Not allowed coaxial parameter");
    assert(K >= 0 && K <= 1 && "Wrong curvature");

    AniMesh am;

    const int nVVert = 4,  nLine = 5,  nSurface = 2;
    double VVert[nVVert*3] = {-2*R, -2*R, -2*R, 2*R, 2*R, 2*R, 0.0};
    int LineD[nLine*3] = {2, 1, 1,   1, 4, 1,  2, 3, 1,  3, 2, 1,       4, 1, 1     };
    int LineP[nLine*2] = {1, 1,      1, 2,     1, 3,     1, 3,          1, 2        }; /*param functions*/
    int exportCurves[] = {2,         1,        4,        4,             1           };
    double LineT[] =     {theta, 1,  0, M_PI,  0, M_PI,  M_PI, 2*M_PI,  M_PI, 2*M_PI};
    int SurfL[] = {
            6, 1, 1, 0, 0,
            2, 1, 1, 0, 0
    };
    int SurfI[] = {
            4, 1,  3, 1,  1, 0,  2, 0,  5, 0, 1, 1,
            3, 0,  4, 0
    };
    double SurfT[4*2] = {-1.1, 1.1, -1.1, 1.1,  -1.1, 1.1, -1.1, 1.1}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [rt = theta](int i, double t, double *pu, double *pv){
        switch(i){
            case 1: { *pu = t   , *pv = 0.0      ; break; }
            case 2: { *pu = cos(t)   , *pv = sin(t)   ; break; }
            case 3: { *pu = rt*cos(t), *pv = rt*sin(t); break; }
        }
        return 1;
    };
    asd.bounsurf = [R, K](int i, double u, double v, double *px, double *py, double *pz) {
        double w0 = max(1 - K*K*(u*u + v*v), 0.0), wc = max(1 - K*K, 0.0);
        double w = max(sqrt(w0) - sqrt(wc), 0.0);
        double z = (w < 100*DBL_EPSILON) ? 0.0 : ((R * w) / K);
        if (fabs(u) < DBL_EPSILON) u = 0;
        if (fabs(v) < DBL_EPSILON) v = 0;

        *px = R * u,  *py = R * v, *pz = ((i % 2) ? 1 : -1) * z;
        return 1;
    };
    makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
    am.face_label.resize(am.faces.size() / 3, 0);
    double max_r = R * theta;
    for (int i = 0; i < am.faces.size() / 3; ++i){
        double cx = 0, cy = 0;
        for (int n = 0; n < 3; ++n) 
            cx += am.vertices[3*(am.faces[3*i+n]-1) + 0], cy += am.vertices[3*(am.faces[3*i+n]-1) + 1];
        cx /= 3, cy /= 3;
        double cr = sqrt(cx*cx + cy*cy);
        if (cr < max_r) am.face_label[i] = 4;
        else am.face_label[i] = 2;
    }

    mark_vertices(am);

    return am;
}
