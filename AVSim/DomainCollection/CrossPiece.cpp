//
// Created by alex on 01.12.2023.
//
#include "AVSim/Core/MeshGen/Helper.h"

static int create_cross_piece(AniMesh& am, std::array<double, 2> a, std::array<double, 2> d, std::array<double, 4> r, double size){
    int nVert = 8+5+8;
    int nLine = 12+4*4;
    int nSurface = 8;
    std::array<double, 8> c = {
      a[1]/2 - r[0], a[0]/2 - r[0], a[0]/2 - r[1], a[1]/2 - r[1],
      a[1]/2 - r[2], a[0]/2 - r[2], a[0]/2 - r[3], a[1]/2 - r[3],
    };
    assert(a[0] > d[0] + size && a[1] > d[1] + size && "Too big central region");
    assert(((a[0] - d[0])/2)*((a[0] - d[0])/2) + ((a[1] - d[1])/2)*((a[1] - d[1])/2) > max({r[0], r[1], r[2], r[3]}) && "Too big radius r");
    double VVert[] = {
         a[0]/2,    c[0], 0,
           c[1],  a[1]/2, 0,
          -c[2],  a[1]/2, 0,
        -a[0]/2,    c[3], 0,
        -a[0]/2,   -c[4], 0,
          -c[5], -a[1]/2, 0,
           c[6], -a[1]/2, 0,
         a[0]/2,   -c[7], 0,

         a[0]/2,    0, 0,
           0,  a[1]/2, 0,
        -a[0]/2,    0, 0,
           0, -a[1]/2, 0,

           0,    0, 0,

         d[0]/2,    0, 0,
         d[0]/2,  d[1]/2, 0,
           0,  d[1]/2, 0,
        -d[0]/2,  d[1]/2, 0,
        -d[0]/2,    0, 0,
        -d[0]/2, -d[1]/2, 0, 
           0, -d[1]/2, 0, 
         d[0]/2, -d[1]/2, 0,

    };
    int LineD[] = {
        1, 2, 1,  3, 4, 1,  5, 6, 1,  7, 8, 1,
        9, 1, 1,  2, 10, 1,  10, 3, 1,  4, 11, 1,  11, 5, 1,  6, 12, 1,  12, 7, 1,  8, 9, 1,
        9, 14, 1,  14, 13, 1,   10, 16, 1,  16, 13, 1,   11, 18, 1,  18, 13, 1,   12, 20, 1,  20, 13, 1,
        14, 15, 1,   15, 16, 1,  16, 17, 1,  17, 18, 1,  18, 19, 1,  19, 20, 1,  20, 21, 1,  21, 14, 1
    };
    int exportCurves[] = {16, 16, 16, 16, 1, 2, 2, 4, 4, 8, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 32, 32, 32, 32, 32, 32, 32, 32};
    int LineP[] = { 1, 1,  1, 2,  1, 3,  1, 4, 
                    0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  
                    0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,
                    0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0, };
    double LineT[] = {3*M_PI_2, M_PI,  2*M_PI, 3*M_PI_2,  M_PI_2, 0,  M_PI, M_PI_2,
                      0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
                      0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
                      0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,}; 
    int SurfL[] = { 7, 0, 1, 0, 0,  7, 0, 1, 0, 0,  7, 0, 1, 0, 0,  7, 0, 1, 0, 0,
                    4, 0, 1, 0, 0,  4, 0, 1, 0, 0,  4, 0, 1, 0, 0,  4, 0, 1, 0, 0};
    int SurfI[] = {
         5,0, 1,0,  6,0, 15,0, 22,1, 21,1, 13,1,
         7,0, 2,0,  8,0, 17,0, 24,1, 23,1, 15,1,
         9,0, 3,0, 10,0, 19,0, 26,1, 25,1, 17,1,
        11,0, 4,0, 12,0, 13,0, 28,1, 27,1, 19,1,
        21,0, 22,0, 16,0, 14,1,
        23,0, 24,0, 18,0, 16,1,
        25,0, 26,0, 20,0, 18,1,
        27,0, 28,0, 14,0, 20,1,
    };
    double A = max(a[0], a[1]);
    double SurfT[] = {  -A, -A, A, A,  -A, -A, A, A,  -A, -A, A, A,  -A, -A, A, A,
                        -A, -A, A, A,  -A, -A, A, A,  -A, -A, A, A,  -A, -A, A, A};
    Ani3dSurfDiscr asd{nVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    asd.bounline = [r, a](int i, double t, double *pu, double *pv){ 
        *pu = ((i == 2 || i == 3) ? -1 : 1) *a[0]/2 + r[i-1] * cos(t);
        *pv = (i > 2 ? -1 : 1)*a[1]/2 + r[i-1] * sin(t);
        return 1;
    };
    asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1;};

    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_cross_piece(std::array<double, 2> a, std::array<double, 2> ac, std::array<double, 4> r, double size){
    AniMesh am;
    create_cross_piece(am, a, ac, r, size);
    mark_vertices(am);

    am.face_label.resize(am.faces.size() / 3, 0);
    for (int i = 0; i < am.faces.size() / 3; ++i){
        double cx = 0, cy = 0;
        for (int n = 0; n < 3; ++n) 
            cx += am.vertices[3*(am.faces[3*i+n]-1) + 0], cy += am.vertices[3*(am.faces[3*i+n]-1) + 1];
        cx /= 3, cy /= 3;
        if (abs(cx) <= ac[0]/2 && abs(cy) <= ac[1]/2) am.face_label[i] = 4;
        else am.face_label[i] = 2;
    }

    return am;
}

AniMesh generate_cross_piece(double a, double ac, double r, double size){ return generate_cross_piece({a, a}, {ac, ac}, {r, r, r, r}, size); }