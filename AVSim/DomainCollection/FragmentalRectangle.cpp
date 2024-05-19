//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"

////mesh generator functions//
//static double global_radius = 1.0;
//static double global_size = 0.1;
//static double g_ssa = 10;
//static double g_ssb = 1;
//static double global_len = 1.0;
//static double global_phi = M_PI_2;
//
//static double special_size(double x, double y, double z) {
//    double a = 2.995732, w = 0.05;
//    std::array<double, 3> xx{x, y, z};
//    auto h = [a, w](std::array<double, 3> p0, std::array<double, 3> p1){
//        double l = (p0[0] - p1[0])*(p0[0] - p1[0]) + (p0[1] - p1[1])*(p0[1] - p1[1]) + (p0[2] - p1[2])*(p0[2] - p1[2]);
//        return 9 * exp(-a*l/(w*w));
//    };
//    return global_size / (1 + h(xx, {-g_ssa/2, -g_ssb/2, 0}) + h(xx, {-g_ssa/2, g_ssb/2, 0}));
//}
//
//static double special_size1(double x, double y, double z) {
//    double a = 2.995732, w = 0.05, pw = 0, pz = 5;
//    std::array<double, 3> xx{x, y, z};
//    auto h = [a, w, pw](std::array<double, 3> p0, std::array<double, 3> p1){
//        double l = (p0[0] - p1[0])*(p0[0] - p1[0]) + (p0[1] - p1[1])*(p0[1] - p1[1]) + (p0[2] - p1[2])*(p0[2] - p1[2]);
//        return pw * exp(-a*l/(w*w));
//    };
//    auto hz = [a, w = 1, pz](double p0, double p1){
//        return pz * exp(-a*(p1 - p0) * (p1 - p0)/(w*w));
//    };
//    return global_size / (1 + hz(xx[0], -5) + h(xx, {-5, -0.5, 0}) + h(xx, {-5, 0.5, 0}));
//}

static int create_rectangle1(AniMesh& am, double a, double b, double size){
    int nVVert = 4;
    int nLine = 4;
    int nSurface = 1;
    double VVert[] = {-a/2,-b/2, 0,   a/2, -b/2, 0,   a/2, b/2, 0,   -a/2, b/2, 0};
    int LineD[] = {1, 2, 1,  2, 3, 1,  3, 4, 1,  4, 1, 1}; /*edges*/
    int LineP[] = {0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0};
    double LineT[] = {0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0,  0, 0};
    int exportCurves[] = {1, 2, 1, 4};
    int SurfL[] = { 6 /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    int SurfI[] = {1, 0 /*direction*/,  2, 0 /*direction*/,  3, 0 /*direction*/,  4, 0,  5,0,  6,0};
    double SurfT[] = {0, 0, 0, 0,
                      0, 0, 0, 0,}; /*parametrization bbox*/

    Ani3dSurfDiscr asd{nVVert, VVert, nLine, LineD, LineP, LineT, exportCurves, nSurface, SurfL, SurfI, SurfT};
    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));

}

static int create_rectangle(AniMesh& am, double a, double b, double size, int fragments = 1){
    std::vector<std::array<double, 3>> VVert;
    double dx = a / fragments;
    for (int i = 0; i <= fragments; ++i){
        VVert.push_back({-a/2 + dx*i, b/2, 0});
        VVert.push_back({-a/2 + dx*i,-b/2, 0});
    }
    int N = VVert.size();
    VVert.push_back({-a/2, 0, 0});
    VVert.push_back({a/2, 0, 0});
    std::vector<std::array<int, 3>> LineD;
    std::vector<int> exportCurves;
    LineD.push_back({1, N+1, 1}); exportCurves.push_back(4);
    LineD.push_back({N+1, 2, 1}); exportCurves.push_back(4);
    for (int i = 1; i < fragments; ++i){
        LineD.push_back({2*i +1, 2*i + 2, 1}); exportCurves.push_back(0);
    }
    LineD.push_back({N-1, N+2, 1}); exportCurves.push_back(2);
    LineD.push_back({N+2, N, 1}); exportCurves.push_back(2);
    int K = LineD.size();
    for (int i = 0; i < fragments; ++i){
        LineD.push_back({2*i+1, 2*(i+1)+1, 1}); exportCurves.push_back(1);
        LineD.push_back({2*(i+1)+2, 2*(i)+2, 1}); exportCurves.push_back(1);
    }
    std::vector<int> LineP(LineD.size()*2, 0);
    std::vector<double> LineT(LineD.size()*2, 0);
    std::vector<std::array<int, 5>> SurfL;
    std::vector<int> SurfI;
    auto back_insert = [](auto& SurfI, auto ins){ std::for_each(ins.begin(), ins.end(), [&SurfI](const auto& i) mutable { SurfI.push_back(i); }); };
    if (fragments == 1){
        SurfL.push_back({6, 0, 1, 0, 0});
        back_insert(SurfI, std::array<int, 12>{2, 1, 1, 1, K+1, 0, 3, 0, 4, 0, K+2, 0 });
    } else {
        SurfL.push_back({5, 0, 1, 0, 0});
        back_insert(SurfI, std::array<int, 10>{2, 1, 1, 1, K+1, 0, 3, 0, K+2, 0 });
        for (int i = 1; i < fragments-1; ++i){
            SurfL.push_back({4, 0, 1, 0, 0});
            back_insert(SurfI, std::array<int, 8>{i+2, 1, K+2*i+1, 0, i+3, 0, K+2*i+2, 0});
        }
        SurfL.push_back({5, 0, 1, 0, 0});
        int KK = LineD.size();
        back_insert(SurfI, std::array<int, 10>{K-2, 1, KK-1, 0, K-1, 0, K, 0, KK, 0});
    }
    std::vector<double> SurfT(4*SurfL.size(), 0.0);
    int nVVert = VVert.size();
    int nLine = LineD.size();
    int nSurface = SurfL.size();

    Ani3dSurfDiscr asd{nVVert, reinterpret_cast<double*>(VVert.data()), nLine, reinterpret_cast<int*>(LineD.data()),
                     LineP.data(), LineT.data(), exportCurves.data(), nSurface, reinterpret_cast<int*>(SurfL.data()), SurfI.data(), SurfT.data()};
    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_rectangle(double a, double b, double size, int fragments){
    if (a < 0.7*size || b < 0.7*size) {
        throw runtime_error("Impossible to generate rectangle with such params a = " + to_string(a) + ", b = " + to_string(b) + ", h = " + to_string(size) + "\n");
    }
    AniMesh am;
    create_rectangle(am, a, b, size, fragments);
    mark_vertices(am);

    set<std::array<int, 2>> bnd;
    for (int i = 0; i < am.edges.size()/2; ++i)
        bnd.insert({am.edges[2*i], am.edges[2*i+1]}), bnd.insert({am.edges[2*i+1], am.edges[2*i]});

    std::array<std::array<int, 2>, 4> wfi = {0};
    std::array<int, 4> wfid = {0};
    for (int f = 0, wff = 0; f < am.faces.size()/3; ++f){
        int ebnd = 0;
        int eind[3] = {0};
        for (int k = 0; k < 3; ++k)
            if (bnd.count({am.faces[3*f+k%3], am.faces[3*f+(k + 1)%3]})) {
                eind[ebnd++] = k;
            }
        if (ebnd >= 2) {
            wfi[wff][0] = f;
            if (eind[0] == 0 && eind[1] == 1) wfid[wff] = 1;
            else if (eind[0] == 1 && eind[1] == 2) wfid[wff] = 2;
            else if (eind[0] == 0 && eind[1] == 2) wfid[wff] = 0;
            wff++;
        }
    }

    for (int i = 0; i < am.faces.size()/3; ++i){
        std::set<int> lf{am.faces[3*i], am.faces[3*i+1], am.faces[3*i+2]};
        for (int j = 0; j < 4; ++j){
            if (    lf.count(am.faces[wfi[j][0]*3 + (wfid[j] + 1)%3])
                    &&  lf.count(am.faces[wfi[j][0]*3 + (wfid[j] + 2)%3])
                    && !lf.count(am.faces[wfi[j][0]*3 + (wfid[j] + 0)%3])){
                wfi[j][1] = i;
            }
        }
    }

    std::set<std::array<int, 4>> tmp;     
    for (int l = 0; l < 4; ++l){
        int lind[4] = {0};
        for (int k = 0; k < 3; ++k){
            lind[k] = am.faces[wfi[l][0]*3 + (wfid[l] + k)%3];
        }
        for (int k = 0; k < 3; ++k) {
            if (am.faces[wfi[l][1]*3 + k] != lind[1] && am.faces[wfi[l][1]*3 + k] != lind[2])
                lind[3] = am.faces[3 * wfi[l][1] + k];
        }
        {   //TODO: work around a strange bug appear for ./model stretch_rectangle 100 100 3 4 0 0 0 1.57079633
            std::array<int, 4> alind{lind[0], lind[1], lind[2], lind[3]};
            auto r0 = tmp.size();
            std::sort(alind.begin(), alind.end());
            tmp.insert(alind);
            if (r0 == tmp.size()) continue;
        }
        am.faces[3 * wfi[l][0] + 0] = lind[0], am.faces[3 * wfi[l][0] + 1] = lind[1], am.faces[3 * wfi[l][0] + 2] = lind[3];
        am.faces[3 * wfi[l][1] + 0] = lind[0], am.faces[3 * wfi[l][1] + 1] = lind[3], am.faces[3 * wfi[l][1] + 2] = lind[2];
    }

    return am;
}