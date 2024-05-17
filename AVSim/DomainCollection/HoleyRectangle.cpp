//
// Created by alex on 01.04.2021.
//
#include "AVSim/Core/MeshGen/Helper.h"

using Holes = vector<pair<array<double, 2>, double>>;
static tuple<bool, Holes, Holes> separate_holes(double a, double b, const Holes& holes){
    auto dif = [](array<double, 2> p1, array<double, 2> p2){
        return array<double, 2>{p1[0] - p2[0], p1[1] - p2[1]};
    };
    auto len = [](array<double, 2> p){
        return sqrt(p[0] * p[0] + p[1] * p[1]);
    };
    for(unsigned i = 0; i < holes.size(); ++i)
        for (unsigned j = i+1; j < holes.size(); ++j){
            const auto& h1 = holes[i], h2 = holes[j];
            const auto& p1 = h1.first, p2 = h2.first;
            if (len(dif(p1, p2)) < h1.second + h2.second)
                return tuple<bool, Holes, Holes>{false, Holes(), Holes()};
        }
    auto check_in = [&a, &b](const pair<array<double, 2>, double>& h){
        const auto& p = h.first;
        const double r = h.second;
        const double eps = min (1e-8, min(a,b) * 1e-5);
        return (p[0]-r > eps && p[1]-r > eps && a - (p[0] + r) > eps && b - (p[1] + r) > eps);
    };
    auto check_intersect = [&a, &b, &len, &dif](const pair<array<double, 2>, double>& h){
        const auto& p = h.first;
        const double r = h.second;
        const double eps = min (1e-8, min(a,b) * 1e-5);
        array<double, 2> v0{0.0, 0.0}, v1{a, b};
        if (p[0] > v0[0] && p[1] > v0[1] && p[0] < v1[0] && p[1] < v1[1])
            return true;
        if (p[0] < v0[0] && p[1] < v0[1]) return (len(dif(p, v0)) < r - eps);
        if (p[0] > v1[0] && p[1] > v1[1]) return (len(dif(p, v1)) < r - eps);
        if (p[0] < v0[0] && p[1] > v1[1]) return (len(dif(p, {v0[0], v1[1]})) < r - eps);
        if (p[0] > v1[0] && p[1] < v0[1]) return (len(dif(p, {v1[0], v0[1]})) < r - eps);
        if (p[0] <= v0[0] + DBL_EPSILON) return fabs(p[0] - v0[0]) < r - eps;
        if (p[0] >= v1[0] - DBL_EPSILON) return fabs(p[0] - v1[0]) < r - eps;
        if (p[1] <= v0[1] + DBL_EPSILON) return fabs(p[1] - v0[1]) < r - eps;
        if (p[1] >= v1[1] - DBL_EPSILON) return fabs(p[1] - v1[1]) < r - eps;
        cout << "check_inserter: something wrong\n";
        return false;
    };
    Holes in, brd;
    for (const auto& h: holes) {
        if (check_in(h)) in.push_back(h);
        else if (check_intersect(h)) brd.push_back(h);
    }
    return tuple<bool, Holes, Holes>{true, in, brd};
}

static int create_rectangle_with_holes(AniMesh& am, double a, double b, tuple<bool, Holes, Holes>& hl, double size)
{
    Holes& Hin = get<1>(hl), Hbrd = get<2>(hl);
    if (Hbrd.size() > 0) throw runtime_error("create_rectangle_with_holes: Boundary intersectors doesn't implemented yet\n");
    int nVVert = 4 + 2 * Hin.size() + 2 * Hbrd.size();
    int nLine = 4 + 2 * Hin.size() + 2 * Hbrd.size();
    int nSurface = 1;
    vector<double> VVert(3*nVVert), LineT(2*nLine);
    vector<int> LineD(3*nLine), LineP(2*nLine), exportCurves(nLine);
    auto to_VVert = [&VVert](double x, double y){
        static int n = 0;
        VVert[n + 0] = x, VVert[n + 1] = y, VVert[n + 2] = 0;
        n+=3;
    };
    auto to_LineD= [&LineD](int v1, int v2){
        static int n = 0;
        LineD[n + 0] = v1, LineD[n + 1] = v2, LineD[n + 2] = 1;
        n+=3;
    };
    auto to_LineP= [&LineP](int func_num, int index){
        static int n = 0;
        LineP[n + 0] = func_num, LineP[n + 1] = index;
        n+=2;
    };
    to_VVert(0, 0); to_LineD(1, 2); to_LineP(0, 0); exportCurves[0] = 2;
    to_VVert(0,b); to_LineD(2, 3); to_LineP(0, 0); exportCurves[1] = 1;
    to_VVert(a, b); to_LineD(3, 4); to_LineP(0, 0); exportCurves[2] = 16;
    to_VVert(a, 0); to_LineD(4, 1); to_LineP(0, 0); exportCurves[3] = 8;
    int off = 5;
    for (const auto& i: Hin){
        auto& p = i.first;
        auto& r = i.second;
        to_VVert(p[0]+r, p[1]); to_LineD(off, off+1); to_LineP(1, (off-5)+1); exportCurves[off-1] = 4;
        to_VVert(p[0]-r, p[1]); to_LineD(off+1, off); to_LineP(1, (off-5)+2); exportCurves[off] = 4;
        off += 2;
    }
    for (const auto& i: Hbrd){
        auto& p = i.first;
        auto& r = i.second;
        to_VVert(p[0]+r, p[1]); to_LineD(off, off+1); to_LineP(1, (off-5)+1); exportCurves[off-1] = 4;
        to_VVert(p[0]-r, p[1]); to_LineD(off+1, off); to_LineP(1, (off-5)+2); exportCurves[off] = 4;
        off += 2;
    }
    for (int i = 0; i < nLine; ++i){
        LineT[2*i] = 0; LineT[2*i + 1] = M_PI;
    }


    int SurfL[] = { nLine /*edges*/, 0 /*param*/, 1, 0 /*backcolor*/, 0 /*direction*/}; /*main surface*/
    vector<int> SurfI(2*nLine);
    for (int i = 1; i <= 4; ++i)
        SurfI[2*(i-1)] = i, SurfI[2*(i-1)+1] = 0;
    for (int i = 5; i <= nLine; ++i)
        SurfI[2*(i-1)] = i, SurfI[2*(i-1)+1] = 0;
    double SurfT[] = {-a, 2*a, -b, 2*b}; /*parametrization bbox*/
    Ani3dSurfDiscr asd{nVVert, VVert.data(), nLine, LineD.data(), LineP.data(), LineT.data(), exportCurves.data(), nSurface, SurfL, SurfI.data(), SurfT};
    asd.bounline = [&g_hl = hl](int i, double t, double *pu, double *pv){
        int id = (i-1)/2, j = (i-1)%2;
        auto& hl = (id > static_cast<int>(get<1>(g_hl).size())) ? get<2>(g_hl)[id - get<1>(g_hl).size()] : get<1>(g_hl)[id];
        auto& p = hl.first;
        double r = hl.second;
        int sign = ((j==0) ? 1 : -1);
        *pu = p[0] + sign*r*cos(t);
        *pv = p[1] + sign*r*sin(t);

        return 1;
    };
    asd.bounsurf = [](int i, double u, double v, double *px, double *py, double *pz){ return *px = u, *py = v, *pz = 0, 1;};

    return makeAft3dSurface(asd, Ani3dFSize(size).fsize, Ani3dMeshOut(am));
}

AniMesh generate_holey_rectangle(double a, double b, double size, const vector<pair<array<double, 2>, double>>& holes){
    tuple<bool, Holes, Holes> hl = separate_holes(a, b, holes);
    if (!get<0>(hl)) throw runtime_error("Mesh shouldn't have intersecting holes\n");
    if (a < 0.5*size || b < 0.5*size) {
        throw runtime_error("Impossible to generate rectangle with such params a = " + to_string(a) + ", b = " + to_string(b) + ", h = " + to_string(size) + "\n");
    }
    AniMesh am;
    create_rectangle_with_holes(am, a, b, hl, size);
    mark_vertices(am);

    return am;
}
