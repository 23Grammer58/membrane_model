//
// Created by alex on 07.05.2021.
//

#ifndef AORTIC_VALVE_FORANI3DFRTPRM_H
#define AORTIC_VALVE_FORANI3DFRTPRM_H

#include <functional>
#include <vector>
#include "../TriangularMeshHelpers.h"

//WARNING: vectors edgeout and edgecolor used or not used at the same time
// and if used then expectNE means expected number of edges and
// edgeout->resize(2*expectNE), edgecolor->resize(expectNE) will be called
// if expectNE will be too small then mesh generating will be done two times
// optional input: exportCurves, vu_map, facecolor, edgeout, edgecolor
int makeAft3dSurface (
        int nVVert, const double *VVertxyz,
        int nLine, const int *LineD, const int *LineP, const double *LineT, const int *exportCurves,
        int nSurface, const int *SurfL, const int *SurfI, const double *SurfT,
        std::function<void(int lbl, double u, double* pv)> v_u,
        std::function<int(int lbl, double u, double v, double* px, double* py, double* pz)> bounsurf,
        std::function<int(int lbl, double t, double* pu, double* pv)> bounline,
        //dir = 0 for u-direction and dir = 1 for v-direction, return period or 0.0 if direction unperiodical
        std::function<double(int iSurf, int dir)> periodicfunction,
        std::function<double(double x, double y, double z)> fsize,
        std::vector<double>& vertexout, std::vector<double>* vu_map,
        std::vector<int>& faceout, std::vector<int>* facecolor,
        int expectNE = 0, std::vector<int>* edgeout = nullptr, std::vector<int>* edgecolor = nullptr,
        int indexshift = 0
);

struct Ani3dSurfDiscr{
    int nVVert = 0;   const double *VVertxyz = nullptr;
    int nLine = 0;    const int *LineD = nullptr, *LineP = nullptr; const double *LineT = nullptr; const int *LineLbl = nullptr;
    int nSurface = 0; const int *SurfL = nullptr, *SurfI = nullptr; const double *SurfT = nullptr;
    std::function<void(int, double, double*)> v_u = nullptr;
    std::function<int(int, double, double, double*, double*, double*)> bounsurf = nullptr;
    std::function<int(int, double, double*, double*)> bounline = nullptr;
    std::function<double(int, int)> periodicfunction = nullptr;
};

struct Ani3dMeshOut{
    std::vector<double>& vertex;
    std::vector<int>& face;
    std::vector<int>* facecolor = nullptr;
    int expectNE = 600000;
    std::vector<int>* edge = nullptr;
    std::vector<int>* edgecolor = nullptr;
    std::vector<double>* vu_map = nullptr;
    Ani3dMeshOut(std::vector<double>& vertex, std::vector<int>& face, std::vector<int>* facecolor = nullptr,
                 int expectNE = 600000, std::vector<int>* edge = nullptr, std::vector<int>* edgecolor = nullptr,
                 std::vector<double>* vu_map = nullptr):
                 vertex(vertex), face(face), facecolor(facecolor),
                 expectNE(expectNE), edge(edge), edgecolor(edgecolor), vu_map(vu_map) {};
    Ani3dMeshOut(AniMesh& am): vertex{am.vertices}, face{am.faces}, facecolor{&am.face_label}, edge(&am.edges), edgecolor(&am.edge_label) {}
};

static int makeAft3dSurface (const Ani3dSurfDiscr& d, std::function<double(double x, double y, double z)> fsize, Ani3dMeshOut out, int indexshift = 1){
    return makeAft3dSurface(d.nVVert, d.VVertxyz,
                     d.nLine, d.LineD, d.LineP, d.LineT, d.LineLbl,
                     d.nSurface, d.SurfL, d.SurfI, d.SurfT,
                     d.v_u, d.bounsurf, d.bounline, d.periodicfunction, std::move(fsize),
                     out.vertex, out.vu_map, out.face, out.facecolor, out.expectNE, out.edge, out.edgecolor, indexshift);
}

typedef int PointID;
typedef int LineID;
typedef int SurfID;
using AniSurfPoint = Point;
struct AniSurfLine{
    PointID st = PointID(), end = PointID();
    bool is_reverse = false;
    bool use_func = false;
    int func_param = 0;
    double st_t = 0, end_t = 0; //parametrization borders
    int lbl = 0;
};
struct AniSurfSurf{
    std::vector<LineID> lines;
    bool use_func = false;
    int func_param = 0;
    int back_color = 0;
    bool is_positive_orient = true;
    double st_u = 0, end_u = 0;
    double st_v = 0, end_v = 0;
};
struct Ani3dSurfDiscrWrap{
    std::vector<double> VVertxyz;
    std::vector<int> LineD, LineP;
    std::vector<double> LineT;
    std::vector<int> LineLbl;
    std::vector<int> SurfL, SurfI;
    std::vector<double> SurfT;
    std::function<void(int, double, double*)> v_u;
    std::function<int(int, double, double*, double*)> line_param;
    std::function<int(int, double, double, double*, double*, double*)>  surf_param;
    std::function<double(int, int)> periodicfunction;

    Ani3dSurfDiscr getAniSurfDiscr() const {
        return Ani3dSurfDiscr{static_cast<int>(VVertxyz.size()/3), VVertxyz.data(),
                            static_cast<int>(LineD.size()/3), LineD.data(), LineP.data(), LineT.data(), LineLbl.data(),
                            static_cast<int>(SurfL.size()/5), SurfL.data(), SurfI.data(), SurfT.data(),
                            v_u, surf_param, line_param, periodicfunction};
    }
    template<typename POINT>
    PointID addPoint(const POINT& p) {
        VVertxyz.resize(VVertxyz.size()+3);
        for (int i = 0; i < 3; ++i) VVertxyz[VVertxyz.size()-3+i] = p[i];
        return static_cast<int>(VVertxyz.size() / 3);
    }
    PointID addPoint() {
        VVertxyz.resize(VVertxyz.size()+3, 0.0);
        return static_cast<int>(VVertxyz.size() / 3);
    }
    LineID addLine(const AniSurfLine& asl){
        LineD.push_back(asl.st), LineD.push_back(asl.end), LineD.push_back(asl.is_reverse ? 0 : 1);
        if (asl.use_func){
            LineP.push_back(1), LineP.push_back(asl.func_param);
            LineT.push_back(asl.st_t), LineT.push_back(asl.end_t);
        } else {
            LineP.push_back(0), LineP.push_back(0);
            LineT.push_back(0), LineT.push_back(0);
        }
        LineLbl.push_back(asl.lbl);
        return static_cast<int>(LineD.size()/3);
    }
    SurfID addSurface(const AniSurfSurf& a){
        SurfL.push_back(a.lines.size());
        if (a.use_func) SurfL.push_back(1), SurfL.push_back(a.func_param);
        else SurfL.push_back(0), SurfL.push_back(0);
        SurfL.push_back(a.back_color);
        SurfL.push_back(a.is_positive_orient ? 0 : 1);
        SurfT.push_back(a.st_u), SurfT.push_back(a.end_u);
        SurfT.push_back(a.st_v), SurfT.push_back(a.end_v);
        for (int i = 0; i < a.lines.size(); ++i)
            if (a.lines[i] >= 0) SurfI.push_back(a.lines[i]), SurfI.push_back(0);
            else SurfI.push_back(-a.lines[i]), SurfI.push_back(1);

        return static_cast<int>(SurfL.size()/5);
    }
    Ani3dSurfDiscrWrap& set_LineParamFcn(const std::function<int(int, double, double*, double*)>& lineParam) { return line_param = lineParam, *this; }
    Ani3dSurfDiscrWrap& set_LineVUFnc(const std::function<void(int, double, double*)>& lineParam) { return v_u = lineParam, *this; }
    Ani3dSurfDiscrWrap& set_SurfParamFcn(const std::function<int(int, double, double, double*, double*, double*)>& surfParam) { return surf_param = surfParam, *this; }
    Ani3dSurfDiscrWrap& set_PeriodicFcn(const std::function<double(int, int)> period){ return periodicfunction = period, *this; }
};

struct Ani3dFSize{
    std::function<double(double x, double y, double z)> fsize;
    Ani3dFSize(std::function<double(double x, double y, double z)> fsize): fsize(fsize) {};
    Ani3dFSize(double size): fsize([size](double, double, double){return size;}) {}
};

static int makeAft3dSurface(const Ani3dSurfDiscrWrap& asdw, Ani3dFSize fsize, AniMesh& am){
    return makeAft3dSurface(asdw.getAniSurfDiscr(), fsize.fsize, Ani3dMeshOut(am), 1);
}


#endif //AORTIC_VALVE_FORANI3DFRTPRM_H
