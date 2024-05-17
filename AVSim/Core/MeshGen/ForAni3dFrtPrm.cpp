//
// Created by alex on 07.05.2021.
//

#include "ForAni3dFrtPrm.h"

#ifdef USE_ANI3D
extern "C" {
#include "aniAFT/src/aniFRT/PRM/common.h"
#include "aniAFT/src/aniFRT/PRM/memory3.h"
#include "aniAFT/src/aniFRT/PRM/region3.h"
}
#include <cstring>

static std::function<void(int, double, double*)> m_v_u;
static std::function<int(int, double, double, double*, double*, double*)> m_bounsurf;
static std::function<int(int, double, double*, double*)> m_bounline;
static std::function<double(int, int)> m_periodicfunction;
static std::function<double(double, double, double)> m_fsize;

static void ssv_u(int lbl, double v, double* pu) { return m_v_u(lbl, v, pu); }
static int ssbounsurf(int lbl, double u, double v, double* px, double* py, double* pz) { return m_bounsurf(lbl, u, v, px, py, pz); }
static int ssbounline(int lbl, double t, double* pu, double* pv) { return m_bounline(lbl, t, pu, pv); }
//dir = 0 for u-direction and dir = 1 for v-direction, return period or 0.0 if direction unperiodical
static double ssperiodicfunction(int iSurf, int dir) { return m_periodicfunction(iSurf, dir); }
static double ssfsize(double x, double y, double z) { return m_fsize(x, y, z); }


static int saveit(surface_mesh& m, std::vector<double>& vertexout, std::vector<double>* vu_map,
                  std::vector<int>& faceout, std::vector<int>* facecolor, int is){
    auto nv = vertexout.size(), nf = faceout.size();
    auto nvu = vu_map ? vu_map->size() : 0, nfc = facecolor ? facecolor->size() : 0;
    vertexout.resize(nv + m.nPoint*3);
    faceout.resize(nf + m.nFace*3);
    if (vu_map) vu_map->resize(nvu + m.nPoint*2);
    if (facecolor) facecolor->resize(nfc + m.nFace);
    for (int i = 0; i < m.nPoint; ++i){
        vertexout[nv+3*i+0] = m.vert[i].x;
        vertexout[nv+3*i+1] = m.vert[i].y;
        vertexout[nv+3*i+2] = m.vert[i].z;
        if (vu_map){
            (*vu_map)[nfc+2*i+0] = m.vert[i].u;
            (*vu_map)[nfc+2*i+1] = m.vert[i].v;
        }
    }
    for (int i = 0; i < m.nFace; ++i){
        faceout[nf+3*i+0] = m.face[i].v1 + is;
        faceout[nf+3*i+1] = m.face[i].v2 + is;
        faceout[nf+3*i+2] = m.face[i].v3 + is;
        if (facecolor){
            (*facecolor)[nfc+i] = m.face[i].color;
        }
    }
    if (m.exportCurves && m.nexpEdge){
        for (int i=0; i< m.nexpEdge; i++) {
            m.expEdge[2*i+0] += is;
            m.expEdge[2*i+1] += is;
        }
    }

    return 0;
}
#endif


int makeAft3dSurface (
        int nVVert, const double *VVertxyz,
        int nLine, const int *LineD, const int *LineP, const double *LineT, const int *exportCurves,
        int nSurface, const int *SurfL, const int *SurfI, const double *SurfT,
        std::function<void(int lbl, double u, double* pv)> v_u,
        std::function<int(int lbl, double u, double v, double* px, double* py, double* pz)> bounsurf,
        std::function<int(int lbl, double t, double* pv, double* pu)> bounline,
        //dir = 0 for u-direction and dir = 1 for v-direction, return period or 0.0 if direction unperiodical
        std::function<double(int iSurf, int dir)> periodicfunction,
        std::function<double(double x, double y, double z)> fsize,
        std::vector<double>& vertexout, std::vector<double>* vu_map,
        std::vector<int>& faceout, std::vector<int>* facecolor,
        int enE, std::vector<int>* edgeout, std::vector<int>* edgecolor,
        int indexshift)
{
#ifdef USE_ANI3D
    m_v_u = std::move(v_u);
    m_bounsurf = std::move(bounsurf);
    m_bounline = std::move(bounline);
    m_periodicfunction = std::move(periodicfunction);
    m_fsize = std::move(fsize);

    surface_mesh mesh;
    int r = 0;

    memset(&mesh, 0, sizeof(mesh));
    mesh.v_u      = m_v_u ? ssv_u : nullptr;
    mesh.bounsurf = m_bounsurf ? ssbounsurf : nullptr;
    mesh.bounline = m_bounline ? ssbounline : nullptr;
    mesh.periodicfunction = m_periodicfunction ? ssperiodicfunction : nullptr;
    mesh.fsize    = m_fsize ? ssfsize : nullptr;
    mesh.exportCurves = const_cast<int*>(exportCurves);
    if (edgeout && edgecolor && exportCurves){
        edgeout->resize(enE*2);
        edgecolor->resize(enE);
        mesh.nexpEdge = 0;
        mesh.expEdge = edgeout->data();
        mesh.expEdgeColor = edgecolor->data();
        mesh.nnexpEdge = enE;
    } else {
        mesh.nexpEdge = -1;
        mesh.expEdge = nullptr;
        mesh.expEdgeColor = nullptr;
        mesh.nnexpEdge = 0;
    }
    initMemory(&mesh);

    auto cVVertxyz = const_cast<double*>(VVertxyz);
    auto cLineD = const_cast<int*>(LineD), cLineP = const_cast<int*>(LineP); auto cLineT = const_cast<double*>(LineT);
    auto cSurfL = const_cast<int*>(SurfL), cSurfI = const_cast<int*>(SurfI); auto cSurfT = const_cast<double*>(SurfT);
    initAFS_(&mesh, &nVVert, cVVertxyz, &nLine, cLineD, cLineP, cLineT, &nSurface, cSurfL, cSurfI, cSurfT);
    if (edgeout && edgecolor && exportCurves && enE < mesh.nexpEdge){
        enE = mesh.nexpEdge;
        freeMemory(&mesh);
        memset(&mesh, 0, sizeof(mesh));
        mesh.v_u      = m_v_u ? ssv_u : nullptr;
        mesh.bounsurf = m_bounsurf ? ssbounsurf : nullptr;
        mesh.bounline = m_bounline ? ssbounline : nullptr;
        mesh.periodicfunction = m_periodicfunction ? ssperiodicfunction : nullptr;
        mesh.fsize    = m_fsize ? ssfsize : nullptr;
        mesh.exportCurves = const_cast<int*>(exportCurves);
        edgeout->resize(enE*2);
        edgecolor->resize(enE);
        mesh.nexpEdge = 0;
        mesh.expEdge = edgeout->data();
        mesh.expEdgeColor = edgecolor->data();
        mesh.nnexpEdge = enE;
        initMemory(&mesh);
        initAFS_(&mesh, &nVVert, cVVertxyz, &nLine, cLineD, cLineP, cLineT, &nSurface, cSurfL, cSurfI, cSurfT);
    }
    r = saveit(mesh, vertexout, vu_map, faceout, facecolor, indexshift);

    freeMemory(&mesh);
    return r;
#else
    throw std::runtime_error("Ani3D library is not included to project");
    return -1;
#endif
}