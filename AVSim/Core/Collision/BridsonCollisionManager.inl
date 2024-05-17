#ifndef AVSIM_BRIDSONCOLLISIONMANAGERINL_H
#define AVSIM_BRIDSONCOLLISIONMANAGERINL_H

using namespace World3d;

#define TESTVFDIST \
                    auto np1 = np;\
                    bool is_removed = !np1.InitVF(fd, vv, k, p);\
                    if (!is_removed){\
                        int nTri = k, nPonTri = p;\
                        auto ft = fd[(nTri+1)%2]->id;\
                        auto vp = vv[nTri][nPonTri];\
                        auto& vt = vv[(nTri+1)%2];\
                        auto& t_x = fd[(nTri+1)%2]->backptr->m_x, &p_x = fd[nTri]->backptr->m_x;\
                        auto& t_dx = fd[(nTri+1)%2]->backptr->m_dx, &p_dx = fd[nTri]->backptr->m_dx;\
                        std::cout << "T: \n" << std::hexfloat <<\
                                    "\tv0: " << t_x[vt[0]] << " -> " << t_dx[vt[0]] << "\n"\
                                    "\tv1: " << t_x[vt[1]] << " -> " << t_dx[vt[1]] << "\n"\
                                    "\tv2: " << t_x[vt[2]] << " -> " << t_dx[vt[2]] << "\n"\
                                    "P: " << p_x[vp] << " -> " << p_dx[vp] << "\n";\
                        DReal ctm = Opt::GetColTime(np1);\
                        Point   ac = t_x[vt[0]] + ctm * t_dx[vt[0]], \
                                bc = t_x[vt[1]] + ctm * t_dx[vt[1]],\
                                cc = t_x[vt[2]] + ctm * t_dx[vt[2]],\
                                pc = p_x[ vp  ] + ctm * p_dx[ vp  ];\
                        auto ww = GeomProjs::BaryCoord(ac, bc, cc, pc);\
                        auto pc_proj = GeomProjs::BaryEval(ac, bc, cc, ww);\
                        auto sqd = (pc_proj - pc).squared_length();\
                        std::cout << "t_eval = " << std::defaultfloat << std::setprecision(16) << Opt::GetColTime(np) << " changed to " << ctm << " where sqd = " << sqd << "\n"\
                                  << "T: \n"\
                                  << "v0: " << ac << "\n"\
                                  << "v1: " << bc << "\n"\
                                  << "v2: " << cc << "\n"\
                                  << "P: " << pc << "\n"\
                                  << "P projection: " << pc_proj << "\n";\
                        std::cout << "In plane test: " << CGAL::cross_product(bc - ac, cc - ac) * (pc - ac) << "\n";\
                        std::cout << "In triangle test: " \
                            << CGAL::cross_product(ac - pc, bc - pc) * CGAL::cross_product(ac - cc, bc - cc) << " " \
                            << CGAL::cross_product(bc - pc, cc - pc) * CGAL::cross_product(bc - ac, cc - ac) << " " \
                            << CGAL::cross_product(cc - pc, ac - pc) * CGAL::cross_product(cc - bc, ac - bc) << "\n"; \
                    }

template<typename NarrowPhase, typename Handler, typename Opt>
int BridsonCollisionManager::_performBroadPhaseSelf(VolumeTreeBase& a){
    NarrowPhase np;
    Handler handler;
    return _performBroadPhaseSelfExt<NarrowPhase, Handler, Opt>(a, np, handler);
}
template<typename NarrowPhase, typename Handler, typename Opt>
int BridsonCollisionManager::_performBroadPhaseSelfExt(VolumeTreeBase& a, NarrowPhase& np, Handler& handler){
    CollisionStatistics ncols;
    VolumeIntersect(a, a, [this, &ncols, &np, &handler](const VolumeTreeBase::LeafIndex l0, void* data0, const VolumeTreeBase::LeafIndex l1, void* data1)->void{
        auto *fd0 = static_cast<DynamicObject3DCI::FaceData*>(data0),
                *fd1 = static_cast<DynamicObject3DCI::FaceData*>(data1);
        if (fd0 == fd1) return;
        DynamicObject3DCI::FaceData* fd[2] = {fd0, fd1};
        std::array<V_ind, 3> vv[2] = {vert_around(fd[0]->backptr->m_obj->m_mesh, fd[0]->id), vert_around(fd[1]->backptr->m_obj->m_mesh, fd[1]->id)};
        std::array<char, 3> link[2];
        int nocommon = 3;
        {
            for (char i = 0; i < 3; ++i)
                link[0][i] = link[1][i] = -1;
            for (char i = 0; i < 3; ++i) {
                for (char j = 0; j < 3; ++j)
                    if (vv[0][i] == vv[1][j]) {
                        link[0][i] = j;
                        link[1][j] = i;
                        --nocommon;
                        break;
                    }
            }
        }
        // double thickness = fd[0]->backptr->m_thickness;
        // double sqr_thickness = thickness * thickness;

        for (int k = 0; k < 2; ++k){
            for (int p = 0; p < 3; ++p)
                if (link[k][p] < 0){
                    ++ncols.nVFBroadCases;
                    bool is_close = np.InitVF(fd, vv, k, p);
                    if (!is_close) continue;
                    handler.InitVF(np, this);
                    if (handler.applyVF())
                        ++ncols.nVFNarrowCases;
                    CHECKCOLIFBEGIN()
                        //TESTVFDIST
                        //assert(is_removed && "VF collision is not removed");
                    CHECKCOLIFEND();
                }
        }
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j){
                if ((nocommon != 3) && (link[0][i] == j || link[0][i] == (j + 1)%3 || link[0][(i+1)%3] == j || link[0][(i+1)%3] == (j + 1)%3)) continue;
                ++ncols.nEEBroadCases;
                bool is_close = np.InitEE(fd, vv, i, j);
                if (!is_close) continue;
                handler.InitEE(np, this);
                if (handler.applyEE())
                    ++ncols.nEENarrowCases;
                // CHECKCOLIF(assert(!np.InitEE(fd, vv, i, j) && "EE collision is not removed"));    
            }
    });
    return ncols.nVFNarrowCases + ncols.nEENarrowCases;
}
template<typename NarrowPhase, typename Handler, typename Opt>
int BridsonCollisionManager::_performBroadPhaseDiff(VolumeTreeBase& a, VolumeTreeBase& b){
    NarrowPhase np;
    Handler handler;
    return _performBroadPhaseDiffExt<NarrowPhase, Handler, Opt>(a, b, np, handler);
}
template<typename NarrowPhase, typename Handler, typename Opt>
int BridsonCollisionManager::_performBroadPhaseDiffExt(VolumeTreeBase& a, VolumeTreeBase& b, NarrowPhase& np, Handler& handler){
    CollisionStatistics ncols;
    VolumeIntersect(a, b, [this, &ncols, &np, &handler](const VolumeTreeBase::LeafIndex l0, void* data0, const VolumeTreeBase::LeafIndex l1, void* data1)->void{
        auto *fd0 = static_cast<DynamicObject3DCI::FaceData*>(data0),
                *fd1 = static_cast<DynamicObject3DCI::FaceData*>(data1);
        if (fd0 == fd1) return;
        DynamicObject3DCI::FaceData* fd[2] = {fd0, fd1};
        std::array<V_ind, 3> vv[2] = {vert_around(fd[0]->backptr->m_obj->m_mesh, fd[0]->id), vert_around(fd[1]->backptr->m_obj->m_mesh, fd[1]->id)};
        // double thickness = (fd[0]->backptr->m_thickness + fd[1]->backptr->m_thickness)/2;
        // double sqr_thickness = thickness * thickness;

        for (int k = 0; k < 2; ++k)
            for (int p = 0; p < 3; ++p){
                ++ncols.nVFBroadCases;
                bool is_close = np.InitVF(fd, vv, k, p);
                if (!is_close) continue;
                handler.InitVF(np, this);
                if (handler.applyVF())
                    ++ncols.nVFNarrowCases;
                CHECKCOLIFBEGIN()
                    //TESTVFDIST
                    //assert(is_removed && "VF collision is not removed");
                CHECKCOLIFEND();
            }

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j){
                ++ncols.nEEBroadCases;
                bool is_close = np.InitEE(fd, vv, i, j);
                if (!is_close) continue;
                handler.InitEE(np, this);
                if (handler.applyEE())
                    ++ncols.nEENarrowCases;
                // CHECKCOLIF(assert(!np.InitEE(fd, vv, i, j) && "EE collision is not removed"));
            }
    });
    return ncols.nVFNarrowCases + ncols.nEENarrowCases;
}

template<typename NarrowPhase, typename Handler, typename Opt>
int BridsonCollisionManager::_performBroadPhase(VolumeTreeBase& a, VolumeTreeBase& b){
    NarrowPhase np;
    Handler handler;
    return _performBroadPhaseExt<NarrowPhase, Handler, Opt>(a, b, np, handler);
}
template<typename NarrowPhase, typename Handler, typename Opt>
int BridsonCollisionManager::_performBroadPhaseExt(VolumeTreeBase& a, VolumeTreeBase& b, NarrowPhase& np, Handler& handler){
    if (&a == &b) return _performBroadPhaseSelfExt<NarrowPhase, Handler, Opt>(a, np, handler);
    else return _performBroadPhaseDiffExt<NarrowPhase, Handler, Opt>(a, b, np, handler);
}

#undef TESTVFDIST

#endif //AVSIM_BRIDSONCOLLISIONMANAGERINL_H