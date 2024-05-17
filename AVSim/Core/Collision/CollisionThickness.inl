
namespace World3d{
    std::unique_ptr<VolumeTreeBase> _CollisionThicknessCompressed_internal_get_volume_tree();

    template<class ThicknessT, typename std::enable_if<std::is_base_of<CollisionThicknessBase, ThicknessT>::value>::type* U>
    void CollisionThicknessCompressed<ThicknessT, U>::setupIndex(int type0, int type1){
        if (type0 < CTB::VERTEX || type0 > CTB::FACE) throw std::runtime_error("Faced unknown primitive type = " + std::to_string(type0));
        if (type1 < CTB::VERTEX || type1 > CTB::FACE) throw std::runtime_error("Faced unknown primitive type = " + std::to_string(type1));
        if (!CTB::m_obj) return;
        auto _t_x = CTB::m_obj->m_mesh.property_map<V_ind, Point>(m_init_x_tag_name);
        if (!_t_x.second) throw std::runtime_error("Access to uninitialized mesh-field = \"" + m_init_x_tag_name + "\"");
        auto& m_x0 = _t_x.first;
        if (type0 > type1) std::swap(type0, type1);
        #define TP(t0, t1) (t0 * 3 + t1 - t0 * (t0 + 1) / 2)
        int tp = TP(type0, type1);
        if (m_distance_map_matrix[tp] != nullptr) return;
        union ComplIDX{
            unsigned int id;
            void* ptr;
        };
        auto fill_bvt = [&](int t, VolumeTreeBase& bvt){
            auto& o = *CTB::m_obj;
            switch(t){
                case CTB::VERTEX:{
                    for (auto v: o.m_mesh.vertices()){
                        auto th = max_compression*operator()(v.idx(), CTB::VERTEX) / 2;
                        auto bv = BoxVol::FromPoints({m_x0[v]}).Expand(th);
                        ComplIDX idx; idx.id = v.idx();
                        bvt.insert(bv, idx.ptr);    
                    }
                    break;
                }
                case CTB::EDGE:{
                    for (auto v: o.m_mesh.edges()){
                        auto vv = vert_around(o.m_mesh, v);
                        auto th = max_compression*operator()(v.idx(), CTB::EDGE) / 2;
                        auto bv = BoxVol::FromPoints({m_x0[vv[0]], m_x0[vv[1]]}).Expand(th);
                        ComplIDX idx; idx.id = v.idx();
                        bvt.insert(bv, idx.ptr);    
                    }
                    break;
                }
                case CTB::FACE:{
                    for (auto v: o.m_mesh.faces()){
                        auto vv = vert_around(o.m_mesh, v);
                        auto th = max_compression*operator()(v.idx(), CTB::FACE) / 2;
                        auto bv = BoxVol::FromPoints({m_x0[vv[0]], m_x0[vv[1]], m_x0[vv[2]]}).Expand(th);
                        ComplIDX idx; idx.id = v.idx();
                        bvt.insert(bv, idx.ptr);    
                    }
                    break; 
                }
            }
        };

        auto bvt0_ = _CollisionThicknessCompressed_internal_get_volume_tree();
        std::unique_ptr<VolumeTreeBase> bvt1_;
        fill_bvt(type0, *bvt0_); 
        bvt0_->optimize();
        VolumeTreeBase& bvt0 = *bvt0_;
        VolumeTreeBase* bvt1p = &bvt0;
        if (type0 != type1) {
            bvt1_ = _CollisionThicknessCompressed_internal_get_volume_tree();
            fill_bvt(type1, *bvt1_);
            bvt1_->optimize();
            bvt1p = bvt1_.get();
        }
        VolumeTreeBase& bvt1 = *bvt1p;
        auto clear_tag = [](auto t){
            if (t.second) return;
            throw std::runtime_error("Object had nonzero v:close_parts_* tag");
            for (auto& i: t.first)
                i.clear();
        };
        using GP = GeomProjs::Proximity;
        switch(tp){
            case TP(CTB::VERTEX, CTB::VERTEX):{
                auto t_ = CTB::m_obj->m_mesh.add_property_map<V_ind, std::map<V_ind, DReal>>("v:close_parts_" + std::to_string(tp));
                clear_tag(t_);
                VolumeIntersect(bvt0, bvt1, [&t = t_.first, &o = *CTB::m_obj, &m_x0 = _t_x.first, this](const auto l1, void* data1, const auto l2, void* data2){
                    ComplIDX id1, id2;
                    id1.ptr = data1, id2.ptr = data2;
                    V_ind v0(id1.id), v1(id2.id);
                    if (v0 > v1) std::swap(v0, v1);
                    DReal max_dist = max_compression*(this->operator()(v0, CTB::VERTEX) + this->operator()(v1, CTB::VERTEX))/2;
                    DReal max_dist_2 = max_dist * max_dist;
                    auto res = GP::Query(GP::Primitive<Vector>(GP::POINT, {m_x0[v0] - CGAL::ORIGIN}), GP::Primitive<Vector>(GP::POINT, {m_x0[v1] - CGAL::ORIGIN}), max_dist_2);
                    if (res.status && res.sqd > 0){
                        t[v0][v1] = sqrt(res.sqd)/max_compression;
                    }
                });
                m_distance_map_matrix[tp] = reinterpret_cast<const void*>(t_.first.data());
                break;
            }
            case TP(CTB::VERTEX, CTB::EDGE):{
                auto t_ = CTB::m_obj->m_mesh.add_property_map<V_ind, std::map<E_ind, DReal>>("v:close_parts_" + std::to_string(tp));
                clear_tag(t_);
                VolumeIntersect(bvt0, bvt1, [&t = t_.first, &o = *CTB::m_obj, &m_x0 = _t_x.first, this](const auto l1, void* data1, const auto l2, void* data2){
                    ComplIDX id1, id2;
                    id1.ptr = data1, id2.ptr = data2;
                    V_ind v0(id1.id); E_ind e1(id2.id);
                    auto vv = vert_around(o.m_mesh, e1);
                    DReal max_dist = max_compression*(this->operator()(v0, CTB::VERTEX) + this->operator()(e1, CTB::EDGE))/2;
                    DReal max_dist_2 = max_dist * max_dist;
                    auto res = GP::Query(GP::Primitive<Vector>(GP::POINT, {m_x0[v0] - CGAL::ORIGIN}), GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv[0]] - CGAL::ORIGIN, m_x0[vv[1]] - CGAL::ORIGIN}), max_dist_2);
                    if (res.status && res.sqd > 0){
                        t[v0][e1] = sqrt(res.sqd)/max_compression;
                    }
                });
                m_distance_map_matrix[tp] = reinterpret_cast<const void*>(t_.first.data());
                break;
            }
            case TP(CTB::VERTEX, CTB::FACE):{
                auto t_ = CTB::m_obj->m_mesh.add_property_map<V_ind, std::map<F_ind, DReal>>("v:close_parts_" + std::to_string(tp));
                clear_tag(t_);
                VolumeIntersect(bvt0, bvt1, [&t = t_.first, &o = *CTB::m_obj, &m_x0 = _t_x.first, this](const auto l1, void* data1, const auto l2, void* data2){
                    ComplIDX id1, id2;
                    id1.ptr = data1, id2.ptr = data2;
                    V_ind v0(id1.id); F_ind f1(id2.id);
                    auto vv = vert_around(o.m_mesh, f1);
                    DReal max_dist = max_compression*(this->operator()(v0, CTB::VERTEX) + this->operator()(f1, CTB::FACE))/2;
                    DReal max_dist_2 = max_dist * max_dist;
                    auto res = GP::Query(GP::Primitive<Vector>(GP::POINT, {m_x0[v0] - CGAL::ORIGIN}), GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv[0]] - CGAL::ORIGIN, m_x0[vv[1]] - CGAL::ORIGIN, m_x0[vv[2]] - CGAL::ORIGIN}), max_dist_2);
                    if (res.status && res.sqd > 0){
                        t[v0][f1] = sqrt(res.sqd)/max_compression;
                    }
                });

                m_distance_map_matrix[tp] = reinterpret_cast<const void*>(t_.first.data());
                break;
            }
            case TP(CTB::EDGE, CTB::EDGE):{
                auto t_ = CTB::m_obj->m_mesh.add_property_map<E_ind, std::map<E_ind, DReal>>("e:close_parts_" + std::to_string(tp));
                clear_tag(t_);
                VolumeIntersect(bvt0, bvt1, [&t = t_.first, &o = *CTB::m_obj, &m_x0 = _t_x.first, this](const auto l1, void* data1, const auto l2, void* data2){
                    ComplIDX id1, id2;
                    id1.ptr = data1, id2.ptr = data2;
                    E_ind e0(id1.id); E_ind e1(id2.id);
                    if (e0 > e1) std::swap(e0, e1);
                    auto vv0 = vert_around(o.m_mesh, e0), vv1 = vert_around(o.m_mesh, e1);
                    DReal max_dist = max_compression*(this->operator()(e0, CTB::EDGE) + this->operator()(e1, CTB::EDGE))/2;
                    DReal max_dist_2 = max_dist * max_dist;
                    auto res = GP::Query(GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv0[0]] - CGAL::ORIGIN, m_x0[vv0[1]] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN}), max_dist_2);
                    if (res.status && res.sqd > 0){
                        t[e0][e1] = sqrt(res.sqd)/max_compression;
                    }
                });
                m_distance_map_matrix[tp] = reinterpret_cast<const void*>(t_.first.data());
                break;
            }
            case TP(CTB::EDGE, CTB::FACE):{
                auto t_ = CTB::m_obj->m_mesh.add_property_map<F_ind, std::map<E_ind, DReal>>("f:close_parts_" + std::to_string(tp));
                clear_tag(t_);
                VolumeIntersect(bvt0, bvt1, [&t = t_.first, &o = *CTB::m_obj, &m_x0 = _t_x.first, this](const auto l1, void* data1, const auto l2, void* data2){
                    ComplIDX id1, id2;
                    id1.ptr = data1, id2.ptr = data2;
                    E_ind e0(id1.id); F_ind f1(id2.id);
                    auto vv0 = vert_around(o.m_mesh, e0); auto vv1 = vert_around(o.m_mesh, f1);
                    DReal max_dist = max_compression*(this->operator()(e0, CTB::EDGE) + this->operator()(f1, CTB::FACE))/2;
                    DReal max_dist_2 = max_dist * max_dist;
                    auto res = GP::Query(GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv0[0]] - CGAL::ORIGIN, m_x0[vv0[1]] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN, m_x0[vv1[2]] - CGAL::ORIGIN}), max_dist_2);
                    if (res.status && res.sqd > 0){
                        t[f1][e0] = sqrt(res.sqd)/max_compression;
                    }
                });
                m_distance_map_matrix[tp] = reinterpret_cast<const void*>(t_.first.data());
                break;
            }
            case TP(CTB::FACE, CTB::FACE):{
                auto t_ = CTB::m_obj->m_mesh.add_property_map<F_ind, std::map<F_ind, DReal>>("f:close_parts_" + std::to_string(tp));
                clear_tag(t_);
                VolumeIntersect(bvt0, bvt1, [&t = t_.first, &o = *CTB::m_obj, &m_x0 = _t_x.first, this](const auto l1, void* data1, const auto l2, void* data2){
                    ComplIDX id1, id2;
                    id1.ptr = data1, id2.ptr = data2;
                    F_ind f0(id1.id); F_ind f1(id2.id);
                    if (f0 > f1) std::swap(f0, f1);
                    auto vv0 = vert_around(o.m_mesh, f0), vv1 = vert_around(o.m_mesh, f1);
                    DReal max_dist = max_compression*(this->operator()(f0, CTB::FACE) + this->operator()(f1, CTB::FACE))/2;
                    DReal max_dist_2 = max_dist * max_dist;
                    auto res = GP::Query(GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv0[0]] - CGAL::ORIGIN, m_x0[vv0[1]] - CGAL::ORIGIN, m_x0[vv0[2]] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN, m_x0[vv1[2]] - CGAL::ORIGIN}), max_dist_2);
                    if (res.status && res.sqd > 0){
                        t[f0][f1] = sqrt(res.sqd)/max_compression;
                    }
                });
                m_distance_map_matrix[tp] = reinterpret_cast<const void*>(t_.first.data());
                break;
            }
        }
        #undef TP
    }
    template<class ThicknessT, typename std::enable_if<std::is_base_of<CollisionThicknessBase, ThicknessT>::value>::type* U>
    void CollisionThicknessCompressed<ThicknessT, U>::removeObj(){
        if (!CTB::m_obj) return;
        #define TP(t0, t1) (t0 * 3 + t1 - t0 * (t0 + 1) / 2)
        for (int t0 = 0; t0 < 3; ++t0)
        for (int t1 = t0; t1 < 3; ++t1){
            auto tp = TP(t0, t1);
            if (m_distance_map_matrix[tp] != nullptr){
                switch (tp){
                    case TP(CTB::VERTEX, CTB::VERTEX):{
                        auto t_ = CTB::m_obj->m_mesh.property_map<V_ind, std::map<V_ind, DReal>>("v:close_parts_" + std::to_string(tp));
                        if (t_.second) CTB::m_obj->m_mesh.remove_property_map(t_.first);
                        break;
                    }
                    case TP(CTB::VERTEX, CTB::EDGE):{
                        auto t_ = CTB::m_obj->m_mesh.property_map<V_ind, std::map<E_ind, DReal>>("v:close_parts_" + std::to_string(tp));
                        if (t_.second) CTB::m_obj->m_mesh.remove_property_map(t_.first);
                        break;
                    }
                    case TP(CTB::VERTEX, CTB::FACE):{
                        auto t_ = CTB::m_obj->m_mesh.property_map<V_ind, std::map<F_ind, DReal>>("v:close_parts_" + std::to_string(tp));
                        if (t_.second) CTB::m_obj->m_mesh.remove_property_map(t_.first);
                        break;
                    }
                    case TP(CTB::EDGE, CTB::EDGE):{
                        auto t_ = CTB::m_obj->m_mesh.property_map<E_ind, std::map<E_ind, DReal>>("e:close_parts_" + std::to_string(tp));
                        if (t_.second) CTB::m_obj->m_mesh.remove_property_map(t_.first);
                        break;
                    }
                    case TP(CTB::EDGE, CTB::FACE):{
                        auto t_ = CTB::m_obj->m_mesh.property_map<F_ind, std::map<E_ind, DReal>>("f:close_parts_" + std::to_string(tp));
                        if (t_.second) CTB::m_obj->m_mesh.remove_property_map(t_.first);
                        break;
                    }
                    case TP(CTB::FACE, CTB::FACE):{
                        auto t_ = CTB::m_obj->m_mesh.property_map<F_ind, std::map<F_ind, DReal>>("f:close_parts_" + std::to_string(tp));
                        if (t_.second) CTB::m_obj->m_mesh.remove_property_map(t_.first);
                        break;
                    }
                }
            }
            m_distance_map_matrix[tp] = nullptr;
        }
        #undef TP
        ThicknessT::removeObj();
    }
    template<class ThicknessT, typename std::enable_if<std::is_base_of<CollisionThicknessBase, ThicknessT>::value>::type* U>
    void CollisionThicknessCompressed<ThicknessT, U>::registerObj(Object3D* obj){
        if (CTB::m_obj != obj) removeObj();
        ThicknessT::registerObj(obj);
    }
    template<class ThicknessT, typename std::enable_if<std::is_base_of<CollisionThicknessBase, ThicknessT>::value>::type* U>
    DReal CollisionThicknessCompressed<ThicknessT, U>::operator()(int id0, int id1, int type0, int type1){
        if (type0 < CTB::VERTEX || type0 > CTB::FACE) throw std::runtime_error("Faced unknown primitive type = " + std::to_string(type0));
        if (type1 < CTB::VERTEX || type1 > CTB::FACE) throw std::runtime_error("Faced unknown primitive type = " + std::to_string(type1));
        if (!CTB::m_obj) throw std::runtime_error("Thickness query on empty object");
        if (type0 > type1) std::swap(id0, id1), std::swap(type0, type1);

        #define TP(t0, t1) (t0 * 3 + t1 - t0 * (t0 + 1) / 2)
        auto tp = TP(type0, type1);
        if (!m_distance_map_matrix[tp]){
            DReal th_base = (operator()(id0, type0) + operator()(id1, type1))/2; 
            DReal init_dist = th_base*max_compression;
            DReal th_base_2 = init_dist*init_dist;
            
            auto _t_x = CTB::m_obj->m_mesh.property_map<V_ind, Point>(m_init_x_tag_name);
            if (!_t_x.second) throw std::runtime_error("Access to uninitialized mesh-field = \"" + m_init_x_tag_name + "\"");
            auto& m_x0 = _t_x.first;
            
            using GP = GeomProjs::Proximity;
            auto& o = *CTB::m_obj;
            switch (tp){
                case TP(CTB::VERTEX, CTB::VERTEX):{
                    V_ind v0(id0); V_ind v1(id1);
                    auto res = GP::Query(GP::Primitive<Vector>(GP::POINT, {m_x0[v0] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::POINT, {m_x0[v1] - CGAL::ORIGIN}), th_base_2);
                    if (res.status) init_dist = sqrt(res.sqd)/max_compression;
                    break;
                }
                case TP(CTB::VERTEX, CTB::EDGE):{
                    V_ind v0(id0); E_ind e1(id1);
                    auto vv1 = vert_around(o.m_mesh, e1);
                    auto res = GP::Query(GP::Primitive<Vector>(GP::POINT, {m_x0[v0] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN}), th_base_2);
                    if (res.status) init_dist = sqrt(res.sqd)/max_compression;
                    break;
                }
                case TP(CTB::VERTEX, CTB::FACE):{
                    V_ind v0(id0); F_ind f1(id1);
                    auto vv1 = vert_around(o.m_mesh, f1);
                    auto res = GP::Query(GP::Primitive<Vector>(GP::POINT, {m_x0[v0] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN, m_x0[vv1[2]] - CGAL::ORIGIN}), th_base_2);
                    if (res.status) init_dist = sqrt(res.sqd)/max_compression;
                    break;
                }
                case TP(CTB::EDGE, CTB::EDGE):{
                    E_ind e0(id0); E_ind e1(id1);
                    auto vv0 = vert_around(o.m_mesh, e0), vv1 = vert_around(o.m_mesh, e1);
                    auto res = GP::Query(GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv0[0]] - CGAL::ORIGIN, m_x0[vv0[1]] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN}), th_base_2);
                    if (res.status) init_dist = sqrt(res.sqd)/max_compression;
                    break;
                }
                case TP(CTB::EDGE, CTB::FACE):{
                    E_ind e0(id0); E_ind f1(id1);
                    auto vv0 = vert_around(o.m_mesh, e0); auto vv1 = vert_around(o.m_mesh, f1);
                    auto res = GP::Query(GP::Primitive<Vector>(GP::SEGMENT, {m_x0[vv0[0]] - CGAL::ORIGIN, m_x0[vv0[1]] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN, m_x0[vv1[2]] - CGAL::ORIGIN}), th_base_2);
                    if (res.status) init_dist = sqrt(res.sqd)/max_compression;
                    break;
                }
                case TP(CTB::FACE, CTB::FACE):{
                    F_ind e0(id0); F_ind e1(id1);
                    auto vv0 = vert_around(o.m_mesh, e0), vv1 = vert_around(o.m_mesh, e1);
                    auto res = GP::Query(GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv0[0]] - CGAL::ORIGIN, m_x0[vv0[1]] - CGAL::ORIGIN, m_x0[vv0[2]] - CGAL::ORIGIN}), 
                                        GP::Primitive<Vector>(GP::TRIANGLE, {m_x0[vv1[0]] - CGAL::ORIGIN, m_x0[vv1[1]] - CGAL::ORIGIN, m_x0[vv1[2]] - CGAL::ORIGIN}), th_base_2);
                    if (res.status) init_dist = sqrt(res.sqd)/max_compression;
                    break;
                }
            }
            return std::min(th_base, init_dist);
        }
        if (type0 == type1 && id0 > id1) std::swap(id0, id1);
        switch (tp)
        {
            case TP(CTB::VERTEX, CTB::VERTEX):{
                auto dat = reinterpret_cast<const std::map<V_ind, DReal>*>(m_distance_map_matrix[tp]);
                
                auto it = dat[id0].find(V_ind(id1));
                if (it != dat[id0].end()) return it->second;
                else return (operator()(id0, type0) + operator()(id1, type1))/2;
                break;
            }
            case TP(CTB::VERTEX, CTB::EDGE):{
                auto dat = reinterpret_cast<const std::map<E_ind, DReal>*>(m_distance_map_matrix[tp]);
                auto it = dat[id0].find(E_ind(id1));
                if (it != dat[id0].end()) return it->second;
                else return (operator()(id0, type0) + operator()(id1, type1))/2;
                break;
            }
            case TP(CTB::VERTEX, CTB::FACE):{
                auto dat = reinterpret_cast<const std::map<F_ind, DReal>*>(m_distance_map_matrix[tp]);
                auto it = dat[id0].find(F_ind(id1));
                if (it != dat[id0].end()) return it->second;
                else return (operator()(id0, type0) + operator()(id1, type1))/2;
                break;
            }
            case TP(CTB::EDGE, CTB::EDGE):{
                auto dat = reinterpret_cast<const std::map<E_ind, DReal>*>(m_distance_map_matrix[tp]);
                auto it = dat[id0].find(E_ind(id1));
                if (it != dat[id0].end()) return it->second;
                else return (operator()(id0, type0) + operator()(id1, type1))/2;
                break;
            }
            case TP(CTB::EDGE, CTB::FACE):{
                auto dat = reinterpret_cast<const std::map<E_ind, DReal>*>(m_distance_map_matrix[tp]);
                auto it = dat[id0].find(E_ind(id1));
                if (it != dat[id0].end()) return it->second;
                else return (operator()(id0, type0) + operator()(id1, type1))/2;
                break;
            }
            case TP(CTB::FACE, CTB::FACE):{
                auto dat = reinterpret_cast<const std::map<F_ind, DReal>*>(m_distance_map_matrix[tp]);
                auto it = dat[id0].find(F_ind(id1));
                if (it != dat[id0].end()) return it->second;
                else return (operator()(id0, type0) + operator()(id1, type1))/2;
                break;
            }
        }
        #undef TP

        return (operator()(id0, type0) + operator()(id1, type1))/2;
    }
}