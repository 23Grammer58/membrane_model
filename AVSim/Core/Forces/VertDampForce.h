#ifndef AVSIM_VERTEXDAMPFORCE_H
#define AVSIM_VERTEXDAMPFORCE_H
#include <functional>
#include "../Object3D.h"
#include "../World.h"

namespace World3d {
    class WVertDampForce: public WorldForceBase{
    public:
        struct IDampFunc{
            World* m_w = nullptr;

            virtual bool registerWorld(World& w) { m_w = &w; return true; };
            virtual void setCurrentState() = 0;
            virtual double operator()(Object3D& obj, const ObjectID id, V_ind v) = 0;
        };
        struct ObjectData{
            Object3D* m_obj;
            Mesh::Property_map<V_ind, Point> m_from_x;
            
            void resetFromX(){ 
                if (m_obj)
                    for (auto v: m_obj->m_mesh.vertices())
                        m_from_x[v] = m_obj->m_x[v];
            }
        };
        struct DampFuncConst: public IDampFunc{
            using DampFunc = std::function<double(Object3D& obj, const ObjectID id, V_ind v)>;
            DampFunc m_damp;

            DampFuncConst(DampFunc damp): m_damp{std::move(damp)} {}
            DampFuncConst& setDampFunc(DampFunc damp) { m_damp = std::move(damp); return *this; }
            DampFuncConst& setDampFunc(DReal damp_coef) { m_damp = [damp_coef](Object3D& obj, const ObjectID id, V_ind v) { return damp_coef; }; return *this; }
            void setCurrentState() override {};
            double operator()(Object3D& obj, const ObjectID id, V_ind v) override {
                auto sigma_coef = m_damp(obj, id, v);
                if (sigma_coef == 0) return 0;
                auto _ms = obj.m_mesh.property_map<V_ind, DReal>("v:mass");
                if (!_ms.second) throw std::runtime_error("On object should be defined \"v:mass\" tag");
                auto ms = _ms.first;
                return sigma_coef * ms[v];
            };
        };
        struct DampFuncEvalDist: public IDampFunc{
            double max_part_ignore = 0.1;
            long max_cnt_ignore = std::numeric_limits<long>::max();
            double max_eval_shift;
            double drop_sigma = 0, max_sigma = std::numeric_limits<double>::max();
            double m_sigma_coef = 0;
            std::vector<double> m_loc_coefs;

            DampFuncEvalDist(double max_eval_shift, double max_part_ignore = 0.1, long max_cnt_ignore = std::numeric_limits<long>::max()):
                max_eval_shift{max_eval_shift}, max_part_ignore{max_part_ignore}, max_cnt_ignore{max_cnt_ignore} {}
            DampFuncEvalDist& setDropSigmaCoef(double _drop_sigma){ drop_sigma = _drop_sigma; return *this; }
            DampFuncEvalDist& setMaxSigmaCoef(double _max_sigma){ max_sigma = _max_sigma; return *this; }
            void setCurrentState() override {
                m_sigma_coef = 0;
                auto Ns = m_w->compStateSize();
                m_w->UpdateForces();
                SparseMatrix sm(Ns);
                m_w->compJacobian(&sm);
                m_loc_coefs.resize(Ns);
                if (m_loc_coefs.empty()) return;
                int N = 0;
                for (auto& o: m_w->objs()){
                    auto& obj = o.second.first;
                    auto _ms = obj.m_mesh.property_map<V_ind, DReal>("v:mass");
                    if (!_ms.second) throw std::runtime_error("On object should be defined \"v:mass\" tag");
                    auto ms = _ms.first;
                    for (auto v: obj.m_mesh.vertices()){
                        auto m = ms[v];
                        for (int i = 0; i < 3; ++i){
                            auto& row = sm[N + v*3 + i];
                            auto it = row.find( N + v*3 + i);
                            DReal val = 0;
                            if (it != row.end()) val = it->second;
                            auto loc_sigma = std::max(abs(obj.m_F[v][i]) / max_eval_shift + val, 0.0);
                            auto loc_coef = loc_sigma / m;
                            m_loc_coefs[N + 3*v + i] = loc_coef;
                        }
                    }
                    N += 3 * obj.m_mesh.num_vertices();
                }
                auto npart = static_cast<long>(N * max_part_ignore);
                if (npart > max_cnt_ignore) npart = max_cnt_ignore;
                if (npart == 0) {
                    m_sigma_coef = *std::max_element(m_loc_coefs.begin(), m_loc_coefs.end());
                } else {
                    std::sort(m_loc_coefs.begin(), m_loc_coefs.end());
                    m_sigma_coef = m_loc_coefs[m_loc_coefs.size() - 1 - npart];
                }
                m_sigma_coef = (m_sigma_coef > drop_sigma) ? m_sigma_coef : 0;
                if (m_sigma_coef > max_sigma) m_sigma_coef = max_sigma;
                std::cout << "Sigma = " << m_sigma_coef << std::endl; //TODO: remove this line
            }
            double operator()(Object3D& obj, const ObjectID id, V_ind v) override{
                if (m_sigma_coef == 0) return 0;
                auto _ms = obj.m_mesh.property_map<V_ind, DReal>("v:mass");
                if (!_ms.second) throw std::runtime_error("On object should be defined \"v:mass\" tag");
                auto ms = _ms.first;
                return m_sigma_coef * ms[v];
            }
        };

        World* m_w = nullptr;
        std::map<ObjectID, ObjectData> m_objs;
        std::shared_ptr<IDampFunc> m_damp;
        
        WVertDampForce(){ type = "WVertDampForce"; }
        WVertDampForce(const WVertDampForce& a) = default;
        WVertDampForce& setDampFunc(std::shared_ptr<IDampFunc> f){ m_damp = std::move(f); return *this; }
        bool addObj(Object3D& obj, const ObjectID& id, bool is_dynamic) override {
            if (is_dynamic){
                ObjectData od;
                od.m_obj = &obj;
                auto r = obj.m_mesh.add_property_map<V_ind, Point>("v:from_x");
                od.m_from_x = r.first;
                if (r.second) od.resetFromX();
                m_objs.insert({id, std::move(od)});
                return true;
            }
            return false;
        }
        void resetCurrentState(){ 
            for (auto& o: m_objs)
                o.second.resetFromX();
            m_damp->setCurrentState();
        }
        void removeObj(const ObjectID& id) override { m_objs.erase(id); }
        bool registerWorld(World& w){ 
            m_w = &w;
            m_objs.clear();
            for (auto& o: m_w->objs())
                addObj(o.second.first, o.first, o.second.second > 0);
            m_damp->registerWorld(w);
            return true;
        }
        int operator()() override{
            for (auto& o: m_objs){
                auto& od = o.second;
                for (auto v: od.m_obj->m_mesh.vertices()){
                    auto sigma = (*m_damp)(*od.m_obj, o.first, v);
                    if (sigma != 0) od.m_obj->m_F[v] += od.m_obj->withBCmask(v, -sigma * (od.m_obj->m_x[v] - od.m_from_x[v])); 
                }
            }
            return 0;
        }
        int fill_matrix(SparseMatrix *sm){
            int N = 0;
            for (auto& o: m_objs){
                auto& od = o.second;
                auto& obj = *od.m_obj;
                for (auto v: obj.m_mesh.vertices()){
                    auto sigma = (*m_damp)(obj, o.first, v);
                    std::array<bool, 3> bc = {obj.is_movable(v, 0), obj.is_movable(v, 1), obj.is_movable(v, 2) };
                    for (int i = 0; i < 3; ++i)
                        if (bc[i] && sigma != 0) (*sm)(N + v*3 + i, N + v*3 + i) -= sigma;
                }
                N += 3 * obj.m_mesh.num_vertices();
            }
            return 0;
        }
        std::unique_ptr<WorldForceBase> clone() override { return std::make_unique<WVertDampForce>(*this); }
    };

    // class VertDampForce: public ForceBase{
    //     public:
    //     struct IDampFunc{
    //         virtual void updateObjData(Object3D* obj) = 0;
    //         virtual double operator()(V_ind v) = 0;
    //     };
    //     using DampFunc = std::function<double(Object3D*, V_ind)>;
    //     struct DampFuncNoUpdate: public IDampFunc{
    //         DampFunc m_damp;
    //         Object3D* m_obj;

    //         DampFuncNoUpdate(DampFunc damp): m_damp{damp} {};
    //         void updateObjData(Object3D* obj) { m_obj = obj; }
    //         double operator()(V_ind v){ return m_damp(m_obj, v); }
    //     };
    //     struct DampFuncEvalDist:  public IDampFunc{
    //         Object3D* m_obj;
    //         std::vector<double> m_F_sqr_norms;
    //         double max_part_ignore = 0.1;
    //         long max_cnt_ignore = std::numeric_limits<long>::max();
    //         double max_eval_shift;
    //         double drop_sigma = 0;
    //         double m_sigma;

    //         DampFuncEvalDist(double max_eval_shift, double max_part_ignore = 0.1, long max_cnt_ignore = std::numeric_limits<long>::max()):
    //             max_eval_shift{max_eval_shift}, max_part_ignore{max_part_ignore}, max_cnt_ignore{max_cnt_ignore} {}
    //         DampFuncEvalDist& setDropSigma(double _drop_sigma){ drop_sigma = _drop_sigma; return *this; }
    //         void updateObjData(Object3D* obj) override{
    //             m_obj = obj;
    //             auto nvert = m_obj->m_mesh.num_vertices();
    //             auto npart = static_cast<long>(nvert * max_part_ignore);
    //             if (npart > max_cnt_ignore) npart = max_cnt_ignore;
    //             double F = std::numeric_limits<double>::max();
    //             if (npart > 0){
    //                 m_F_sqr_norms.resize(nvert);
    //                 int i = 0;
    //                 for (auto v: m_obj->m_mesh.vertices())
    //                     m_F_sqr_norms[i++] = m_obj->m_F[v].squared_length();
    //                 std::sort(m_F_sqr_norms.begin(), m_F_sqr_norms.end());
    //                 F = sqrt(m_F_sqr_norms[npart]);
    //             } else {
    //                 for (auto v: m_obj->m_mesh.vertices()){
    //                     auto newF =  m_obj->m_F[v].squared_length();
    //                     if (newF < F) F = newF;
    //                 }
    //                 F = sqrt(F);
    //             }
    //             m_sigma = F / max_eval_shift;
    //             if (m_sigma < drop_sigma) m_sigma = 0;
    //         }
    //         double operator()(V_ind v) override { return m_sigma; }
    //     };
    //     Object3D* m_obj;
    //     std::shared_ptr<IDampFunc> m_damp;
    //     Mesh::Property_map<V_ind, Point> m_from_x;

    //     VertDampForce(){ type = "VertDampForce"; }
    //     VertDampForce(DampFunc f): m_damp(new DampFuncNoUpdate(std::move(f))) { type = "VertDampForce"; } 
    //     VertDampForce(const VertDampForce& a) = default;
    //     VertDampForce& setDampFunc(DampFunc f){ m_damp = std::shared_ptr<IDampFunc>(new DampFuncNoUpdate(std::move(f))); return *this; }
    //     VertDampForce& setDampFunc(std::shared_ptr<IDampFunc> f){ m_damp = std::move(f); return *this; }
    //     void resetFromX(){ 
    //         if (m_obj)
    //             for (auto v: m_obj->m_mesh.vertices())
    //                m_from_x[v] = m_obj->m_x[v];
    //     }
    //     void registerObj(Object3D* obj){
    //         m_obj = obj;
    //         auto r = obj->m_mesh.add_property_map<V_ind, Point>("v:from_x");
    //         m_from_x = r.first;
    //         if (r.second) resetFromX();
    //     }
    //     int operator()(Object3D &obj) override{
    //         if (&obj != m_obj) registerObj(&obj);
    //         m_damp->updateObjData(m_obj);
    //         for (auto v: m_obj->m_mesh.vertices()){
    //             auto sigma = (*m_damp)(v);
    //             if (sigma != 0) m_obj->m_F[v] += m_obj->withBCmask(v, -sigma * (m_obj->m_x[v] - m_from_x[v])); 
    //         }
    //         return 0;
    //     }
    //     std::unique_ptr<ForceBase> clone() override{ return std::make_unique<VertDampForce>(*this); }
    //     int fill_matrix(Object3D &obj, SparseMatrix::IndexFrom &m) override{
    //         for (auto v: m_obj->m_mesh.vertices()){
    //             auto sigma = (*m_damp)(v);
    //             std::array<bool, 3> bc = { obj.is_movable(v, 0), obj.is_movable(v, 1), obj.is_movable(v, 2) };
    //             for (int i = 0; i < 3; ++i)
    //                 if (bc[i] && sigma != 0) m(v*3 + i, v*3 + i) -= sigma;
    //         }
    //         return 0;
    //     }
    // };
};

#endif //AVSIM_VERTEXDAMPFORCE_H