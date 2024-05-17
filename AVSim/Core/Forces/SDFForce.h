#ifndef AORTIC_VALVE_SDFFORCE_H
#define AORTIC_VALVE_SDFFORCE_H
#include "AVSim/Core/Object3D.h"
#include "AVSim/Core/SignedDistanceField.h"

namespace World3d {
    class SDFForce : public ForceBase {
        UpdatableProperty<F_ind, Vector> _S;
        Mesh::Property_map<V_ind, SignedDistanceField::SDF> _sdf_data;
        std::shared_ptr<SignedDistanceField> m_field;
        Object3D *_obj = nullptr;
        using Func1D = std::function<double(double)>;
        Func1D m_dist_func, m_dist_deriv;
    public:
        SDFForce(const SDFForce&) = default;
        SDFForce(std::shared_ptr<SignedDistanceField> field, DReal P, DReal half_penetration_depth, DReal dist_shift = 0);
        SDFForce(std::shared_ptr<SignedDistanceField> field, Func1D dist_f, Func1D dist_df);
        template<typename ForceT = SignedDistanceField>
        ForceT *field() { 
            if (!m_field || !m_field.get()) return nullptr;
            return static_cast<ForceT *>(m_field.get()); 
        }
        void set_dist_funcs(DReal P, DReal half_penetration_depth, DReal dist_shift = 0);
        void set_dist_funcs(std::function<double(double d)> dist_func, std::function<double(double d)> dist_deriv);
        void registerObj(Object3D *obj);
        int operator()(Object3D &obj) override;
        virtual ~SDFForce(){ clear_sdf_property(); }
        std::unique_ptr<ForceBase> clone() override{ auto res = std::make_unique<SDFForce>(*this); res->_obj = nullptr; return res;}
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override;
        static std::pair<Func1D, Func1D> get_default_dist_func(DReal P, DReal half_penetration_depth, DReal dist_shift = 0);
        static std::pair<Func1D, Func1D> get_default_dist_func(DReal* P, DReal* half_penetration_depth, DReal* dist_shift = nullptr);
    private:
        void create_sdf_property();
        void clear_sdf_property();
        void fill_sdf_data();
    };
}

#endif