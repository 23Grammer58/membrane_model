//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_MSM_H
#define AORTIC_VALVE_MSM_H

#include "../Object3D.h"

namespace World3d {
    class SimpleMassSpringModel : public ForceBase {
        DReal _E;
        DReal _h;
        ConstProperty <E_ind, DReal> _l0;
        UpdatableProperty <E_ind, DReal> _l;
        Object3D *_obj = nullptr;
        typedef DReal SavedData;
        ConstProperty <E_ind, SavedData> _k;
    public:
        SimpleMassSpringModel(double E, double h) : _E{E}, _h{h} { type = "SimpleMassSpringModel"; }

        void registerObj(Object3D *obj);
        int operator()(Object3D &obj) override;
        std::unique_ptr <ForceBase> clone() override;
        int element_matrix(Object3D &obj, LocMatrix &matr, F_ind element) override;

        int element_matrix(Object3D &obj, LocMatrix &matr, E_ind edge);
    };
}


#endif //AORTIC_VALVE_MSM_H
