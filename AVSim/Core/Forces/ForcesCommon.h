//
// Created by alex on 29.01.2021.
//

#ifndef AORTIC_VALVE_FORCESCOMMON_H
#define AORTIC_VALVE_FORCESCOMMON_H
#include "../Object3D.h"
#include <Eigen/Dense>

namespace World3d{
    struct LocalCartesian{
        Object3D& obj;
        bool flat;
        Eigen::Matrix<DReal, 2, 3> init;
        explicit LocalCartesian(Object3D& obj, bool flat = false);
        Eigen::Matrix<DReal, 2, 3> operator()(F_ind f);
    };

    ConstProperty<F_ind, std::array<DReal, 2 * 3>> set_S2d(Object3D& obj, bool flat, bool forced = false);
    ConstProperty<F_ind, std::array<DReal, 2 * 3>> set_L(Object3D& obj, bool forced = false);
    ConstProperty<F_ind, std::array<DReal, 4 * 2 * 3>> set_N_vecs(Object3D& obj, bool forced = false);
    ConstProperty<F_ind, std::array<DReal, 2 * 2 * 3>> set_JBend(Object3D& obj, bool forced = false);
}


#endif //AORTIC_VALVE_FORCESCOMMON_H
