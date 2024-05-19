#include "DataDriven.h"

using namespace World3d;

std::array<double, 3> DataDrivenMaterial::PK2_tensor(std::array<double, 3> E, F_ind f){
    Interp& interp = *(interps[interpolantChooser(this, _obj, f)]);
    std::array<double, 3> xi;
    Response::fromStrain(xi, E); 
    std::array<double, 3> r = interp.evaluate(xi);
    std::array<std::array<double, 3>, 3> dxi_dE;
    Response::derivByStrain(xi, dxi_dE);
    std::array<double, 3> pk2 = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        pk2[i] = r[0] * dxi_dE[0][i] + r[1] * dxi_dE[1][i] + r[2] * dxi_dE[2][i];
    return pk2;    
}
std::pair<std::array<double, 3>, std::array<double, 6>> DataDrivenMaterial::PK2_dPK2_tensor(std::array<double, 3> E, F_ind f){
    Interp& interp = *(interps[interpolantChooser(this, _obj, f)]);
    std::array<double, 3> xi;
    Response::fromStrain(xi, E); 
    std::array<double, 3> r = interp.evaluate(xi);
    std::array<std::array<double, 3>, 3> dr_ = interp.gradient(xi);
    std::array<std::array<double, 3>, 3> dxi_dE;
    std::array<std::array<double, 6>, 3> ddxi_dE;
    Response::secondDerivByStrain(xi, dxi_dE, ddxi_dE);
    std::array<double, 3> pk2 = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        pk2[i] = r[0] * dxi_dE[0][i] + r[1] * dxi_dE[1][i] + r[2] * dxi_dE[2][i];
    Eigen::Matrix<double, 3, 3> dr, dX;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j){
        dr(i, j) = dr_[i][j];
        dX(i, j) = dxi_dE[i][j];
    }
    std::array<Eigen::Matrix<double, 3, 3>, 3> ddxi;
    for (int k = 0; k < 3; ++k)
        ddxi[k] <<  ddxi_dE[k][0], ddxi_dE[k][3], ddxi_dE[k][4],
                    ddxi_dE[k][3], ddxi_dE[k][1], ddxi_dE[k][5],
                    ddxi_dE[k][4], ddxi_dE[k][5], ddxi_dE[k][1];
    Eigen::Matrix<double, 3, 3> dpk2_dE = dX.transpose() * dr * dX + r[0] * ddxi[0] + r[1] * ddxi[1] + r[2] * ddxi[2];
    
    std::array<double, 6> dpk2 = { dpk2_dE(0, 0), dpk2_dE(1, 1), dpk2_dE(2, 2), dpk2_dE(0, 1), dpk2_dE(0, 2), dpk2_dE(1, 2) };

    return {pk2, dpk2}; 
}