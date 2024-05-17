#ifndef AORTIC_VALVE_SDF_H
#define AORTIC_VALVE_SDF_H

#include "AVMesh.h"

namespace World3d
{
    class SignedDistanceField{
    public:
        struct SDF{
            DReal sdf;  ///< signed distance
            Vector grad_sdf; ///< gradient of signed distance field
            std::array<DReal, 6> _grad_grad_sdf = {0}; ///< hessian of signed distance field
            DReal& grad_grad_sdf(int i, int j){
                if (i > j) std::swap(i, j);
                int k = j + 2*i - i/2;
                return _grad_grad_sdf[k];
            }
            DReal grad_grad_sdf(int i, int j) const {
                if (i > j) std::swap(i, j);
                int k = j + 2*i - i/2;
                return _grad_grad_sdf[k];
            }
            std::string to_string() const{
                std::stringstream os;
                auto& g = _grad_grad_sdf;
                os << "SDF{ sdf = " << sdf << ", grad_sdf = (" << grad_sdf << ")^T, grad_grad_sdf = ("
                << "[ " << g[0] << ", " << g[1] << ", "  << g[2] << "], " 
                << "[ " << g[1] << ", " << g[3] << ", "  << g[4] << "], " 
                << "[ " << g[2] << ", " << g[4] << ", "  << g[5] << "]) }";
                return os.str(); 
            }
        };

        virtual SDF operator()(const Vector& x) const = 0;
    };
} // namespace World3d


#endif