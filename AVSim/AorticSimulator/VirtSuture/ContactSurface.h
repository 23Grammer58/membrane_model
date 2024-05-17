#ifndef VIRTSUTURE_CONTACT_SURFACE_H
#define VIRTSUTURE_CONTACT_SURFACE_H

#include "AVSim/Core/MeshGen/Bezier.h"
#include "AVSim/Core/SignedDistanceField.h"

namespace World3d {
    struct CilindricContactSurface{
        struct LeafContactSDF: public SignedDistanceField{
            DReal m_lambda = 0;
            std::array<DReal, 2> m_sigma;
            Vector m_n;
            std::array<BezierCurve<3, Vector>, 2> m_main_lines;
            std::array<BezierCurve<2, Vector>, 2> m_main_derivatives;
            Vector opposit_point;
            
            std::array<BezierCurve<3, Vector>, 2> m_current_lines;
            std::array<BezierCurve<2, Vector>, 2> m_current_derivatives;
            std::array<BezierCurve<1, Vector>, 2> m_current_dderivs;
            std::array<BezierCurve<5, DReal>, 2> m_current_dist_func0;
            
            void make_shift(const Vector& v);
            void set_lambda(DReal lambda);
            SignedDistanceField::SDF operator()(const Vector& X) const override;
        };
        struct LeafPlaneContactSDF: public SignedDistanceField{
            DReal m_lambda = 0;
            Vector m_n;
            Vector m_A, m_C, m_B, m_Cs;

            void set_lambda(DReal lambda);
            SignedDistanceField::SDF operator()(const Vector& X) const override;
        };
        struct DilationPointRes{
            Vector C;
            DReal v;
            bool on_border;
        };

        std::array<Vector, 3> m_A;
        std::array<DReal, 3> m_l;

        Vector m_n, m_C;
        std::array<Vector, 3> m_B;
        std::array<DReal, 3> m_sigma;
        std::array<BezierCurve<3, Vector>, 3> m_lines;

        CilindricContactSurface(){}
        CilindricContactSurface(const std::array<Vector, 3>& commissure_points, const std::array<DReal, 3>& leaf_width);
        CilindricContactSurface(std::array<Point, 3> commissure_points, const std::array<DReal, 3>& leaf_width);
        void setup();
        void setPointsAndWeights(const std::array<Vector, 3>& commissure_points, const std::array<DReal, 3>& leaf_width);
        std::shared_ptr<LeafContactSDF> getLeafContactSdf(int ileaf);
        std::shared_ptr<LeafPlaneContactSDF> getLeafPlaneContactSdf(int ileaf);
        
        static std::pair<Vector,DReal> computeNormalAndArea(const std::array<Vector, 3>& A);
        static DilationPointRes computeDilationPointFromLengths(const std::array<Point, 3>& A, const std::array<DReal, 3>& l);
        static DilationPointRes computeDilationPoint(const std::array<Vector, 3>& A, const std::array<DReal, 3>& c);
        static std::array<Vector, 3> compute_B_points(const std::array<Vector, 3>& A, const Vector& C);
        static std::pair<BezierCurve<3, Vector>, DReal> compute_line(const Vector& A, const Vector& B, const Vector& C, DReal c, DReal rel_err = 1e-7);
    };
    

}

#endif