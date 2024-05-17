//
// Created by alex on 04.05.2021.
//

#ifndef AORTIC_VALVE_STRETCHFIGURE_H
#define AORTIC_VALVE_STRETCHFIGURE_H

#include "ForAni3dFrtPrm.h"
#include "../Object3D.h"

class StretchFigure {
public:
    Vector m_dir = Vector(0, 0, 1);
    std::function<Point(double t)> m_curve = [](double t) { return Point(t, 0, 0); };
    std::function<double(double t)> m_dirLen = [](double t) { return 1.0; };
    std::vector<double> m_tPnts = {0.0, 1.0};
    std::vector<int> m_eLbls = {2, 2, 2, 1};

    StretchFigure() = default;
    StretchFigure(const StretchFigure&) = default;
    StretchFigure(StretchFigure&&) noexcept = default;
    StretchFigure(Vector dir, const std::function<Point(double t)>& curve, const std::function<double(double t)>& dirLen, const std::vector<double>& t_pnts):
        m_dir(dir), m_curve(curve), m_dirLen(dirLen), m_tPnts(t_pnts) {}
    StretchFigure& setStretchDirection(Vector dir) { return m_dir = dir / sqrt(dir.squared_length()), *this; }
    StretchFigure& setCurveFunc(const std::function<Point(double t)>& curve) { return m_curve = curve, *this; }
    StretchFigure& setStretchLengthFunc(const std::function<double(double t)> & len) { return m_dirLen = len, *this; }
    StretchFigure& setPntParams(const std::vector<double>& t_pnts);
    StretchFigure& setPntParams(const std::vector<double>& t_pnts, const std::vector<int>& eLbls);

    [[nodiscard]] Ani3dSurfDiscrWrap prepareFigure() const;
    [[nodiscard]] std::function<Point_2(double, double)> makeParametricFlatter() const;

    [[nodiscard]] Object3D generateMesh(std::function<double(double, double, double)> fsize);
};

struct QuickLineIntegral{
public:
    struct QuickIntegrData{
            double t;
            double l;
            Point p;
    };
    std::vector<QuickIntegrData> quick_index;

    std::function<Point(double)> m_curve;
    double Tst = 0, Tend = 1;
    double rel = 1e-4;
    int Ninit = 100;
    int max_refine_depth = 300;

    double operator()(double t) const;

    void setup();
};

struct PiecewieceLineIntegral{
    std::vector<QuickLineIntegral> pieces;
    double operator()(double t) const;
    void setup();
};

int StretchFigureTest(string dir);


#endif //AORTIC_VALVE_STRETCHFIGURE_H
