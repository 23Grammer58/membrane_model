#ifndef MESHGEN_BEZIER_H
#define MESHGEN_BEZIER_H

#include "../AVMesh.h"
#include "StretchFigure.h"

namespace World3d {
    template<int N, typename Field = DReal>
    struct BezierCurve{
        Field B[N+1];

        BezierCurve() = default;
        BezierCurve(std::array<Field, N+1> points){ std::copy(points.data(), points.data() + N + 1, B); }
        BezierCurve(const Field* points) { std::copy(points, points + N + 1, B); }
        Field eval(DReal t) const;
        std::pair<BezierCurve<N, Field>, BezierCurve<N, Field>> splitAt(DReal t_split) const;
        BezierCurve<N-1, Field> derivative() const;
        std::array<Field, N+1> polynomial() const;
        template<int K, typename std::enable_if< (K>N) >::type* = nullptr>
        BezierCurve<K, Field> raise_order() const;
        BezierCurve<N, Field>& Add(const BezierCurve<N, Field>& a) { for (int i = 0; i < N + 1; ++i) B[i] += a.B[i]; return *this; }
        BezierCurve<N, Field>& Add(const Field& a) { for (int i = 0; i < N + 1; ++i) B[i] += a; return *this; }
        BezierCurve<N, Field>& Subtract(const BezierCurve<N, Field>& a) { for (int i = 0; i < N + 1; ++i) B[i] -= a.B[i]; return *this; }
        BezierCurve<N, Field>& Subtract(const Field& a) { for (int i = 0; i < N + 1; ++i) B[i] -= a; return *this; }
        BezierCurve<N, Field>& Scal(DReal a) { for (int i = 0; i < N + 1; ++i) B[i] *= a; return *this; }
        BezierCurve<N, Field>& Op_axpy(DReal a, const BezierCurve<N, Field>& x){ for (int i = 0; i < N + 1; ++i) B[i] += a*x.B[i]; return *this; }
    };
    template<int N, int K, typename Field, typename MakeFieldProduct>
    BezierCurve<N+K, decltype(std::declval<MakeFieldProduct>().operator()(std::declval<Field>(), std::declval<Field>()))>
        PerformBezierProduct(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b, MakeFieldProduct prod = MakeFieldProduct());

    template<int N, int K, typename Field = DReal>
    BezierCurve<N+K, decltype(std::declval<Field>()*std::declval<Field>())> operator*(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b);
    template<int N, typename Field = DReal>
    DReal ArcLength(const BezierCurve<N, Vector>& x, DReal t_st = 0, DReal t_end = 1.0, DReal err = 1e-7);

    template<int N, typename Field = DReal> BezierCurve<N, Field> operator+(const BezierCurve<N, Field>& a, const Field& b){ return BezierCurve<N, Field>(a).Add(b); }
    template<int N, typename Field = DReal> BezierCurve<N, Field> operator+(const Field& a, const BezierCurve<N, Field>& b){ return b+a; }
    template<int N, typename Field = DReal> BezierCurve<N, Field> operator-(const BezierCurve<N, Field>& a, const Field& b){ return BezierCurve<N, Field>(a).Subtract(b); }
    template<int N, typename Field = DReal> BezierCurve<N, Field> operator-(const Field& a, BezierCurve<N, Field>& b){ std::array<Field, N+1> B; for (int i = 0; i < N + 1; ++i) B[i] = a - b.B[i]; return {B}; }
    template<int N, typename Field = DReal> BezierCurve<N, Field> operator+(const BezierCurve<N, Field>& a, const BezierCurve<N, Field>& b){ return BezierCurve<N, Field>(a).Add(b); }
    template<int N, int K, typename Field = DReal, typename std::enable_if< (N>K) >::type* = nullptr> BezierCurve<N, Field> operator+(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b){ return BezierCurve<N, Field>(a).Add(b.template raise_order<N>()); }
    template<int N, int K, typename Field = DReal, typename std::enable_if< (N<K) >::type* = nullptr> BezierCurve<K, Field> operator+(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b){ return b+a; }
    template<int N, typename Field = DReal> BezierCurve<N, Field> operator-(const BezierCurve<N, Field>& a, const BezierCurve<N, Field>& b){ return BezierCurve<N, Field>(a).Subtract(b); }
    template<int N, int K, typename Field = DReal, typename std::enable_if< (N>K) >::type* = nullptr> BezierCurve<N, Field> operator-(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b){ return BezierCurve<N, Field>(a).Subtract(b.template raise_order<N>()); }
    template<int N, int K, typename Field = DReal, typename std::enable_if< (N<K) >::type* = nullptr> BezierCurve<K, Field> operator-(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b){ return BezierCurve<K, Field>(a.template raise_order<K>()).Subtract(b); }
    template<int N, typename Field1, typename Field2> BezierCurve<N, decltype(declval<Field1>()*declval<Field2>())> operator*(const Field1& a, const BezierCurve<N, Field2>& b);
    template<int N, typename Field1, typename Field2> BezierCurve<N, decltype(declval<Field1>()*declval<Field2>())> operator*(const BezierCurve<N, Field1>& a, const Field2& b);
    template<int N, typename Field = DReal> std::ostream& operator<<(std::ostream& out, const BezierCurve<N, Field>& a);
}

#include "Bezier.inl"

#endif