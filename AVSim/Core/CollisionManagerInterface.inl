//
// Created by Liogky Alexey on 31.05.2022.
//

#ifndef AVSIM_COLLISIONMANAGER_INL
#define AVSIM_COLLISIONMANAGER_INL

namespace World3d{
    template<typename InputIterator>
    void BoxVol::GetBoxFromPP<InputIterator, true>::execute(InputIterator beg, InputIterator end, BoxVol &bv) {
        assert(end != beg && "List of points should not be empty");
        std::array<DReal , 3> st, ed;
        for (int k = 0; k < 3; ++k){
            st[k] = (*beg)[k];
            ed[k] = (*beg)[k];
        }
        for (auto it = ++beg; it < end; ++it){
            const auto& vec = *it;
            for (int k = 0; k < 3; ++k){
                if (st[k] > vec[k]) st[k] = vec[k];
                if (ed[k] < vec[k]) ed[k] = vec[k];
            }
        }
        bv.st = Vector(st[0], st[1], st[2]);
        bv.ed = Vector(ed[0], ed[1], ed[2]);
    }

    inline BoxVol BoxVol::FromCE(const Vector& center, const Vector& extent){
        BoxVol box;
        box.st = center - extent;
        box.ed = center + extent;
        return box;
    }

    template<class PointIterator>
    BoxVol BoxVol::FromPoints(PointIterator first, PointIterator last){
        BoxVol bv;
        using ValType = typename std::iterator_traits<PointIterator>::value_type;
        GetBoxFromPP<PointIterator, std::is_same<ValType, Vector>::value || std::is_same<ValType, Point>::value>::execute(first, last, bv);
        return bv;
    }

    inline BoxVol& BoxVol::SignedExpand(const Vector& e)
    {
        std::array<double, 3> a, b;
        for (int k = 0; k < 3; ++k){
            if (e[k] > 0){
                a[k] = st[k];
                b[k] = ed[k] + e[k];
            } else {
                a[k] = st[k] + e[k];
                b[k] = ed[k];
            }
        }
        st = Vector(a[0], a[1], a[2]), ed = Vector(b[0], b[1], b[2]);
        return *this;
    }

    inline BoxVol& BoxVol::ScaleExpand(const Vector& e){
        Vector c = Center();
        auto l0 = Extents();
        Vector l(l0[0]*e[0], l0[1]*e[1], l0[2]*e[2]);
        st = c - l; ed = c + l;
        return *this;
    }

    inline bool BoxVol::Contain(const BoxVol &a) const {
        return ((st[0] <= a.st[0]) &&
                (st[1] <= a.st[1]) &&
                (st[2] <= a.st[2]) &&
                (ed[0] >= a.ed[0]) &&
                (ed[1] >= a.ed[1]) &&
                (ed[2] >= a.ed[2]));
    }

    inline int BoxVol::Classify(const Vector& n, DReal o, int s) const
    {
        std::array<DReal, 3> pi, px;
        for (int k = 0; k < 3; ++k)
            if (s & (1 << k)){
                pi[k] = ed[k];
                px[k] = st[k];
            } else {
                pi[k] = st[k];
                px[k] = ed[k];
            }
        DReal di = 0, dx = 0;
        for (int k = 0; k < 3; ++k) di += n[k] * px[k], dx = n[k] * pi[k];
        if ((dx + o) < 0) return (-1);
        if ((di + o) >= 0) return (+1);
        return (0);
    }

    inline BoxVol BoxVol::Merge(const BoxVol &a, const BoxVol &b) {
        std::array<DReal, 3> pi, px;
        for (int i = 0; i < 3; ++i)
        {
            if (a.st[i] < b.st[i])
                pi[i] = a.st[i];
            else
                pi[i] = b.st[i];
            if (a.ed[i] > b.ed[i])
                px[i] = a.ed[i];
            else
                px[i] = b.ed[i];
        }
        BoxVol r;
        r.st = Vector(pi[0], pi[1], pi[2]);
        r.ed = Vector(px[0], px[1], px[2]);
        return r;
    }

    namespace GeomProjs {
        template<typename NT>
        static inline NT Clamp_internal(const NT &x, const NT &l, const NT &h) {
            return (x < l ? l : x > h ? h : x);
        }

        template<class Vector3d, typename NT>
        static inline void ProjectOrigin(const Vector3d &a,
                                         const Vector3d &b,
                                         Vector3d &prj,
                                         NT &sqd) {
            const Vector3d d = b - a;
            auto squared_length = [](const Vector3d &d) { return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; };
            const NT m2 = squared_length(d);
            if (m2 > (4 + m2) * std::numeric_limits<NT>::epsilon()) {
                const NT t = Clamp_internal<NT>(-a * d / m2, 0, 1);
                const Vector3d p = a + d * t;
                const NT l2 = squared_length(p);
                if (l2 < sqd) {
                    prj = p;
                    sqd = l2;
                }
            } else {
                const NT l2 = squared_length(a);
                if (l2 < sqd) {
                    prj = a;
                    sqd = l2;
                }
            }
        }

        template<class Vector3d, typename NT>
        static inline void ProjectOriginW(const Vector3d &a,
                                         const Vector3d &b,
                                         Vector3d &prj,
                                         NT& w,
                                         NT &sqd) {
            const Vector3d d = b - a;
            auto squared_length = [](const Vector3d &d) { return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; };
            const NT m2 = squared_length(d);
            if (m2 > (squared_length(b) + squared_length(a) + m2) * std::numeric_limits<NT>::epsilon()) {
                const NT t = Clamp_internal<NT>(-a * d / m2, 0, 1);
                const Vector3d p = a + d * t;
                const NT l2 = squared_length(p);
                if (l2 < sqd) {
                    prj = p;
                    sqd = l2;
                    w = t;
                }
            } else {
                const NT l2 = squared_length(a);
                if (l2 < sqd) {
                    prj = a;
                    sqd = l2;
                    w = 0;
                }
            }
        }

        template<class Vector3d, typename NT>
        static inline bool ProjectOrigin(const Vector3d &a,
                                         const Vector3d &b,
                                         const Vector3d &c,
                                         Vector3d &prj,
                                         NT &sqd, NT* md) {
            auto squared_length = [](const Vector3d &d) { return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; };
            auto cross_product = [](const Vector3d &a, const Vector3d &b) {
                return Vector3d(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
            };
            const Vector3d q = cross_product(b - a, c - a);
            const double m2 = squared_length(q);
            if (m2 > (squared_length(b-a)*squared_length(c-a)) * std::numeric_limits<NT>::epsilon()) {
                const double aq = a * q;
                const double k2 = aq * aq / m2;
                if (md && (k2 > (*md)) ) return false;
                if (k2 < sqd) {
                    const Vector3d p = aq / m2 * q;
                    if ((cross_product(a - p, b - p) * q > 0) &&
                        (cross_product(b - p, c - p) * q > 0) &&
                        (cross_product(c - p, a - p) * q > 0)) {
                        prj = p;
                        sqd = k2;
                        return true;
                    } else {
                        ProjectOrigin(a, b, prj, sqd);
                        ProjectOrigin(b, c, prj, sqd);
                        ProjectOrigin(c, a, prj, sqd);
                        return (!md || (sqd < *md));
                    }
                }
            } else {
                ProjectOrigin(a, b, prj, sqd);
                return (!md || (sqd < *md));
            }
        }

        template<class Vector3d, typename NT>
        struct ProjectOnSegment_internals {
            static bool
            exec(const Vector3d &a, const Vector3d &b, const Vector3d &o, Vector3d &prj, NT &w, NT &sqd, NT* md) {
                sqd = DBL_MAX;
                ProjectOriginW(a - o, b - o, prj, w, sqd);
                prj = prj + o;
                return (!md || sqd < *md);
            }
        };

        template<typename NT>
        struct ProjectOnSegment_internals<typename CGAL::Simple_cartesian<NT>::Point_3, NT> {
            using Vector3d = typename CGAL::Simple_cartesian<NT>::Point_3;

            static inline bool
            exec(const Vector3d &a, const Vector3d &b, const Vector3d &o, Vector3d &prj, NT& w, NT &sqd, NT* md) {
                typename CGAL::Simple_cartesian<NT>::Segment_3 t(a, b);
                using Kernel = CGAL::Simple_cartesian<NT>;
                Kernel ker;
                typename Kernel::Construct_projected_point_3 projector = ker.construct_projected_point_3_object();
                prj = projector(t, o);
                sqd = (prj - o).squared_length();
                w = (prj - a) * (b - a) / (b - a).squared_length();
                return (!md || sqd < *md);
            }
        };

        template<class Vector3d, typename NT>
        static inline bool
        ProjectOnSegment(const Vector3d &a, const Vector3d &b, const Vector3d &o, Vector3d &prj, NT &w, NT &sqd, NT* md) {
            return ProjectOnSegment_internals<Vector3d, NT>::exec(a, b, o, prj, w, sqd, md);
        }

        template<class Vector3d, typename NT>
        struct ProjectOnTriangle_internals {
            static bool
            exec(const Vector3d &a, const Vector3d &b, const Vector3d &c, const Vector3d &o, Vector3d &prj, NT &sqd, NT* md) {
                sqd = DBL_MAX;
                auto vprj = prj - o;
                if (!ProjectOrigin(a - o, b - o, c - o, vprj, sqd, md)) return false;
                prj = o + vprj;
                return true;
            }
        };

        template<typename NT>
        struct ProjectOnTriangle_internals<typename CGAL::Simple_cartesian<NT>::Point_3, NT> {
            using Vector3d = typename CGAL::Simple_cartesian<NT>::Point_3;

            static inline bool
            exec(const Vector3d &a, const Vector3d &b, const Vector3d &c, const Vector3d &o, Vector3d &prj, NT &sqd, NT* md) {
                typename CGAL::Simple_cartesian<NT>::Triangle_3 t(a, b, c);
                if (md && (squared_distance(t.supporting_plane(), o) >= *md)) return false;
                using Kernel = CGAL::Simple_cartesian<NT>;
                Kernel ker;
                typename Kernel::Construct_projected_point_3 projector = ker.construct_projected_point_3_object();
                prj = projector(t, o);
                sqd = (prj - o).squared_length();
                return !md || sqd < *md;
            }
        };

        template<typename NT>
        struct ProjectOnTriangle_internals<typename CGAL::Simple_cartesian<NT>::Vector_3, NT> {
            using Vector3d = typename CGAL::Simple_cartesian<NT>::Vector_3;

            static inline bool
            exec(const Vector3d &a, const Vector3d &b, const Vector3d &c, const Vector3d &o, Vector3d &prj, NT &sqd, NT* md) {
                typename CGAL::Simple_cartesian<NT>::Triangle_3 t(CGAL::ORIGIN + a, CGAL::ORIGIN + b, CGAL::ORIGIN + c);
                if (md && (squared_distance(t.supporting_plane(), CGAL::ORIGIN+ o) >= *md)) return false;
                using Kernel = CGAL::Simple_cartesian<NT>;
                Kernel ker;
                typename Kernel::Construct_projected_point_3 projector = ker.construct_projected_point_3_object();
                prj = projector(t, CGAL::ORIGIN+o) - CGAL::ORIGIN;
                sqd = (prj - o).squared_length();
                return !md || sqd < *md;
            }
        };

        template<class Vector3d, typename NT>
        static inline bool ProjectOnTriangle(const Vector3d &a,
                                             const Vector3d &b,
                                             const Vector3d &c,
                                             const Vector3d &o,
                                             Vector3d &prj,
                                             NT &sqd, NT* md) {
            return ProjectOnTriangle_internals<Vector3d, NT>::exec(a, b, c, o, prj, sqd, md);
        }

        template<class Vector3d>
        static inline Vector3d BaryCoord(const Vector3d &a,
                                         const Vector3d &b,
                                         const Vector3d &c,
                                         const Vector3d &p) {
            auto cross_product = [](const auto &a, const auto &b) {
                return Vector3d(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
            };
            auto norm = [](const auto &a) { return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); };
            const auto w0 = norm(cross_product(a - p, b - p)), w1 = norm(cross_product(b - p, c - p)), w2 = norm(
                    cross_product(c - p, a - p));
            const auto isum = 1 / (w0 + w1 + w2);
            return (Vector3d(w1 * isum, w2 * isum, w0 * isum));
        }

        template<typename T, typename R>
        struct BaryEval_internals{
               static T exec(const T &a, const T &b, const T &c, const R &coord){
                   return (a * coord[0] + b * coord[1] + c * coord[2]);
               }
        };
        template<typename R>
        struct BaryEval_internals<typename CGAL::Simple_cartesian<DReal>::Point_3, R>{
            using T = CGAL::Simple_cartesian<DReal>::Point_3;
            static T exec(const T &a, const T &b, const T &c, const R &coord){
                return CGAL::ORIGIN + ((a - CGAL::ORIGIN) * coord[0] + (b - CGAL::ORIGIN) * coord[1] + (c - CGAL::ORIGIN) * coord[2]);
            }
        };

        template<typename T, typename R>
        static inline T BaryEval(const T &a,
                                 const T &b,
                                 const T &c,
                                 const R &coord) {
            return BaryEval_internals<T, R>::exec(a, b, c, coord);
        }

        template<typename T, typename R>
        struct MakeLerp{
              static T exec(const T& a, const T& b, const R t) { return a * (1 - t) + b * t; }
        };
        template<typename R>
        struct MakeLerp<Point, R>{
            using T = Point;
            static T exec(const T& a, const T& b, const R t) { return Point(a[0]*(1-t) + b[0]*t, a[1]*(1-t) + b[1]*t, a[2]*(1-t) + b[2]*t); }
        };
        template<typename R>
        struct MakeLerp<Point_2, R>{
            using T = Point_2;
            static T exec(const T& a, const T& b, const R t) { return Point_2(a[0]*(1-t) + b[0]*t, a[1]*(1-t) + b[1]*t); }
        };
        template<typename T, typename R>
        static inline T Lerp(const T& a, const T& b, const R t){
            return MakeLerp<T, R>::exec(a, b, t);
        }

        //(1-w0) * a0 + w0 * a1
        //(1-w1) * b0 + w1 * b1
        template<class Vector3d, typename NT>
        static inline bool
        SegmentsProjects(const Vector3d &a0, const Vector3d &a1, const Vector3d &b0, const Vector3d &b1,
                         Vector3d &prj0, Vector3d &prj1, NT &w0, NT &w1, NT &sqd, NT* md) {
            auto squared_length = [](const auto &d) { return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; };
            auto dot_prod = [](const auto &a, const auto &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
            auto abs_vec = [](const auto &a){ return Vector3d(abs(a[0]), abs(a[1]), abs(a[2])); };
            auto cross_product = [](const auto &a, const auto &b) { return Vector3d(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]); };
            auto max_ind = [](const auto &a){ 
                int ind = 0; 
                auto m = abs(a[0]); 
                for (int i = 1; i < 3; ++i) if (abs(a[i]) > m) { m = abs(a[i]); ind = i; }
                return ind;
            };
            auto dr = a0 - b0;
            auto l0 = a1 - a0, l1 = b1 - b0;
            auto l0_2 = squared_length(l0), l1_2 = squared_length(l1), l01 = dot_prod(l0, l1);
            auto cross_l = cross_product(l0, l1);
            auto det = squared_length(cross_l);//l0_2 * l1_2 - l01 * l01;
            if (det <= (l0_2 * l1_2 + l01 * l01) * std::numeric_limits<NT>::epsilon()) {
                auto choose_part = [&](double* r0, double* r1, const Vector3d& n0, double db0){
                    if (r1[1] <= r0[0]){
                        sqd = squared_length(a0 - b1);
                        w0 = 0, w1 = 1;
                        prj0 = a0, prj1 = b1;
                    } else if (r1[0] >= r0[1]){
                        sqd = squared_length(a1 - b0);
                        w0 = 1, w1 = 0;
                        prj0 = a1, prj1 = b0;
                    } else {
                        double r = (std::max(r0[0], r1[0]) + std::min(r0[1], r1[1])) / 2;
                        sqd = squared_length(Vector3d(-db0*n0[0], -db0*n0[1], -db0*n0[2]) - dr);
                        w0 = (r - r0[0]) / (r0[1] - r0[0]);
                        if (std::isnan(w0) || w0 > 1 || w0 < 0) w0 = 0.5;
                        w1 = (r - r1[0]) / (r1[1] - r1[0]);
                        if (std::isnan(w1) || w1 > 1 || w1 < 0) w1 = 0.5;
                        prj0 = a0 + w0 * (a1 - a0),  prj1 = b0 + w1 * (b1 - b0);
                    }
                };
                if (l0_2 > 2 * dot_prod(abs_vec(a0) + abs_vec(a1), Vector3d(1, 1, 1)) * std::numeric_limits<NT>::epsilon()){
                    double ll0 = sqrt(l0_2), ll1 = sqrt(l1_2);
                    Vector3d n0 = l0 / ll0;
                    double db0 = dot_prod(-dr, n0);
                    double r0[2] = {0, ll0}, r1[2] = {db0, db0 + ll1};
                    choose_part(r0, r1, n0, db0);
                } else if (l1_2 > 2 * dot_prod(abs_vec(b0) + abs_vec(b1), Vector3d(1, 1, 1)) * std::numeric_limits<NT>::epsilon()){
                    double ll0 = sqrt(l0_2), ll1 = sqrt(l1_2);
                    Vector3d n0 = l1 / ll1;
                    double db0 = dot_prod(-dr, n0);
                    double r0[2] = {0, ll0}, r1[2] = {db0, db0 + ll1};
                    choose_part(r0, r1, n0, db0);
                } else {
                    prj0 = (a0 + a1) / 2, prj1 = (b0 + b1) / 2;
                    w0 = w1 = 0.5;
                    sqd = squared_length(prj0 - prj1);
                }
            }
            else {
                auto drx = dr - dot_prod(dr, cross_l) * cross_l / det;
                auto cross_0 = cross_product(drx, l1), cross_1 = cross_product(drx, l0);
                auto ind = max_ind(cross_l);
                auto t0r = Clamp_internal<NT>(-cross_0[ind] / cross_l[ind], 0, 1),
                     t1r = Clamp_internal<NT>(-cross_1[ind] / cross_l[ind], 0, 1);
                bool in0 = t0r > 0 && t0r < 1, in1 = t1r > 0 && t1r < 1;
                auto sqd1 = sqd;
                if (!in0 || !in1) {
                    auto tmp = a0;
                    if (in0){
                        ProjectOnSegment(a0, a1, b0 + l1 * t1r, tmp, t0r, sqd1);
                    } else if (in1){
                        ProjectOnSegment(b0, b1, a0 + l0 * t0r, tmp, t1r, sqd1);
                    } else {
                        sqd1 = squared_length((a0 + l0 * t0r) - (b0 + l1 * t1r));
                        auto sqd2 = sqd1;
                        auto w = t0r;
                        if (ProjectOnSegment(a0, a1, b0, tmp, w, sqd2, &sqd1)){
                            sqd1 = sqd2;
                            t0r = w;
                            t1r = 0;
                        }
                        if (ProjectOnSegment(a0, a1, b1, tmp, w, sqd2, &sqd1)){
                            sqd1 = sqd2;
                            t0r = w;
                            t1r = 1;
                        }
                        if (ProjectOnSegment(b0, b1, a0, tmp, w, sqd2, &sqd1)){
                            sqd1 = sqd2;
                            t1r = w;
                            t0r = 0;
                        }
                        if (ProjectOnSegment(b0, b1, a1, tmp, w, sqd2, &sqd1)){
                            sqd1 = sqd2;
                            t1r = w;
                            t0r = 1;
                        }
                    }
                    prj0 = a0 + l0 * t0r;
                    prj1 = b0 + l1 * t1r;
                } else {
                    prj0 = a0 + l0 * t0r;
                    prj1 = b0 + l1 * t1r;
                    sqd1 = squared_length(prj0 - prj1);
                }
                sqd = sqd1;
                w0 = t0r, w1 = t1r;
            }
            return (!md || sqd < *md);
        }

        template<class Vector3d, typename NT> 
        static inline bool 
        FaceEdgeProjects(const Vector3d &f0, const Vector3d &f1, const Vector3d &f2, const Vector3d &e0, const Vector3d &e1, 
                         Vector3d &prj_f, Vector3d &prj_e, NT* w_f/*[3]*/, NT* w_e/*[2]*/, NT& sqd, NT* max_sqr_dist){
            auto squared_length = [](const auto &d) { return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; };
            auto dot_prod = [](const auto &a, const auto &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
            auto abs_vec = [](const auto &a){ return Vector3d(abs(a[0]), abs(a[1]), abs(a[2])); };
            auto cross_product = [](const auto &a, const auto &b) { return Vector3d(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]); };

            std::array<NT, 5> dists = {-1, -1, -1, -1, -1};
            std::array<NT, 3> lw_f, lw_e;
            std::array<Vector3d, 5> prjs_f, prjs_e;
            std::array<const Vector3d*, 3> fp = {&f0, &f1, &f2};

            //check whether triangle is degenerate
            auto l0 = f1 - f0, l1 = f2 - f0;
            auto l0_2 = squared_length(l0), l1_2 = squared_length(l1), l01 = dot_prod(l0, l1);
            auto cross_l = cross_product(l0, l1);
            auto det = squared_length(cross_l);//l0_2 * l1_2 - l01 * l01;
            if (det <= (l0_2 * l1_2 + l01 * l01) * std::numeric_limits<NT>::epsilon()){ //triangle is degenerate
                auto l2_2 = squared_length(f2 - f1);
                unsigned char max_ind = 0;
                auto max_l2 = l0_2;
                if (l1_2 > max_l2) max_l2 = l1_2, max_ind = 2;
                if (l2_2 > max_l2) max_l2 = l2_2, max_ind = 1;
                if (SegmentsProjects(*(fp[max_ind]), *(fp[(max_ind+1)%3]), e0, e1, prjs_f[max_ind], prjs_e[max_ind], lw_f[max_ind], lw_e[max_ind], dists[max_ind], max_sqr_dist)){
                    prj_f = prjs_f[max_ind], prj_e = prjs_e[max_ind];
                    if (w_f != nullptr){
                        w_f[max_ind] = 1 - lw_f[max_ind];
                        w_f[(max_ind+1)%3] = lw_f[max_ind];
                        w_f[(max_ind+2)%3] = 0;
                    }
                    if (w_e != nullptr){
                        w_e[0] = 1 - lw_e[max_ind];
                        w_e[1] = lw_e[max_ind];
                    }
                    sqd = dists[max_ind];
                    return true;
                }
                return false;
            }

            //check whether line segment is degenerate
            auto l3 = e1 - e0;
            auto l3_2 = squared_length(l3);
            if (l3_2 <= (l0_2 + l1_2) * std::numeric_limits<NT>::epsilon()){ //edge is degenerate
                Vector3d p = e0 + (e1 - e0)/2;
                if(ProjectOnTriangle(f0, f1, f2, p, prjs_f[3], sqd, max_sqr_dist)){
                    prj_f = prjs_f[3], prj_e = p;   
                    if (w_f != nullptr){
                        auto w = GeomProjs::BaryCoord(f0, f1, f2, prjs_f[3]);
                        for (unsigned char i = 0; i < 3; ++i) w_f[i] = w[i];
                    }  
                    if (w_e != nullptr){
                        w_e[0] = w_e[1] = 0.5;
                    }                   
                    return true;
                }
                return false;
            }

            //whether line segment intersects interior of the triangle
            auto dp = dot_prod(cross_l, l3);
            if (dp * dp >  det * l3_2 * std::numeric_limits<NT>::epsilon()){ //triangle and edge are not parallel
                auto t = -dot_prod(cross_l, e0 - f0) / dp;
                if (t >= 0 && t <= 1){
                    Vector3d p = e0 + t * l3;
                    if ((dot_prod(cross_product(f0 - p, f1 - p), cross_l) > 0) &&
                        (dot_prod(cross_product(f1 - p, f2 - p), cross_l) > 0) &&
                        (dot_prod(cross_product(f2 - p, f0 - p), cross_l) > 0)){
                            prj_f = p, prj_e = p;
                            if (w_f != nullptr){
                                auto w = GeomProjs::BaryCoord(f0, f1, f2, p);
                                for (unsigned char i = 0; i < 3; ++i) w_f[i] = w[i];
                            }  
                            if (w_e != nullptr){
                                w_e[0] = 1 - t;
                                w_e[1] = t;
                            }
                            sqd = 0;
                            return true;  
                        }
                }
            }
            //either 
            //1. nearest point on the edge is its end 
            //2. nearest point on the triangle on its boundary
            unsigned char is_near_enough_e = 0, is_near_enough_p = 0;
            for (unsigned char i = 0; i < 3; ++i)
                if (!SegmentsProjects(*(fp[i]), *(fp[(i+1)%3]), e0, e1, prjs_f[i], prjs_e[i], lw_f[i], lw_e[i], dists[i], max_sqr_dist)) 
                    dists[i] = -1;
                else 
                    ++is_near_enough_e;
            if(!ProjectOnTriangle(f0, f1, f2, e0, prjs_f[3], dists[3], max_sqr_dist)) dists[3] = -1, prjs_e[3] = e0; else ++is_near_enough_p;
            if(!ProjectOnTriangle(f0, f1, f2, e1, prjs_f[4], dists[4], max_sqr_dist)) dists[4] = -1, prjs_e[4] = e1; else ++is_near_enough_p;  
            if (is_near_enough_e == 0 && is_near_enough_p == 0) return false;

            unsigned char s = 0; //shift in remap
            std::array<unsigned char, 5> r = {0, 1, 2, 3, 4}; //remap
            std::sort(r.begin(), r.end(), [&dists](unsigned char a, unsigned char b) { return dists[a] < dists[b]; });
            for (int i = 0; i < r.size(); ++i) if (r[i] >= 0) { s = i; break;}
            
            prj_f = prjs_f[r[s]], prj_e = prjs_e[r[s]];
            if (w_f != nullptr){
                if (r[s] < 3){
                    w_f[r[s]] = 1 - lw_f[r[s]];
                    w_f[(r[s]+1)%3] = lw_f[r[s]];
                    w_f[(r[s]+2)%3] = 0;
                } else {
                    auto w = GeomProjs::BaryCoord(f0, f1, f2, prjs_f[r[s]]);
                    for (unsigned char i = 0; i < 3; ++i) w_f[i] = w[i];
                }
            }
            if (w_e != nullptr){
                if (r[s] < 3){
                    w_e[0] = 1 - lw_e[r[s]];
                    w_e[1] = lw_e[r[s]];
                } else {
                    if (r[s] == 3) w_e[0] = 1, w_e[1] = 0;
                    else w_e[0] = 0, w_e[1] = 1;
                }
            }
            sqd = dists[r[s]];

            return true;
        }
        template<class Vector3d, typename NT> 
        static inline bool 
        FaceFaceProjects(const Vector3d &f0, const Vector3d &f1, const Vector3d &f2, const Vector3d &q0, const Vector3d &q1, const Vector3d &q2,
                         Vector3d &prj_f, Vector3d &prj_q, NT* w_f/*[3]*/, NT* w_q/*[3]*/, NT& sqd, NT* max_sqr_dist){
            auto squared_length = [](const auto &d) { return d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; };
            auto dot_prod = [](const auto &a, const auto &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
            auto abs_vec = [](const auto &a){ return Vector3d(abs(a[0]), abs(a[1]), abs(a[2])); };
            auto cross_product = [](const auto &a, const auto &b) { return Vector3d(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]); };

            std::array<Vector3d, 6> vp{f0, f1, f2, q0, q1, q2};
            std::array<NT*, 2> ww{w_f, w_q};
            std::array<Vector3d, 6> lp{f1-f0, f2-f1, f0-f2, q1-q0, q2-q1, q0-q2};
            std::array<NT, 6> lp_2;
            std::array<Vector3d, 2> nf; std::array<NT, 2> nf_2;
            //check whether there is degenerate triangle
            for (unsigned char k = 0; k < 2; ++k){
                for (unsigned char i = 0; i < 3; ++i) lp_2[k*3 + i] = squared_length(lp[k*3 + i]); 
                nf[k] = cross_product(lp[3*k + 2], lp[3*k + 0]);
                nf_2[k] = squared_length(nf[k]);
                auto l20 = dot_prod(lp[3*k + 2], lp[3*k + 0]);
                if (nf_2[k] <= (lp_2[3*k + 0] * lp_2[3*k + 2] + l20 * l20) * std::numeric_limits<NT>::epsilon()){ //the k-th triangle is degenerate
                    NT* it = std::max_element(lp_2.data() + 3*k, lp_2.data() + 3*(k+1));
                    unsigned char j = std::distance(lp_2.data() + 3*k, it);
                    std::array<NT, 6> lw;
                    if (FaceEdgeProjects(vp[3*((k+1)%2)+0], vp[3*((k+1)%2)+1], vp[3*((k+1)%2)+2], vp[3*k+j], vp[3*k+(j+1)%3], prj_q, prj_f, lw.data(), lw.data() + 3, sqd, max_sqr_dist)){
                        if (k == 1) std::swap(prj_q, prj_f);
                        if (ww[k] != nullptr){
                             ww[k][j] = lw[3+0];
                             ww[k][(j+1)%3] = lw[3+1];
                             ww[k][(j+2)%3] = 0;
                        }
                        if (ww[(k+1)%2] != nullptr) std::copy(lw.data(), lw.data()+3, ww[(k+1)%2]);
                        return true;
                    }
                    return false;
                }
            }
            //check every edge pair
            bool is_guarant_separated = false;
            unsigned int ief = 0, ieq = 0;
            NT esqd = 0;
            std::array<NT, 2> elw;
            std::array<Vector3d, 2> eprjs;
            {
                auto separate_dir = [&cross_product, &dot_prod, &squared_length](const Vector3d& a, const Vector3d& b, NT w0, const Vector3d& f0,
                                                                const Vector3d& c, const Vector3d& d, NT w1, const Vector3d& f1){
                    auto get_case = [](NT w){ return (w >= 1) ? 2 : ((w <= 0) ? 0 : 1); };
                    Vector prj0 = (1 - w0) * a + w0 * b, prj1 = (1 - w1) * c + w1 * d;

                    int cs = get_case(w0) + 3*get_case(w1);
                    #define ID(i, j) (3*i + j)
                    switch (cs){
                        case ID(1,1): {
                            auto s = cross_product(b - a, d - c);
                            if (squared_length(s) > (squared_length(b-a) * squared_length(d-c))*std::numeric_limits<NT>::epsilon()){
                                if (dot_prod(s, c-a) < 0) s *= -1;
                                return s;
                            } else {//edges are parallel
                                Vector3d l = b - a;
                                NT l_2 = squared_length(l);
                                Vector3d e0 = f0 - (f0*l)*l/l_2, e1 = f1 - (f1*l)*l/l_2;
                                Vector3d s0 = cross_product(e0, l), s1 = cross_product(e1, l);
                                s0 /= sqrt(squared_length(s0)); s1 /= sqrt(squared_length(s1));
                                Vector3d s = (s0 + s1) / 2;
                                if (abs(squared_length(s)) <= 8*std::numeric_limits<NT>::epsilon()){
                                    s = s0;
                                } else s /= sqrt(squared_length(s0));
                                if (dot_prod(s, c-a) < 0) s *= -1;
                                return s;
                            }
                        }
                        case ID(0, 1): return cross_product(b - a, cross_product(prj1 - prj0, b - a));
                        case ID(1, 0): return cross_product(d - c, cross_product(prj1 - prj0, d - c));
                        case ID(1, 2): return cross_product(d - c, cross_product(prj1 - prj0, d - c));
                        case ID(2, 1): return cross_product(b - a, cross_product(prj1 - prj0, b - a));
                        case ID(0, 0): return c - a;
                        case ID(0, 2): return c - b;
                        case ID(2, 0): return d - a;
                        case ID(2, 2): return d - b;
                    }
                    #undef ID
                    return a;
                };
                std::array<Vector3d, 2> lprjs; 
                std::array<NT, 2> lw;
                NT lsqd;
                NT lmin_sqr_dist = squared_length(q0 - f0)+1;
                for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j){
                    if (SegmentsProjects(vp[i], vp[(i+1)%3], vp[3 + j], vp[3 + (j+1)%3], lprjs[0], lprjs[1], lw[0], lw[1], lsqd, &lmin_sqr_dist)){
                        ief = i, ieq = j;
                        esqd = lsqd;
                        elw = lw;
                        eprjs = lprjs;
                        if (lsqd < lmin_sqr_dist * std::numeric_limits<NT>::epsilon()){ //edge-edge intersection
                            prj_f = lprjs[0], prj_q = lprjs[1];
                            if (w_f != nullptr) w_f[i] = (1 - lw[0]), w_f[(i+1)%3] = lw[0], w_f[(i+2)%3] = 0;
                            if (w_q != nullptr) w_q[j] = (1 - lw[1]), w_q[(j+1)%3] = lw[1], w_q[(j+2)%3] = 0;
                            sqd = 0;
                            return (!max_sqr_dist || sqd < *max_sqr_dist);
                        }
                        lmin_sqr_dist = lsqd;
                        auto s = separate_dir(vp[i], vp[(i+1)%3], lw[0], vp[(i+2)%3] - vp[i], vp[3 + j], vp[3 + (j+1)%3], lw[1], vp[3+(j+2)%3] - vp[3 + j]);
                        auto s0 = dot_prod(vp[(i+2)%3] - lprjs[0], s),
                             s1 = dot_prod(vp[3+(j+2)%3] - lprjs[1], s);
                        if (s0 <= 0 && s1 >= 0){ //edge-edge nearest points are the nearest points between triangles 
                            prj_f = lprjs[0], prj_q = lprjs[1];
                            if (w_f != nullptr) w_f[i] = (1 - lw[0]), w_f[(i+1)%3] = lw[0], w_f[(i+2)%3] = 0;
                            if (w_q != nullptr) w_q[j] = (1 - lw[1]), w_q[(j+1)%3] = lw[1], w_q[(j+2)%3] = 0;
                            sqd = lsqd;
                            return (!max_sqr_dist || sqd < *max_sqr_dist);
                        }
                        auto d = dot_prod(lprjs[1] - lprjs[0], s);
                        if (s0 < 0) s0 = 0;
                        if (s1 > 0) s1 = 0;
                        if (d - s0 + s1 > 0) is_guarant_separated = true;
                    }
                }
            }
            
            for (unsigned char k = 0; k < 2; ++k){
                NT Tp[3];
                for (int i = 0; i < 3; ++i) Tp[i] = dot_prod(vp[3*k+i] - vp[3*((k+1)%2)], nf[(k+1)%2]);
                int ip = -1;
                if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0)){
                    if (Tp[0] < Tp[1]) ip = 0; else ip = 1;
                    if (Tp[2] < Tp[ip]) ip = 2;
                } else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0)){
                    if (Tp[0] > Tp[1]) ip = 0; else ip = 1;
                    if (Tp[2] > Tp[ip]) ip = 2;
                }
                if (ip >= 0){
                    is_guarant_separated = true;
                    if (   dot_prod(vp[3*k + ip] - vp[3*((k+1)%2)+0],cross_product(nf[(k+1)%2],lp[3*((k+1)%2)+0])) > 0 
                        && dot_prod(vp[3*k + ip] - vp[3*((k+1)%2)+1],cross_product(nf[(k+1)%2],lp[3*((k+1)%2)+1])) > 0
                        && dot_prod(vp[3*k + ip] - vp[3*((k+1)%2)+2],cross_product(nf[(k+1)%2],lp[3*((k+1)%2)+2])) > 0){
                        prj_f = vp[3*k + ip] + Tp[ip]*nf[(k+1)%2] / nf_2[(k+1)%2]; 
                        prj_q = vp[3*k + ip];
                        if (ww[(k+1)%2]){
                            auto wf = GeomProjs::BaryCoord(vp[3*((k+1)%2)], vp[3*((k+1)%2)+1], vp[3*((k+1)%2)+2], prj_f);
                            for (int i = 0; i < 3; ++i) ww[(k+1)%2][i] = wf[i];
                        }
                        if (ww[k]) ww[k][ip] = 1, ww[k][(ip+1)%3] = 0, ww[k][(ip+2)%3] = 0;
                        sqd = squared_length(prj_f - prj_q);
                        if (k != 0) std::swap(prj_f, prj_q);
                        return (!max_sqr_dist || sqd < *max_sqr_dist);
                    }
                }
            }

            if (!is_guarant_separated){ //triangles intersects
                sqd = 0;
                NT _sqd = 0;
                int _k = 0;
                std::array<NT, 3> lw[2];
                std::array<Vector3d, 2> prjs;
                for (int k = 0; k < 2; ++k)
                for (int i = 0; i < 3; ++i){
                    auto& l3 = lp[3*((k+1)%2)+i];
                    auto l3_2 = lp_2[3*((k+1)%2)+i];
                    auto dp = dot_prod(nf[k], l3);
                    auto check_point = [&](Vector3d p, NT t) -> bool{
                        if ((dot_prod(cross_product(vp[3*k+0] - p, vp[3*k+1] - p), nf[k]) > 0) &&
                            (dot_prod(cross_product(vp[3*k+1] - p, vp[3*k+2] - p), nf[k]) > 0) &&
                            (dot_prod(cross_product(vp[3*k+2] - p, vp[3*k+0] - p), nf[k]) > 0)){
                            prj_q = prj_f = p;
                            if (ww[k]){
                                auto wf = GeomProjs::BaryCoord(vp[3*k], vp[3*k+1], vp[3*k+2], prj_f);
                                for (int i = 0; i < 3; ++i) ww[k][i] = wf[i];
                            }
                            if (ww[(k+1)%2]) ww[(k+1)%2][i] = 1 - t, ww[(k+1)%2][(i+1)%3] = t, ww[(k+1)%2][(i+2)%3] = 0; 
                            return true;
                        }
                        else return false;
                    };
                    if (dp * dp >  nf_2[k] * l3_2 * std::numeric_limits<NT>::epsilon()){ //triangle and edge are not parallel
                        auto t = -dot_prod(nf[k], vp[3*((k+1)%2)+i] - vp[3*k]) / dp;
                        if (t >= 0 && t <= 1){
                            Vector3d p = vp[3*((k+1)%2)+i] + t * l3;
                            if (check_point(p, t)) return true;
                        }
                    } else {
                        auto ddp = dot_prod(nf[k], vp[3*((k+1)%2)+i] - vp[3*k]);
                        auto r_2 = squared_length(vp[3*((k+1)%2)+i] - vp[3*k]);
                        if (ddp*ddp < r_2*nf_2[k]*std::numeric_limits<NT>::epsilon()){//edge belongs to triangle plane
                            Vector3d p = vp[3*((k+1)%2)+i];
                            if (check_point(p, 0)) return true;
                        }
                    }
                }
            }

            //triangles separated
            prj_f = eprjs[0], prj_q = eprjs[1];
            if (w_f != nullptr) w_f[ief] = (1 - elw[0]), w_f[(ief+1)%3] = elw[0], w_f[(ief+2)%3] = 0;
            if (w_q != nullptr) w_q[ieq] = (1 - elw[1]), w_q[(ieq+1)%3] = elw[1], w_q[(ieq+2)%3] = 0;
            sqd = is_guarant_separated ? esqd : 0.0;
            return (!max_sqr_dist || sqd < *max_sqr_dist);
        }

        template<class Vector3d, typename NT> 
        Proximity::PQR<Vector3d, NT> Proximity::Query(const Primitive<Vector3d>& p0, const Primitive<Vector3d>& p1, NT max_sqr_dist){
            PQR<Vector3d, NT> res;
            if (p0.m_t < p1.m_t) {
                res = Query(p1, p0, max_sqr_dist);
                std::swap(res.query[0], res.query[1]);
                return res;
            }
            if (p0.m_t == NPRIMITIVES || p1.m_t == NPRIMITIVES)
                throw std::runtime_error("Wring type of primitive");

            res.query[0].prim = p0;
            res.query[1].prim = p1;
            NT* max_sqr_dist_p = max_sqr_dist > 0 ? &max_sqr_dist : nullptr; 
            #define ID(T0, T1) (static_cast<unsigned>(T0)*NPRIMITIVES + static_cast<unsigned>(T1))
            switch (ID(p0.m_t, p1.m_t)){
                case ID(POINT, POINT):{
                    res.query[0].prj = p0.m_v[0];
                    res.query[0].w[0] = 1;
                    res.query[1].prj = p1.m_v[0];
                    res.query[1].w[0] = 1;
                    res.sqd = (p0.m_v[0] - p1.m_v[0]).squared_length();
                    res.status = (max_sqr_dist <= 0 || res.sqd < max_sqr_dist);
                    break;
                }
                case ID(SEGMENT, POINT):{
                    auto &e = res.query[0].prim.m_v, &p = res.query[1].prim.m_v;
                    res.status = ProjectOnSegment(e[0], e[1], p[0], 
                                                    res.query[0].prj, res.query[0].w[1], res.sqd, max_sqr_dist_p);
                    res.query[0].w[0] = 1 - res.query[0].w[1];
                    res.query[1].w[0] = 1;
                    res.query[1].prj = p[0];
                    break;
                }
                case ID(TRIANGLE, POINT):{
                    auto &f = res.query[0].prim.m_v, &p = res.query[1].prim.m_v;
                    res.status = ProjectOnTriangle(f[0], f[1], f[2], p[0], res.query[0].prj, res.sqd, max_sqr_dist_p);
                    auto w = BaryCoord(f[0], f[1], f[2], res.query[0].prj);
                    for (int i = 0; i < 3; ++i) res.query[0].w[i] = w[i];
                    res.query[1].w[0] = 1;
                    res.query[1].prj = p[0];
                    break;
                }
                case ID(SEGMENT, SEGMENT):{
                    auto &e0 = res.query[0].prim.m_v, &e1 = res.query[1].prim.m_v;
                    res.status = SegmentsProjects(e0[0], e0[1], e1[0], e1[1],
                                        res.query[0].prj, res.query[1].prj, res.query[0].w[1], res.query[1].w[1], res.sqd, max_sqr_dist_p);
                    res.query[0].w[0] = 1 - res.query[0].w[1], res.query[1].w[0] = 1 - res.query[1].w[1];
                    break;
                }
                case ID(TRIANGLE, SEGMENT):{
                    auto &f = res.query[0].prim.m_v, &e = res.query[1].prim.m_v;
                    res.status = FaceEdgeProjects(f[0], f[1], f[2], e[0], e[1], 
                                        res.query[0].prj, res.query[1].prj, res.query[0].w.data(), res.query[1].w.data(), res.sqd, max_sqr_dist_p);
                    break;
                }
                case ID(TRIANGLE, TRIANGLE):{
                    auto &f0 = res.query[0].prim.m_v, &f1 = res.query[1].prim.m_v;
                    res.status = FaceFaceProjects(f0[0], f0[1], f0[2], f1[0], f1[1], f1[2], 
                                        res.query[0].prj, res.query[1].prj, res.query[0].w.data(), res.query[1].w.data(), res.sqd, max_sqr_dist_p);
                    break;
                }
            };
            #undef ID
            return res;
        }

        template<class Vector3d, typename NT> 
        Proximity::PQR<Vector3d, NT> Proximity::Query(PrimType t0, const std::array<Vector3d, MaxPrimitiveVertices> v0, 
                                        PrimType t1, const std::array<Vector3d, MaxPrimitiveVertices> v1,  NT max_sqr_dist){
            Primitive<Vector3d> p0; p0.m_t = t0, p0.m_v = v0;
            Primitive<Vector3d> p1; p1.m_t = t1, p1.m_v = v1;
            return Query(p0, p1, max_sqr_dist);
        }

        template<class Vector3d> 
        Proximity::Primitive<Vector3d>::Primitive(int t, std::array<Vector3d, MaxPrimitiveVertices> v): m_v{v} {
            if (t < 0 || t > NPRIMITIVES) throw std::runtime_error("Wrong Primitive Type");
            m_t = static_cast<PrimType>(t);
        }
        template<class Vector3d> 
        Proximity::Primitive<Vector3d>::Primitive(int t, const std::initializer_list<Vector3d>& v){
            if (t < 0 || t > NPRIMITIVES) throw std::runtime_error("Wrong Primitive Type");
            m_t = static_cast<PrimType>(t);
            unsigned long int cnt = v.size() < 3 ? v.size() : 3;
            std::copy(v.begin(), v.begin() + cnt, m_v.begin());
        }
        template<class Vector3d> 
        Proximity::Primitive<Vector3d>::Primitive(PrimType t, const std::initializer_list<Vector3d>& v): m_t{t} {
            unsigned long int cnt = v.size() < 3 ? v.size() : 3;
            std::copy(v.begin(), v.begin() + cnt, m_v.begin());
        }

        template<class Vector3d, typename NT>
        std::ostream& operator<<(std::ostream& out, const Proximity::PQR<Vector3d, NT>& r){
             auto print_vector = [](std::ostream& out, const Vector3d& v) -> std::ostream&{
                return out << v[0] << ", " << v[1] << ", " << v[2];
             };
             auto print_prim = [&print_vector](std::ostream& out, const Proximity::Primitive<Vector3d>& p, std::string prefix = "") -> std::ostream&{
                static const std::array<std::string, 4> PrimTypesName = {"POINT", "SEGMENT", "TRIANGLE", "UNKNOWN"};
                out << "Primitive{\n";
                std::string type_name = PrimTypesName.back();
                if (p.m_t < PrimTypesName.size()-1) type_name = PrimTypesName[static_cast<int>(p.m_t)];
                if (p.m_t == Proximity::NPRIMITIVES) type_name = "NPRIMITIVES";
                out << prefix << "\ttype = \"" << type_name << "\"\n";
                int nPrim = Proximity::Primitive<Vector3d>::PrimVertNum(p.m_t);
                if (nPrim < 0) nPrim = 0;
                out << prefix << "\tv[" << nPrim << "] = { ";
                for (int i = 0; i < nPrim-1; ++i)
                    { out << "( "; print_vector(out, p.m_v[i]); out << " ), "; }
                if (nPrim > 0) out << "( "; print_vector(out, p.m_v[nPrim-1]); out << " )";
                out << " }\n";
                out << prefix << "}";
                return out;
             };
             auto print_proximity_res = [&print_vector, &print_prim](std::ostream& out, const typename Proximity::PQR<Vector3d, NT>::PrimProximityRes& p, std::string prefix = "") -> std::ostream&{
                int nPrim = Proximity::Primitive<Vector3d>::PrimVertNum(p.prim.m_t);
                if (nPrim < 0) nPrim = 0;
                out << "PrimProximityRes{\n";
                out << prefix << " "; print_prim(out, p.prim, prefix + " ") << "\n";
                out << prefix << " " << "project = ("; print_vector(out, p.prj) << ")\n";
                out << prefix << " " << "w[" << nPrim << "] = { ";
                for (int i = 0; i < nPrim -1; ++i)
                    out << p.w[i] << ", ";
                if (nPrim > 0) out << p.w[nPrim-1];
                out << " }\n";
                out << prefix << "}";
                return out;
             };
             out << "PQR {\n";
             out << " status = " << (r.status ? "true" : "false") << "\n";
             out << " sqd = " << r.sqd << "\n";
             out << " query{\n";
             out << "  [0]:" ;
             print_proximity_res(out, r.query[0], "      ") << "\n";
             out << "  [1]:";
             print_proximity_res(out, r.query[1], "      ") << "\n";
             out << " }\n";
             out << "}";
             return out; 
        }
    };
}

#endif //AVSIM_COLLISIONMANAGER_INL
