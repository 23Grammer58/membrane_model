//
// Created by Liogky Alexey on 19.05.2022.
//

#ifndef AVSYM_COLLISIONMANAGERINTERFACE_H
#define AVSYM_COLLISIONMANAGERINTERFACE_H
#include "Object3D.h"

namespace World3d{
    class BoxVol{
        template <typename InputIterator, bool IsPoint>
        struct GetBoxFromPP { };

        template <typename InputIterator>
        struct GetBoxFromPP<InputIterator, true> {
            static inline void execute(InputIterator beg, InputIterator end, BoxVol& bv);
        };

        Vector  st = Vector(NAN, NAN, NAN),
                ed = Vector(NAN, NAN, NAN);
    public:
        inline Vector Center() const { return ((st + ed) / 2); }
        inline Vector Lengths() const { return (ed - st); }
        inline Vector Extents() const { return ((ed - st) / 2); }
        inline const Vector& Mins() const { return (st); }
        inline const Vector& Maxs() const { return (ed); }
        static inline BoxVol FromCE(const Vector& center, const Vector& extent);
        static inline BoxVol FromCR(const Vector& c, DReal r) { return FromCE(c, Vector{r, r, r}); }
        static inline BoxVol FromMM(const Vector& st, const Vector& ed){ BoxVol bv; bv.st = st, bv.ed = ed; return bv; }
        template<class PointIterator>
        static inline BoxVol FromPoints(PointIterator first, PointIterator last);
        static inline BoxVol FromPoints(const std::initializer_list<Point>& points){ return FromPoints(points.begin(), points.end()); }
        inline BoxVol& Expand(DReal e) { return Expand(Vector(e, e, e)); }
        inline BoxVol& Expand(const Vector& e){ st -= e, ed += e; return *this; }
        inline BoxVol& SignedExpand(const Vector& e);
        inline BoxVol& ScaleExpand(const DReal& e) { return ScaleExpand(Vector(e, e, e)); }
        inline BoxVol& ScaleExpand(const Vector& e);
        inline bool Contain(const BoxVol& a) const;
        inline int Classify(const Vector& n, DReal o, int s = 0) const;
        inline static BoxVol Merge(const BoxVol& a, const BoxVol& b);
    };


    class VolumeTreeBase{
    public:
        struct LeafIndex {
            long leaf_id = -1;
            void *ptr = nullptr;
        };
        using CollideProcessor = std::function<void(const LeafIndex l1, void* data1, const LeafIndex l2, void* data2)>;

        virtual void clear() {}
        virtual bool empty() const = 0;
        virtual LeafIndex insert(const BoxVol& box, void* data) = 0;
        virtual void update(LeafIndex id, BoxVol& new_volume) = 0;
        virtual void remove(LeafIndex id) = 0;
        virtual void optimize() {};
        [[nodiscard]] virtual std::shared_ptr<VolumeTreeBase> copy_empty() const = 0;

        int drive_id = -1; ///< mark of specific implementation
        friend void VolumeIntersect(VolumeTreeBase& a, VolumeTreeBase& b, VolumeTreeBase::CollideProcessor cp);
    private:
        virtual void collide(VolumeTreeBase& other, const VolumeTreeBase::CollideProcessor& cp) = 0;
    };
    inline void VolumeIntersect(VolumeTreeBase& a, VolumeTreeBase& b, VolumeTreeBase::CollideProcessor cp){
        assert(a.drive_id == b.drive_id && "Intersection is defined only for same implementations");
        a.collide(b, cp);
    }
    using VTFactry = std::function<std::shared_ptr<VolumeTreeBase>()>;

    namespace GeomProjs {
        ///@param [in] a,b - ends of the segment
        ///@param [in] o - reference point
        ///@param [out] prj - projection of point on the object
        ///@param [out] w - such value that prj = (1-w)*a + w*b
        ///@param [out] sqd - square distance from point to projection
        ///@param [in] max_sqr_dist - square distance limit for search
        ///@return true if square distance from object to point less than max_sqr_dist and fill prj and sqd values
        ///otherwise return false and fill prj and sqd values by waste
        template<class Vector3d, typename NT>
        static inline bool ProjectOnSegment(const Vector3d &a, const Vector3d &b, const Vector3d &o, Vector3d &prj, NT &w, NT &sqd, NT* max_sqr_dist = nullptr);

        ///@param [in] a,b,c - vertices of the triangle
        ///@param [in] o - reference point
        ///@param [out] prj - projection of point on the object
        ///@param [out] sqd - square distance from point to projection
        ///@param [in] max_sqr_dist - square distance limit for search
        ///@return true if square distance from object to point less than max_sqr_dist and fill prj and sqd values
        ///otherwise return false and fill prj and sqd values by waste
        template<class Vector3d, typename NT>
        static inline bool ProjectOnTriangle(const Vector3d &a, const Vector3d &b, const Vector3d &c,
                                             const Vector3d &o, Vector3d &prj, NT &sqd, NT* max_sqr_dist = nullptr);

        template<class Vector3d>
        static inline Vector3d BaryCoord(const Vector3d &a, const Vector3d &b, const Vector3d &c, const Vector3d &p);

        template<typename T, typename R>
        static inline T BaryEval(const T &a, const T &b, const T &c, const R &coord);

        template<typename T, typename R>
        static inline T Lerp(const T& a, const T& b, const R t);

        ///@param [in] a0,a1 - ends of first segment
        ///@param [in] b0,b1 - ends of second segment
        ///@param [out] w0,w1 - such numbers that prj0 = (1 - w0)*a0 + w0*a1, prj1 = (1 - w1)*b0 + w1*b1
        ///@param [out] prj0,prj1 - the nearest points between segments on first and second segment respectively
        ///@param [out] sqd - square distance from segments
        ///@param [in] max_sqr_dist - square distance limit for search
        ///@return true if square distance between segments less than max_sqr_dist and fill prj and sqd values
        ///otherwise return false and fill prj and sqd values by waste
        template<class Vector3d, typename NT>
        static inline bool SegmentsProjects(const Vector3d &a0, const Vector3d &a1, const Vector3d &b0, const Vector3d &b1,
                                            Vector3d &prj0, Vector3d &prj1, NT &w0, NT &w1, NT &sqd, NT* max_sqr_dist = nullptr);
        ///@param [in] f0,f1,f2 - vertices of the triangle face
        ///@param [in] e0,e1 - ends of the line segment
        ///@param [out] w_f[3],w_e[2] - baricentric coordinates of nearest points on the triangle and the line segment
        ///@param [out] prj_f,prj_e - the nearest points on the triangle and the line segment \f(prj_f = \sum_{i=0}^3 w_f[i] \cdot f\{i\}\f), \f(prj_e = \sum_{i=0}^2 w_e[i] \cdot e\{i\}\f)
        ///@param [out] sqd - square distance from primitives
        ///@param [in] max_sqr_dist - square distance limit for search
        ///@return true if square distance between segments less than max_sqr_dist and fill prj and sqd values
        ///otherwise return false and fill prj and sqd values by waste
        template<class Vector3d, typename NT> 
        static inline bool FaceEdgeProjects(const Vector3d &f0, const Vector3d &f1, const Vector3d &f2, const Vector3d &e0, const Vector3d &e1, 
                                            Vector3d &prj_f, Vector3d &prj_e, NT* w_f/*[3]*/, NT* w_e/*[2]*/, NT& sqd, NT* max_sqr_dist = nullptr);
        ///@param [in] f0,f1,f2 - vertices of the first triangle face
        ///@param [in] q0,q1,q2 - vertices of the second triangle face
        ///@param [out] w_f[3],w_q[3] - baricentric coordinates of nearest points on the triangle and the line segment
        ///@param [out] prj_f,prj_q - the nearest points on the triangle and the line segment \f(prj_f = \sum_{i=0}^3 w_f[i] \cdot f\{i\}\f), \f(prj_q = \sum_{i=0}^3 w_q[i] \cdot q\{i\}\f)
        ///@param [out] sqd - square distance from primitives
        ///@param [in] max_sqr_dist - square distance limit for search
        ///@return true if square distance between segments less than max_sqr_dist and fill prj and sqd values
        ///otherwise return false and fill prj and sqd values by waste
        template<class Vector3d, typename NT> 
        static inline bool FaceFaceProjects(const Vector3d &f0, const Vector3d &f1, const Vector3d &f2, const Vector3d &q0, const Vector3d &q1, const Vector3d &q2,
                                            Vector3d &prj_f, Vector3d &prj_q, NT* w_f/*[3]*/, NT* w_q/*[3]*/, NT& sqd, NT* max_sqr_dist = nullptr);
        ///Find continous collision between triangle T and point P using uniform linear interpolation between start and end positions
        ///@param a0,b0,c0,a1,b1,c1 - start and end positions of T vertices
        ///@param p0,p1 - start and end positions of point P
        ///@param[out] t - [optional] compute and return time in [0, 1] of collision if occur (v(t) = v0*(1-t) + v1*t)
        ///@return true if collision occur
        bool PerformVFCCD(const Vector& a0, const Vector& b0, const Vector& c0, const Vector& p0,
                          const Vector& a1, const Vector& b1, const Vector& c1, const Vector& p1, DReal* t = nullptr);
        ///Find continous collision between edge AB and edge CD using uniform linear interpolation between start and end positions
        ///@param a0,b0,a1,b1 - start and end positions of AB vertices
        ///@param c0,d0,c1,d1 - start and end positions of CD vertices
        ///@param[out] t - [optional] compute and return time in [0, 1] of collision if occur (v(t) = v0*(1-t) + v1*t)
        ///@return true if collision occur
        bool PerformEECCD(const Vector& a0, const Vector& b0, const Vector& c0, const Vector& d0,
                          const Vector& a1, const Vector& b1, const Vector& c1, const Vector& d1, DReal* t = nullptr);
        struct Proximity{
            enum PrimType{
                POINT = 0,
                SEGMENT,
                TRIANGLE,
                NPRIMITIVES
            };
            static const int MaxPrimitiveVertices = 3;
            ///To store description of primitives
            template<class Vector3d> 
            struct Primitive{ 
                static int PrimVertNum(PrimType t){
                    static const int PrimitiveVertices[NPRIMITIVES+1] = {1, 2, 3, -1};
                    return PrimitiveVertices[static_cast<int>(t)];
                }

                PrimType m_t = NPRIMITIVES;
                std::array<Vector3d, MaxPrimitiveVertices> m_v;
                Primitive() = default;
                Primitive(PrimType t, std::array<Vector3d, MaxPrimitiveVertices> v): m_t{t}, m_v{v} {}
                Primitive(int t, std::array<Vector3d, MaxPrimitiveVertices> v);
                Primitive(int t, const std::initializer_list<Vector3d>& v);
                Primitive(PrimType t, const std::initializer_list<Vector3d>& v);
            };
            ///To store Proximity Query Result
            template<class Vector3d, typename NT> 
            struct PQR{
                struct PrimProximityRes{
                    Primitive<Vector3d> prim;
                    Vector3d prj; ///< nearest point on the primitive
                    std::array<NT, MaxPrimitiveVertices> w; ///< weights of primitive vertices
                };
                std::array<PrimProximityRes, 2> query;
                NT sqd; ///< square distance
                bool status = false; ///< true if distance found successfully
            };
            ///Function to find nearest points between primitives
            ///@param p0,p1 primitives between nearest points will be founded
            ///@param max_sqr_dist {optional} if positive then will set PQR::status to false 
            ///if squared distance between primitives more than max_sqr_dist and then other parts of PQR will not be evaluated 
            template<class Vector3d, typename NT = DReal> 
            static PQR<Vector3d, NT> Query(const Primitive<Vector3d>& p0, const Primitive<Vector3d>& p1, NT max_sqr_dist = -1);

            template<class Vector3d, typename NT = DReal> 
            static PQR<Vector3d, NT> Query(PrimType t0, const std::array<Vector3d, MaxPrimitiveVertices> v0, 
                                           PrimType t1, const std::array<Vector3d, MaxPrimitiveVertices> v1,  NT max_sqr_dist = -1);
        };

        template<class Vector3d, typename NT>
        std::ostream& operator<<(std::ostream& out, const Proximity::PQR<Vector3d, NT>& r);
    };

    class CollisionManagerBase{
    public:
        virtual void addObj(Object3D& obj, ObjectID id, bool is_dynamic = true, double margin = 0.1) {
            (void) obj, (void) id, (void) is_dynamic, (void) margin;
        }
        virtual void removeObj(ObjectID id) { (void) id; }
        virtual void findCollisions() {}
        virtual void solveCollisions() {}
    };
};

#include "CollisionManagerInterface.inl"

#endif //AVSYM_COLLISIONMANAGERINTERFACE_H
