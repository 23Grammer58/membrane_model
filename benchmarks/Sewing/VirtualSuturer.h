//
// Created by Liogky Alexey on 01.08.2022.
//

#ifndef AVSYM_VIRTUALSUTURER_H
#define AVSYM_VIRTUALSUTURER_H

#include <cmath>
#include <functional>
#include <array>
#include "AVSim/LeafTemplates/TemplatesCollection.h"
#include "AVSim/Core/TriangularMeshHelpers.h"
#include "AVSim/AorticSimulator/AVSimulator.h"
#include "AVSim/Core/MeshGen/StretchFigure.h"
#include "AVSim/Core/Forces/Forces.h"
#include "AVSim/Core/ForceAppliers/ForceAppliers.h"
#include "AVSim/Core/NSWorldWrapper.h"
#include "AVSim/Solvers/NonLinearSolverKinsol.h"
#include <AVSim/Core/Collision/BridsonCollisionManager.h>
#include "AVSim/AorticSimulator/VirtSuture/ContactSurface.h"
#include "AVSim/Core/Forces/SDFForce.h"
#include <AVSim/Core/Forces/ContactForce.h>
#include <regex>

namespace World3d{
    class SutureMapCFDI{
    protected:
        DReal J0 = 1, l_l = 1;
    public:
        DReal getJ0() const { return J0; }
        DReal getLl() const { return l_l; }
        virtual bool setup(DReal _J0, DReal _l_l) = 0;
        virtual DReal operator()(DReal s) const  = 0;
        virtual SutureMapCFDI* clone() const = 0;
        virtual bool operator()(std::vector<double>& val) const;
        virtual bool operator()(std::vector<double>& val, double _J0, double _l_l);
    };
    class SutureMapCFD{
    private:
        std::unique_ptr<SutureMapCFDI> m_invoker;
    public:
        SutureMapCFD() = default;
        SutureMapCFD(SutureMapCFD&&) = default;
        template<typename CFDT, typename std::enable_if<std::is_base_of<SutureMapCFDI, CFDT>::value>::type* = nullptr>
        SutureMapCFD(CFDT c): m_invoker{new CFDT(c)} {}
        SutureMapCFD(const SutureMapCFD& a){ if (a.m_invoker) m_invoker.reset(a.m_invoker->clone()); }
        SutureMapCFD& operator=(SutureMapCFD&&) = default;
        SutureMapCFD& operator=(const SutureMapCFD& a);
        template<typename CFDT>
        CFDT *target() { return static_cast<CFDT *>(m_invoker.get()); }
        operator bool(){ return m_invoker != nullptr; }

        DReal getJ0() const { return m_invoker ? m_invoker->getJ0() : -1; }
        DReal getLl() const { return m_invoker ? m_invoker->getLl() : -1; }
        bool setup(DReal _J0, DReal _l_l);
        DReal operator()(DReal s) const{ return m_invoker ? (*m_invoker)(s) : -1; }
        bool operator()(std::vector<double>& val) const { return m_invoker ? (*m_invoker)(val) : false; }
        bool operator()(std::vector<double>& val, double _J0, double _l_l){ return m_invoker ? (*m_invoker)(val, _J0, _l_l) : false; }
    };
    class SutureMapCenterConstant: public SutureMapCFDI{
    private:
        DReal d = 0.0;
    public:
        DReal m_a = 3;
        SutureMapCenterConstant() = default;
        SutureMapCenterConstant(DReal alpha): m_a{alpha} {}
        bool setup(DReal _J0, DReal _l_l);
        DReal operator()(DReal x) const;
        SutureMapCFDI* clone() const { return new SutureMapCenterConstant(*this); }
    };
    class SutureMapCenterExp: public SutureMapCFDI{
    private:
        DReal c, w, sqrt_pi = -1, erf_1d2w;
    public:
        DReal m_a = 3;
        SutureMapCenterExp() = default;
        SutureMapCenterExp(DReal alpha): m_a{alpha} {}
        bool setup(DReal _J0, DReal _l_l);
        DReal operator()(DReal x) const;
        SutureMapCFDI* clone() const { return new SutureMapCenterExp(*this); }
    };

    class VirtualSuturerCil{
        public:
        using Curve = std::function<World3d::Point (double)>;
        using Field = std::function<World3d::Vector (double)>;
        ///Storage for information associated with leaflet
        struct LeafInfo{
            enum ContactSurfType{
                PLANE,
                BEZIER_CILINDRIC
            };
            ObjectID m_id = ObjectID();
            DReal m_l = NAN; ///< lengths of leaflet
            DReal m_P = NAN; ///< pressure
            DReal m_Ht = NAN; ///< thickness of leaflets
            DReal m_E = NAN; ///< Young modulus
            DReal m_Ps = 450; ///< dimensionless pressure P*l / (E * Ht)
            Curve m_suture_line; ///< aortic suture line, [0, 1] -> R^3
            Field m_clamp_suture_b; ///< aortic suture clamped vector, [0, 1] -> S^2, m_clamp_suture_b(t) \perp m_suture_line(t)
            bool m_use_cilidric_clamp = false;
            std::string suture_b_tag_name = "";
            ContactSurfType m_cst = BEZIER_CILINDRIC;
            SutureMapCFD m_smap;
            

            LeafInfo& setLength(DReal l){ m_l = l; return *this; }
            LeafInfo& setPressure(DReal P){ m_P = P; return *this; }
            LeafInfo& setThickness(DReal Ht){ m_Ht = Ht; return *this; }
            LeafInfo& setDlessPressure(DReal Ps){ m_Ps = Ps; return *this; }
            LeafInfo& setYoungModulus(DReal E){ m_E = E; return *this; }
            LeafInfo& setSutureLine(Curve l){ m_suture_line = std::move(l); return *this; }
            LeafInfo& setClampVectorField(Field b){ m_clamp_suture_b = std::move(b); return *this; }
            LeafInfo& setCilindricClampField(){ m_use_cilidric_clamp = true; return *this; }
            LeafInfo& setContactSurfType(ContactSurfType ct){ m_cst = ct; return *this; }
            LeafInfo& setSutureMapDistrib(SutureMapCFD cfd){ m_smap = std::move(cfd); return *this; }
            LeafInfo& setSutureBFieldTagName(std::string name = "v:bdata"){ suture_b_tag_name = name; return *this; }
        };
        struct ViewWindowContext{
            bool with_view = false;
            int argc;
            char** argv;
        };
        struct SolverCtx{
            DReal lin_droptol = 8e-3, lin_reusetol = 8e-4;
            DReal nonlin_ftol = 1e-7, nonlin_scaled_steptol = 1e-5;
            int maxits = 50;
            int verbose_level = 2;
        };
        struct SolverData{
            NSWorldWrapper m_nsww;
            NonLinearSolverKinsol m_nlsp;
            LinearSolver m_ls;
            DReal R, abs_err_init;
            World3d::Timer* m_time;

            void initialize(World& w, DReal problem_length, World3d::Timer& time);
            void setContext(const SolverCtx& ctx);
        };
        struct InitialShiftContext{
            DReal c0 = 1, c1 = 0.1;
            DReal cKcil = 0;
        };
        struct PresForceCtx{
            DReal cP = 1.0;
        };
        struct ElasticForceCtx{
            DReal cE = 1.0, cE_Ht = 1.0;
        };
        struct BendForceCtx{
            DReal cE = 1.0, cE_Ht = 1.0; 
            bool use_clamped_bc = true;
        };
        struct FreeEdgeForceCtx{
            DReal cFE_PS = 1e3, cFE_Ht = 0.2, FE_Ht = -1;
        };
        struct DampForceCtx{
            DReal csigma = 0.0, sigma = -1;
            DampForceCtx(DReal csigma = 0.0): csigma{csigma}{}
        };
        struct SDFForceCtx{
            DReal cP = 2, cHt = 0.2, Ht = -1;
            DReal cshift = 0;
        };
        struct AorticForceCtx: public SDFForceCtx{
            AorticForceCtx(): SDFForceCtx(){ cP = 2, cHt = 0.2, Ht = -1; }
        };
        struct SelfCollideCtx{
            DReal cHt = 0.2, Ht = -1;
            DReal cPspr = 250;
            DReal Dt = 0.1;
        };
        using InitialShift = std::array<InitialShiftContext, 3>;
        using PresForce = std::array<PresForceCtx, 3>;
        using ElasticForce = std::array<ElasticForceCtx, 3>;
        using BendForce = std::array<BendForceCtx, 3>;
        using FreeEdgeForce = std::array<FreeEdgeForceCtx, 3>;
        using ContactForce = std::array<SDFForceCtx, 3>;
        using AorticForce = std::array<AorticForceCtx, 3>;
        struct FreeMembraneCxt{
            PresForce m_pres;
            ElasticForce m_elast;
            FreeEdgeForce m_free_edge;
            std::array<DReal,3> suture_dt = {1.0, 1.0, 1.0};
            std::vector<DampForceCtx> m_damps;
            SolverCtx m_solve_ctx;
            FreeMembraneCxt();
        };
        struct FreeShellCxt{
            PresForce m_pres;
            ElasticForce m_elast;
            BendForce m_bend;
            FreeEdgeForce m_free_edge;
            DampForceCtx m_damps;
            SolverCtx m_solve_ctx;
        };
        struct ContactSutureCxt{
            PresForce m_pres;
            ElasticForce m_elast;
            BendForce m_bend;
            FreeEdgeForce m_free_edge;
            DampForceCtx m_damps;
            std::vector<ContactForce> m_contact;
            std::array<DReal,3> contact_dt = {1.0, 1.0, 1.0};
            SolverCtx m_solve_ctx;
            
            ContactSutureCxt();
        };
        struct AorticSutureCxt{
            PresForce m_pres;
            ElasticForce m_elast;
            BendForce m_bend;
            FreeEdgeForce m_free_edge;
            DampForceCtx m_damps;
            ContactForce m_contact;
            AorticForce m_aortic;
            SolverCtx m_solve_ctx;
        };
        struct CustomCxt{
            std::pair<PresForce, bool> m_pres;
            std::pair<ElasticForce, bool> m_elast;
            std::pair<BendForce, bool> m_bend;
            std::pair<FreeEdgeForce, bool> m_free_edge;
            std::pair<DampForceCtx, bool> m_damps;
            std::pair<ContactForce, bool> m_contact;
            std::pair<AorticForce, bool> m_aortic;
            SolverCtx m_solve_ctx;
        };
        struct ObjForces{
            Force_ID elast = -1, bend = -1, pres = -1, free_edge = -1;
            Force_ID contact_sdf = -1, aorta_sdf = -1;
        };
        struct AppliedForces{
            std::array<ObjForces, 3> sf;
            WorldForce_ID damp = -1;
        };
        struct SutureAlgoParams{
            AppliedForces m_forces;
            SelfCollideCtx m_collide;
            InitialShift m_shifts; 
            FreeMembraneCxt m_fmc;
            std::vector<FreeShellCxt> m_fsc;
            ContactSutureCxt m_csc;
            AorticSutureCxt m_asc;
            std::vector<CustomCxt> m_custom;

            SutureAlgoParams();
        };
        
        struct CommissureGeom{
            std::vector<Point> m_A; ///<commissure points
            Point m_O; ///< circumcenter
            Point m_M; ///< midcenter
            Vector m_n; ///< blood flow direction
            DReal m_R; ///< circumradius
            DReal m_area; ///< area of commissure triangle
            std::array<DReal, 3> m_a; ///< sides of commissure triangle

            void setup();
        };
        struct CilDist: public SignedDistanceField{
            Vector C;
            Vector n;
            DReal R;
            SDF operator()(const Vector& x) const override;
            CilDist(Vector _C, Vector _n, DReal _R): C{_C}, n{_n}, R{_R} {}
        };
        struct LeafBoundaries{
            Object3D* obj;
            std::vector<V_ind> v_suture_bnd;
            std::vector<V_ind> v_free_bnd;
            std::vector<E_ind> e_suture_bnd;
            std::vector<E_ind> e_free_bnd;
            std::map<E_ind, BendingForce::BoundaryData::Entry> clamp_data;
            DReal sutureT = 0.0;
            std::vector<Vector> suture_shift;
            bool addShift(DReal dt);
            DReal init_length(){ return sqrt((obj->m_x0[v_suture_bnd.back()] - obj->m_x0[v_suture_bnd.front()]).squared_length()); }
        };
        struct SaveCxt{
            enum StepSign{
                NUMERIC,
                NAMEABLE
            };
            std::string save_directory;
            std::string prefix;
            StepSign sign_stat = NUMERIC;
            void save_step(Object3D& obj, int ileaf, int step, std::string step_name) const;
        };
        
        std::vector<Object3D> m_objs;
        //World m_w; ///< container for leaflets
        CommissureGeom m_commissure;
        //std::vector<Point> m_commissure_points; ///< commissure points (m_commissure_points[i] = 0.5 * [m_linfo[(i+1)%3].m_suture_line(1) + m_linfo[(i+2)%3].m_suture_line(0)])
        std::array<LeafInfo, 3> m_linfo;
        bool m_use_cilindric_aorta = false;
        std::shared_ptr<SignedDistanceField> m_aorta; ///< signed distance field generated by aortic surface
        ViewWindowContext m_vw_ctx;
        int m_nleafs = 0;
        SutureAlgoParams m_p;
        bool m_regenerate_elastic_force = false;
        bool m_regenerate_bending_force = false;
        std::string m_gen_force_dir = "../../../generated";
        World3d::Timer m_sym_time;
        SolverData m_solver;
        std::thread* front_end;
        SaveCxt m_saver;

        VirtualSuturerCil();
        VirtualSuturerCil& setCommissurePoints(std::array<Point, 3> p);
        VirtualSuturerCil& setLeafLengths(std::array<DReal, 3> l);
        VirtualSuturerCil& setViewCtx(bool with_view = false, int argc = 0, char** argv = nullptr);
        VirtualSuturerCil& setAorticSDF(std::shared_ptr<SignedDistanceField> aorta);
        VirtualSuturerCil& setCilindricAortic();
        VirtualSuturerCil& setSaveDirectory(std::string dir, std::string prefix = "", SaveCxt::StepSign stat = SaveCxt::NUMERIC);
        LeafInfo& pushLeaf(const Mesh& m, Curve suture_line = nullptr, Field clamp_suture_b = nullptr);
        void setup();
        LeafBoundaries geometrical_shifting_leaf(int ileaf);
        bool geometrical_shifting(int save_step = 0, std::array<LeafBoundaries, 3>* bnds = nullptr);
        // bool suture();
        bool suture_leaf(int ileaf);//< you should perform setup before call!
        void set_bdata_tag(int ileaf, LeafBoundaries& bnd);
        void save_leafs(int step, std::string step_name);
        // void start_renderer();
        // void stop_renderer();
        // void apply_aortic_constraint(std::array<LeafBoundaries, 3>& bnds);
        // void apply_aortic_constraint_forces(const AorticSutureCxt& cxt, std::array<LeafBoundaries, 3>& bnds);
        // void apply_contact_suture(std::array<LeafBoundaries, 3>& bnds);
        // void apply_contact_suture_forces(ContactSutureCxt& cxt, std::array<LeafBoundaries, 3>& bnds);
        // void apply_free_shell_suture(std::array<LeafBoundaries, 3>& bnds);
        void initialize_collider(World& w);
        bool set_collission_object_thickness(World& w, ObjectID w_id, int ileaf);
        void apply_free_shell_forces(FreeShellCxt& cxt, std::array<LeafBoundaries, 3>& bnd);
        void apply_free_membrane_suture(std::array<LeafBoundaries, 3>& bnd);
        bool solve_problem(SolverCtx ctx);
        void set_free_membrane_forces();
        void set_sdf_lambda(int ileaf, DReal lambda);
        void set_sdf_force(int ileaf, const SDFForceCtx& ctx, DReal lambda = 1.0);
        void set_elastic_force(int ileaf, const ElasticForceCtx& ctx);
        void set_damp_force(World& w, const DampForceCtx& ctx);
        void set_bending_force(int ileaf, const BendForceCtx& ctx);
        void set_clamped_bc(int ileaf, const std::map<E_ind, BendingForce::BoundaryData::Entry>& bdata);
        void set_pressure_force(int ileaf, const PresForceCtx& ctx);
        void set_free_edge_force(int ileaf, const FreeEdgeForceCtx& ctx);
        void set_aortic_force(int ileaf, const AorticForceCtx& ctx);
        void set_free_membrane_forces(int ileaf);
        void set_free_membrane_forces(Object3D& obj, const FreeMembraneCxt& ctx, int ileaf);
        void set_free_shell_forces(Object3D& obj, const FreeShellCxt& cxt, int ileaf, LeafBoundaries& bnd);
        void set_contact_suture_forces(Object3D& obj, const ContactSutureCxt& cxt, int ileaf, LeafBoundaries& bnd, int ncontact);
        void set_aortic_constraint_forces(Object3D& obj, const AorticSutureCxt& cxt, int ileaf, LeafBoundaries& bnd);
        void set_custom_forces(Object3D& obj, const CustomCxt& ctx, int ileaf, LeafBoundaries& bnd);
        static std::function<Vector(Point, Vector)> getCilinderSuturer(Point C, Vector n);
        
        static void apply_free_membrane_suture(World& w, LeafBoundaries* bnds, int nbnds, const FreeMembraneCxt& ctx,
                                               SolverData& s, WorldForce_ID& damp, double default_damp_sigma);
        static void apply_free_shell_suture(World& w, const FreeShellCxt& ctx,
                                               SolverData& s, WorldForce_ID& damp, double default_damp_sigma); 
        static void apply_contact_suture(World& w, Object3D* objs[], LeafInfo::ContactSurfType tp[], Force_ID contact_sdf[], DReal contact_dt[], int nobjs, const ContactSutureCxt& ctx,
                                               SolverData& s, WorldForce_ID& damp, double default_damp_sigma, bool zero_solve = true); 
        static void apply_aortic_constraint(World& w, const AorticSutureCxt& ctx,
                                               SolverData& s, WorldForce_ID& damp, double default_damp_sigma);    
        static void apply_custom_suture(World& w, const CustomCxt& ctx,
                                               SolverData& s, WorldForce_ID& damp, double default_damp_sigma);                                                                                                                                    
        static void start_renderer(World& w, const ViewWindowContext& ctx, std::thread* &front_end);
        static void stop_renderer(World& w, const ViewWindowContext& ctx, std::thread* &front_end);
        static void solve_collisions(SolverData& s);
        static bool solve_problem(SolverData& s, const SolverCtx& ctx, WorldForce_ID damp_id);
        static void set_elastic_force(Object3D& leaf, const LeafInfo& linfo, const ElasticForceCtx& ctx, Force_ID& fid, 
                                        const std::string& gendir, bool& regen_elast_force);
        static void set_bending_force(Object3D& leaf, const LeafInfo& linfo, const BendForceCtx& ctx, Force_ID& fid, 
                                        const std::string& gendir, bool& regen_elast_force, bool& regen_bend_force); 
        static void set_clamped_bc(Object3D& obj, Force_ID& bend_id, const std::map<E_ind, BendingForce::BoundaryData::Entry>& bdata);                                                              
        static void set_pressure_force(Object3D& leaf, const LeafInfo& linfo, const PresForceCtx& ctx, Force_ID& fid);   
        static void set_free_edge_force(Object3D& leaf, const LeafInfo& linfo, const FreeEdgeForceCtx& ctx, Force_ID& fid,
                                        Vector commissure_normal, Point commissure_origin); 
        static void set_damp_force(World& w, const DampForceCtx& ctx, WorldForce_ID& id, double default_sigma);                                                             
        static void set_sdf_force(Object3D& leaf, const std::array<LeafInfo, 3>& linfo, int ileaf, const std::array<Point, 3>& commissure_points, 
                                    const SDFForceCtx& ctx, Force_ID& fid, DReal lambda);                  
        static void set_sdf_lambda(Object3D& obj, LeafInfo::ContactSurfType tp, Force_ID& id, DReal lambda);
        static void set_aortic_force(Object3D& leaf, const LeafInfo& linfo, const AorticForceCtx& ctx, Force_ID& fid, const std::shared_ptr<SignedDistanceField>& aorta); 
        private:
        void initial_shift_and_rotate_leaf(int ileaf);
        void initial_shift_and_rotate();
        void cilindrize(int ileaf, Vector n, Vector r, Point new_origin);
        bool is_obj_intersect_suture_line(int ileaf);

        static void align_to_vectors(Object3D* leaf, Vector right, Vector up);
        std::array<LeafBoundaries, 3> map_cusps_to_suture_lines();
        static LeafBoundaries map_cusp_to_aortic_suture(Object3D& leaf, Curve& line, SutureMapCFD& cfd, Field& b);
        void set_dirichlet_bc_leaf(int ileaf);
        void set_dirichlet_bc();
        
        static DReal get_free_boundary(Object3D& leaf, std::vector<V_ind>& bnd, std::vector<E_ind>& ebnd);
        static DReal get_suture_boundary(Object3D& leaf, std::vector<V_ind>& bnd, std::vector<E_ind>& ebnd);
    };
};

#endif //AVSYM_VIRTUALSUTURER_H