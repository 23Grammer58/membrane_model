#include "VirtualSuturer.h"

bool SutureMapCFDI::operator()(std::vector<double>& val) const {
    for (auto& x: val){
        x =  (*this)(x);
        if (x < 0) return false;
    }
    return true;
}

bool SutureMapCFDI::operator()(std::vector<double>& val, double _J0, double _l_l){
    if (!setup(_J0, _l_l)) return false;
    return (*this)(val);
}

SutureMapCFD& SutureMapCFD::operator=(const SutureMapCFD& a){
    if (&a == this) return *this; 
    if (a.m_invoker) m_invoker.reset(a.m_invoker->clone()); 
    return *this;
}

bool SutureMapCFD::setup(DReal _J0, DReal _l_l) {
    if (!m_invoker) return false;
    return m_invoker->setup(_J0, _l_l);
}

bool SutureMapCenterConstant::setup(DReal _J0, DReal _l_l){ 
    J0 = _J0, l_l = _l_l; 
    if (J0*m_a < 1) m_a = 1/J0;
    d = std::max(m_a * (1 - J0) / (2 * (m_a - 1)), DReal(0.0)); 
    return true; 
}
DReal SutureMapCenterConstant::operator()(DReal x) const{
    DReal y = 0;
    if (x < 0.5 - d) y = x/J0;
    else if (x < 0.5 + d) {
        DReal x0 = 0.5 - d;
        y = x0 / J0 + (x - x0) / (m_a * J0);
    }
    else 
        y = (0.5 - d)/J0 + 2*d / (m_a * J0) + (x - (0.5 + d)) / J0;
    return y;
}

bool SutureMapCenterExp::setup(DReal _J0, DReal _l_l){ 
    //p(x) = 1/J0 - c * e^((x-0.5)^2 / w^2)
    //p(0.5) = 1/(alpha*J0)
    //int p(x) dx from 0 to 1 = 1
    J0 = _J0, l_l = _l_l; 
    c = (m_a - 1) / (m_a * J0);
    w = 2.5;
    if (sqrt_pi < 0) sqrt_pi = sqrt(M_PI);
    {
        double q = 2*(1 - J0)*m_a / ((m_a - 1) * sqrt_pi);
        double m = 0, M = 10, tol = 1e-7;
        while (M - m > tol){
            double x = (m+M)/2;
            double f = x * erf(1/x);
            if (f > q) M = x;
            else m = x;
        }
        w = (M + m) / 4;
    }
    erf_1d2w = erf(1/(2*w));
    return true; 
}
DReal SutureMapCenterExp::operator()(DReal x) const{
    return x / J0 - c * sqrt_pi / 2 * w * (erf((2*x - 1) / (2*w)) + erf_1d2w);
}

void VirtualSuturerCil::SolverData::initialize(World& w, DReal problem_length, World3d::Timer& time){
    m_time = &time;
    R = problem_length;
    abs_err_init = w.getWorldNorm();
    m_nsww.resetWorld(w);
    NLProblem nlp = NSWorldWrapper::makeProblem(m_nsww);
    m_nlsp.SetNLProblem(nlp);
    m_ls = LinearSolver("inner_mptiluc");
    m_ls.setVerbosityLevel(1);
    m_ls.SetParameterReal("relative_tolerance", 1e-20);
    m_nlsp.SetLinearSolver(&m_ls);
    static_cast<SUNLinearSolverContent_Com>(static_cast<SUNLinSolCustom>(m_nlsp.GetLinearSolver()->content)->content)->verbosity = 1;
    m_nlsp.SetInfoHandlerFn([](const char *module, const char *function, char *msg){
        std::cout << "[" << module << "] " << function << "\n   " << msg << std::endl;
    });
    m_nlsp.SetVerbosityLevel(1);
    m_nlsp.setInterIterationMonitorFunc([&nsww = m_nsww, &time](const double * x){
        nsww.setX(x);
        std::cout << "time = " << time.elapsed() << std::endl;
        nsww.RenderIteration();
        return 0;
    });
    m_nlsp.SetParameterInt("MaxSetupCalls", 1);
    m_nlsp.SetParameterInt("MaxSubSetupCalls", 1);
    m_nlsp.SetSolveStrategy(NonLinearSolverKinsol::LINESEARCH);
    m_nlsp.SetMaxBetaFails(100);
    m_nlsp.SetMaxNewtonStep(R / 10 * sqrt(nlp.m_dofs));

    m_nlsp.SetInitialGuess(m_nsww.getCurrentX());
}
void VirtualSuturerCil::SolverData::setContext(const SolverCtx& ctx){
    m_ls.SetParameterReal("drop_tolerance", ctx.lin_droptol);
    m_ls.SetParameterReal("reuse_tolerance", ctx.lin_reusetol);
    m_ls.setVerbosityLevel(ctx.verbose_level-1);
    m_nlsp.SetNumMaxIters(ctx.maxits);
    m_nlsp.SetScaledStepTol(ctx.nonlin_scaled_steptol);
    m_nlsp.SetFuncNormTol(ctx.nonlin_ftol * abs_err_init);
    m_nlsp.SetVerbosityLevel(ctx.verbose_level);
    SUNLinSolContent_Com(m_nlsp.LS)->verbosity = ctx.verbose_level;
    if (ctx.verbose_level >= 1){
        m_nlsp.SetInfoHandlerFn([](const char *module, const char *function, char *msg){
            std::cout << "[" << module << "] " << function << "\n   " << msg << std::endl;
        });
        m_nlsp.setInterIterationMonitorFunc([&nsww = m_nsww, time = m_time](const double * x){
            nsww.setX(x);
            if (time) std::cout << "time = " << time->elapsed() << std::endl;
            nsww.RenderIteration();
            return 0;
        });
    } else {
        m_nlsp.SetInfoHandlerFn([](const char *module, const char *function, char *msg){});
        m_nlsp.setInterIterationMonitorFunc([&nsww = m_nsww](const double * x){
            nsww.setX(x);
            nsww.RenderIteration();
            return 0;
        });
    }
}

VirtualSuturerCil::FreeMembraneCxt::FreeMembraneCxt(){
    m_damps.reserve(3);
    m_damps.emplace_back(1.7);
    m_damps.emplace_back(1e-2);
    m_damps.emplace_back(1e-4);
    for (int i = 0; i < 3; ++i) m_pres[i].cP = 1;
}

VirtualSuturerCil::ContactSutureCxt::ContactSutureCxt(){ 
    m_solve_ctx.maxits = 20; 
    m_contact.resize(4);
    for (int i = 0; i < 3; ++i){
        m_contact[0][i].cP = 0.1;
        m_contact[0][i].cHt = 0.5;
        m_contact[1][i].cP = 0.1;
        m_contact[1][i].cHt = 0.25;
        m_contact[2][i].cP = 0.5;
        m_contact[2][i].cHt = 0.25;
        m_contact[3][i].cP = 1;
        m_contact[3][i].cHt = 0.25;
    }
    
}

VirtualSuturerCil::SutureAlgoParams::SutureAlgoParams(){
    m_fsc.resize(2);
    for (int i = 0; i < 3; ++i) m_fsc[0].m_pres[i].cP = 1;
    m_fsc[0].m_damps.csigma = 4e-2;
    for (int i = 0; i < 3; ++i) m_fsc[1].m_pres[i].cP = 1;
    m_fsc[1].m_damps.csigma = 1e-4;
}

void VirtualSuturerCil::CommissureGeom::setup(){
    auto& A = m_A;
    m_n = CGAL::cross_product(A[1] - A[0], A[2] - A[0]);
    m_area = sqrt(m_n.squared_length()) / 2;
    m_n /= (2*m_area);
    std::array<DReal, 3> a;
    for (int i = 0; i < 3; ++i) a[i] = (A[(i+1)%3] - A[(i+2)%3]).squared_length();
    std::array<DReal, 3> l;
    DReal s = 0;
    for (int i = 0; i < 3; ++i){
        l[i] = a[i] * (a[(i+1)%3] + a[(i+2)%3] - a[i]);
        s += l[i];
    }
    for (int i = 0; i < 3; ++i) m_a[i] = sqrt(a[i]);
    for (int i = 0; i < 3; ++i) l[i] /= s;
    m_O = CGAL::ORIGIN + l[0] * (A[0]-CGAL::ORIGIN) + l[1] * (A[1]-CGAL::ORIGIN) + l[2] * (A[2]-CGAL::ORIGIN);
    m_M = CGAL::ORIGIN + ((A[0]-CGAL::ORIGIN) + (A[1]-CGAL::ORIGIN) + (A[2]-CGAL::ORIGIN))/3;
    m_R = sqrt((A[0] - m_O).squared_length());
}

SignedDistanceField::SDF VirtualSuturerCil::CilDist::operator()(const Vector& x) const{
    SDF d;
    Vector rv = (x - C) - ((x - C)*n)*n;
    DReal r = sqrt(rv.squared_length());
    d.sdf = R - r;
    d.grad_sdf = -rv / r;
    Vector f = CGAL::cross_product(n, d.grad_sdf);
    for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j)
        d.grad_grad_sdf(i, j) = f[i]*f[j] / r;
    return d;
}

bool VirtualSuturerCil::LeafBoundaries::addShift(DReal dt){
    bool res = ((1.0 - sutureT) <= dt);
    if (res) dt = 1.0 - sutureT;
    sutureT += dt;
    for (int i = 0; i < v_suture_bnd.size(); ++i)
        obj->m_x[v_suture_bnd[i]] += dt*suture_shift[i];
    return res;
}

VirtualSuturerCil::VirtualSuturerCil() { m_commissure.m_A.reserve(3); }
VirtualSuturerCil& VirtualSuturerCil::setCommissurePoints(std::array<Point, 3> p){ 
    m_commissure.m_A.resize(3); 
    std::copy(p.begin(), p.end(), m_commissure.m_A.begin()); 
    return *this;
}
VirtualSuturerCil& VirtualSuturerCil::setLeafLengths(std::array<DReal, 3> l){
    for (int i = 0; i < 3; ++i) m_linfo[i].m_l = l[i];
    return *this;
}
VirtualSuturerCil& VirtualSuturerCil::setViewCtx(bool with_view, int argc, char** argv){
    m_vw_ctx.with_view = with_view;
    if (argc > 0) m_vw_ctx.argc = argc;
    if (argv != nullptr) m_vw_ctx.argv = argv;
    return *this;
}
VirtualSuturerCil& VirtualSuturerCil::setAorticSDF(std::shared_ptr<SignedDistanceField> aorta){ 
    m_aorta = std::move(aorta); 
    return *this; 
}
VirtualSuturerCil& VirtualSuturerCil::setCilindricAortic(){
    m_use_cilindric_aorta = true;
    return *this;
}
VirtualSuturerCil::LeafInfo& VirtualSuturerCil::pushLeaf(const Mesh& m, Curve suture_line, Field clamp_suture_b){
    if (m_nleafs >= 3) throw std::runtime_error("Can't insert new leaflet");
    m_linfo[m_nleafs].m_id = m_nleafs;
    m_objs.push_back(Object3D(m));
    //m_linfo[m_nleafs].m_id = m_w.addObject3D(Object3D(m), 1);
    if (suture_line) m_linfo[m_nleafs].m_suture_line = suture_line;
    if (clamp_suture_b) m_linfo[m_nleafs].m_clamp_suture_b = clamp_suture_b;
    return m_linfo[m_nleafs++];
}
VirtualSuturerCil& VirtualSuturerCil::setSaveDirectory(std::string dir, std::string prefix, VirtualSuturerCil::SaveCxt::StepSign stat){
    m_saver.save_directory = dir;
    m_saver.prefix = prefix;
    m_saver.sign_stat = stat;
    return *this;
}
void VirtualSuturerCil::setup(){
    if (m_commissure.m_A.size() != 3) {
        for (int i = m_commissure.m_A.size(); i < 3; ++i){
            auto &l0 = m_linfo[(i+1)%3].m_suture_line, &l1 = m_linfo[(i+2)%3].m_suture_line;
            if (l0 && l1) m_commissure.m_A.push_back(CGAL::ORIGIN + ((l0(1) - CGAL::ORIGIN) + (l1(0) - CGAL::ORIGIN))/2); 
            else throw std::runtime_error("3 commissure points must be defined");
        }
    }
    m_commissure.setup();
    Vector n = m_commissure.m_n;
    Point C = m_commissure.m_O;
    DReal R = m_commissure.m_R;
    for (int i = 0; i < 3; ++i){
        if (std::isnan(m_linfo[i].m_l)){
            if (m_linfo[i].m_id < 0) throw std::runtime_error("3 leaf lengths must be defined");
            else {
                //auto& o = m_w.obj(m_linfo[i].m_id);
                auto& o = m_objs[m_linfo[i].m_id];
                Point p[2];
                for (auto v: o.m_mesh.vertices()){
                    if (o.m_boundary[v] == (2 | 1) )
                        p[0] = o.m_x[v];
                    else if (o.m_boundary[v] == (2 | 4) )
                        p[1] = o.m_x[v];
                }
                m_linfo[i].m_l = sqrt((p[1] - p[0]).squared_length());
            }
        }
        if (!(m_linfo[i].m_id < 0)){
            if (!m_linfo[i].m_suture_line)
                throw std::runtime_error("Suture line for " + to_string(i) + "th leaf is not defined");
            if (std::isnan(m_linfo[i].m_P) || std::isnan(m_linfo[i].m_Ht) || std::isnan(m_linfo[i].m_E)){
                if (std::isnan(m_linfo[i].m_Ps)) throw std::runtime_error("Wrong DlessPressure for " + to_string(i) + "th leaf");
                if (std::isnan(m_linfo[i].m_P)){
                    if (std::isnan(m_linfo[i].m_Ht)) m_linfo[i].m_Ht = 0.5;
                    if (std::isnan(m_linfo[i].m_E)) m_linfo[i].m_E = 1e3;
                    m_linfo[i].m_P = m_linfo[i].m_Ps * m_linfo[i].m_E * m_linfo[i].m_Ht / m_linfo[i].m_l;
                } else {
                    if (std::isnan(m_linfo[i].m_E)){
                        if (std::isnan(m_linfo[i].m_Ht)) m_linfo[i].m_Ht = 0.5;
                        m_linfo[i].m_E = m_linfo[i].m_P * m_linfo[i].m_l / (m_linfo[i].m_Ps * m_linfo[i].m_Ht);
                    } else{
                        m_linfo[i].m_Ht = m_linfo[i].m_P * m_linfo[i].m_l / (m_linfo[i].m_Ps * m_linfo[i].m_E);
                    }
                }
            }
            if (m_linfo[i].m_use_cilidric_clamp && !m_linfo[i].m_clamp_suture_b){
                m_linfo[i].m_clamp_suture_b = [s = getCilinderSuturer(C, n), l = m_linfo[i].m_suture_line](double t){
                    Vector l0, l1;
                    Point lp = l(t);
                    if (t - 1e-7 >= 0) l0 = l(t - 1e-7) - CGAL::ORIGIN;
                    else l0 = lp - CGAL::ORIGIN;
                    if (t + 1e-7 <= 1) l1 = l(t + 1e-7) - CGAL::ORIGIN;
                    else l1 = lp - CGAL::ORIGIN;
                    Vector dir = l1 - l0;
                    dir /= sqrt(dir.squared_length());
                    return s(lp, dir);
                };
            }
            if (!m_linfo[i].m_smap) m_linfo[i].m_smap = SutureMapCenterConstant();
        }
    }
    if (m_vw_ctx.with_view && (m_vw_ctx.argc < 1 || m_vw_ctx.argv == nullptr))
        throw std::runtime_error("Window's context is not complete");
    
    if (m_use_cilindric_aorta && !m_aorta){
        m_aorta = std::make_shared<CilDist>(C-CGAL::ORIGIN, n, R);
    }
}
VirtualSuturerCil::LeafBoundaries VirtualSuturerCil::geometrical_shifting_leaf(int ileaf){
    if (ileaf >= m_nleafs) throw std::runtime_error("Access to uncreated leaflet");
    initial_shift_and_rotate_leaf(ileaf);
    auto res = map_cusp_to_aortic_suture(m_objs[m_linfo[ileaf].m_id], m_linfo[ileaf].m_suture_line, m_linfo[ileaf].m_smap, m_linfo[ileaf].m_clamp_suture_b);
    set_bdata_tag(ileaf, res);
    set_dirichlet_bc_leaf(ileaf);
    return res;
}
bool VirtualSuturerCil::geometrical_shifting(int save_step, std::array<LeafBoundaries, 3>* pbnds){
    if (m_nleafs <= 0) return true;
    initial_shift_and_rotate();
    auto bnds = map_cusps_to_suture_lines();
    for (int i = 0; i < m_nleafs; ++i) set_bdata_tag(i, bnds[i]);
    set_dirichlet_bc();
    if (save_step >= 0) save_leafs(save_step, "initial_shift_and_rotate");
    if (pbnds) *pbnds = std::move(bnds);
    return true;
}
bool VirtualSuturerCil::suture_leaf(int ileaf){
    if (ileaf >= m_nleafs) throw std::runtime_error("Access to uncreated leaflet");
    auto bnd = geometrical_shifting_leaf(ileaf);
    double default_sigma = m_linfo[ileaf].m_E*m_linfo[ileaf].m_Ht * m_commissure.m_R;

    World lw;
    auto lid = lw.addObject3D(std::move(m_objs[ileaf]), 1);
    Object3D& obj = lw.obj(lid);
    bnd.obj = &obj;
    start_renderer(lw, m_vw_ctx, front_end);
    m_solver.initialize(lw, m_commissure.m_R, m_sym_time);
    initialize_collider(lw);
    set_collission_object_thickness(lw, lid, ileaf);
    
    set_free_membrane_forces(obj, m_p.m_fmc, ileaf);
    apply_free_membrane_suture(lw, &bnd, 1, m_p.m_fmc, m_solver, m_p.m_forces.damp, default_sigma);
    m_saver.save_step(obj, ileaf, 1, "free_memrane_suture");
    reinterpret_cast<BridsonCollisionManager*>(lw.getCollider())->set_MaxRelaxIts(5);

    for (auto& ctx: m_p.m_fsc){
        set_free_shell_forces(obj, ctx, ileaf, bnd);
        apply_free_shell_suture(lw, ctx, m_solver, m_p.m_forces.damp, default_sigma);
    }
    if (!m_p.m_fsc.empty()) m_saver.save_step(obj, ileaf, 2, "free_shell_suture");

    Object3D* pobj = &obj;
    auto tmp = m_p.m_csc.contact_dt;
    for (int k = 0; k < m_p.m_csc.m_contact.size(); ++k){
        set_contact_suture_forces(obj, m_p.m_csc, ileaf, bnd, k);
        apply_contact_suture(lw, &pobj, &m_linfo[ileaf].m_cst, &m_p.m_forces.sf[ileaf].contact_sdf, &m_p.m_csc.contact_dt[ileaf], 1, m_p.m_csc,
                            m_solver, m_p.m_forces.damp, default_sigma, k==0);
        std::fill(m_p.m_csc.contact_dt.begin(), m_p.m_csc.contact_dt.end(), 1.0);                    
    }
    m_p.m_csc.contact_dt = tmp;
    if (!m_p.m_csc.m_contact.empty()) m_saver.save_step(obj, ileaf, 3, "contact_surface_suture");

    if (m_aorta){
        set_aortic_constraint_forces(obj, m_p.m_asc, ileaf, bnd);
        apply_aortic_constraint(lw, m_p.m_asc, m_solver, m_p.m_forces.damp, default_sigma);
        m_saver.save_step(obj, ileaf, 4, "aortic_surface_suture");
    }

    for (auto& c: m_p.m_custom){
        set_custom_forces(obj, c, ileaf, bnd);
        apply_custom_suture(lw, c, m_solver, m_p.m_forces.damp, default_sigma);
    }
    if (!m_p.m_custom.empty()) m_saver.save_step(obj, ileaf, 5, "aortic_custom_suture");

    stop_renderer(lw, m_vw_ctx, front_end);
    m_objs[ileaf] = lw.releaseObj(lid);

    return true;
}

// bool VirtualSuturerCil::suture(){
//     if (m_nleafs <= 0) return true;
//     setup();
//     std::array<LeafBoundaries, 3> bnds;
//     geometrical_shifting(0, &bnds);
//     set_free_membrane_forces();
//     initialize_collider(m_w);
//     for (int i = 0; i < 3; ++i) set_collission_object_thickness(m_w, m_linfo[i].m_id, i);
//     m_solver.initialize(m_w, m_commissure.m_R, m_sym_time);
//     start_renderer(m_w, m_vw_ctx, front_end);

//     apply_free_membrane_suture(bnds);
//     save_leafs(1, "free_memrane_suture");
//     reinterpret_cast<BridsonCollisionManager*>(m_w.getCollider())->set_MaxRelaxIts(5);
//     apply_free_shell_suture(bnds);
//     save_leafs(2, "free_shell_suture");
//     apply_contact_suture(bnds);
//     save_leafs(3, "contact_surface_suture");
//     apply_aortic_constraint(bnds);
//     if (m_aorta) save_leafs(4, "aortic_surface_suture");

//     stop_renderer(m_w, m_vw_ctx, front_end);
//     return true;
// }
void VirtualSuturerCil::save_leafs(int step, std::string step_name){
    for (int i = 0; i < m_nleafs; ++i)
        m_saver.save_step(m_objs[m_linfo[i].m_id], i, step, step_name);
}
// void VirtualSuturerCil::apply_aortic_constraint(std::array<LeafBoundaries, 3>& bnds){
//     if (!m_aorta) return;
//     apply_aortic_constraint_forces(m_p.m_asc, bnds);
//     solve_problem(m_p.m_asc.m_solve_ctx);
// }
// void VirtualSuturerCil::apply_contact_suture(std::array<LeafBoundaries, 3>& bnds){
//     apply_contact_suture_forces(m_p.m_csc, bnds);
//     bool is_sutured = false;
//     std::array<DReal,3> lambdas = {0, 0, 0}; 
//     for (int i = 0; i < m_nleafs; ++i) set_sdf_lambda(i, lambdas[i]);
//     solve_problem(m_p.m_csc.m_solve_ctx);
//     do{
//         is_sutured = true;
//         for (int i = 0; i < m_nleafs; ++i) {
//             if (lambdas[i] >= 1) continue;
//             lambdas[i] += m_p.m_csc.contact_dt[i];
//             if (lambdas[i] >= 1.0) {
//                 lambdas[i] = 1.0;
//             } else is_sutured &= false;
//             set_sdf_lambda(i, lambdas[i]);
//         }
//         solve_problem(m_p.m_csc.m_solve_ctx);
//     } while (!is_sutured);
// }
// void VirtualSuturerCil::apply_aortic_constraint_forces(const AorticSutureCxt& cxt, std::array<LeafBoundaries, 3>& bnds){
//     for (int nobj = 0; nobj < m_nleafs; ++nobj){
//         set_elastic_force(nobj, cxt.m_elast[nobj]);
//         set_bending_force(nobj, cxt.m_bend[nobj]);
//         set_clamped_bc(nobj, bnds[nobj].clamp_data);
//         set_pressure_force(nobj, cxt.m_pres[nobj]);
//         set_free_edge_force(nobj, cxt.m_free_edge[nobj]);
//         set_sdf_force(nobj, cxt.m_contact[nobj]);
//         set_aortic_force(nobj, cxt.m_aortic[nobj]);
//     }
//     set_damp_force(cxt.m_damps);
// }
// void VirtualSuturerCil::apply_contact_suture_forces(ContactSutureCxt& cxt, std::array<LeafBoundaries, 3>& bnds){
//     for (int nobj = 0; nobj < m_nleafs; ++nobj){
//         set_elastic_force(nobj, cxt.m_elast[nobj]);
//         set_bending_force(nobj, cxt.m_bend[nobj]);
//         set_clamped_bc(nobj, bnds[nobj].clamp_data);
//         set_pressure_force(nobj, cxt.m_pres[nobj]);
//         set_free_edge_force(nobj, cxt.m_free_edge[nobj]);
//         set_sdf_force(nobj, cxt.m_contact[nobj], 0.0);
//     }
//     set_damp_force(cxt.m_damps);
// }
// void VirtualSuturerCil::apply_free_shell_suture(std::array<LeafBoundaries, 3>& bnds){
//     for (auto& ctx: m_p.m_fsc){
//         apply_free_shell_forces(ctx, bnds);
//         solve_problem(ctx.m_solve_ctx);
//     }
// }
bool VirtualSuturerCil::set_collission_object_thickness(World& w, ObjectID w_id, int ileaf){
    BridsonCollisionManager* pcol = reinterpret_cast<BridsonCollisionManager*>(w.getCollider());
    if (!pcol) return false;
    DReal Ht = m_p.m_collide.Ht;
    if (Ht < 0) Ht = m_p.m_collide.cHt * m_linfo[ileaf].m_Ht;

    CollisionThicknessCompressed<CollisionThicknessConst> cthickness{CollisionThicknessConst(Ht)};
    CollisionThickness Htf(cthickness);
    pcol->set_thickness(w_id, Htf);
    return true;
}
void VirtualSuturerCil::initialize_collider(World& w){
    unique_ptr<CollisionManagerBase> colMan = std::make_unique<BridsonCollisionManager>();
    w.setCollider(std::move(colMan));
    BridsonCollisionManager* pcol = reinterpret_cast<BridsonCollisionManager*>(w.getCollider());
    // for (int i = 0; i < m_nleafs; ++i){
    //     DReal Ht = m_p.m_collide.Ht;
    //     if (Ht < 0) Ht = m_p.m_collide.cHt * m_linfo[i].m_Ht;
    //     pcol->set_thickness(m_linfo[i].m_id, Ht);
    // }
    DReal P = 0;
    for (int i = 0; i < m_nleafs; ++i)
        P += m_linfo[i].m_P;
    P /= m_nleafs;
    pcol->set_Kspring(P * m_p.m_collide.cPspr);
    pcol->set_Dt(m_p.m_collide.Dt);
    pcol->setCollisionStatus(BridsonCollisionManager::SelfCollision);
    pcol->setAsInitialNonCollidingState();
    pcol->set_MaxRelaxIts(0);
    pcol->set_RelaxEps(0.5);
}
void VirtualSuturerCil::set_free_shell_forces(Object3D& obj, const FreeShellCxt& ctx, int ileaf, LeafBoundaries& bnd){
    set_elastic_force(obj, m_linfo[ileaf], ctx.m_elast[ileaf], m_p.m_forces.sf[ileaf].elast, m_gen_force_dir, m_regenerate_elastic_force);
    set_bending_force(obj, m_linfo[ileaf], ctx.m_bend[ileaf], m_p.m_forces.sf[ileaf].bend, m_gen_force_dir, m_regenerate_elastic_force, m_regenerate_bending_force);
    if (ctx.m_bend[ileaf].use_clamped_bc) set_clamped_bc(obj, m_p.m_forces.sf[ileaf].bend, bnd.clamp_data);
    set_pressure_force(obj, m_linfo[ileaf], ctx.m_pres[ileaf], m_p.m_forces.sf[ileaf].pres);
    set_free_edge_force(obj, m_linfo[ileaf], ctx.m_free_edge[ileaf], m_p.m_forces.sf[ileaf].free_edge, m_commissure.m_n, m_commissure.m_O);
}
void VirtualSuturerCil::set_contact_suture_forces(Object3D& obj, const ContactSutureCxt& ctx, int ileaf, LeafBoundaries& bnd, int ncontact){
    set_elastic_force(obj, m_linfo[ileaf], ctx.m_elast[ileaf], m_p.m_forces.sf[ileaf].elast, m_gen_force_dir, m_regenerate_elastic_force);
    set_bending_force(obj, m_linfo[ileaf], ctx.m_bend[ileaf], m_p.m_forces.sf[ileaf].bend, m_gen_force_dir, m_regenerate_elastic_force, m_regenerate_bending_force);
    if (ctx.m_bend[ileaf].use_clamped_bc) set_clamped_bc(obj, m_p.m_forces.sf[ileaf].bend, bnd.clamp_data);
    set_pressure_force(obj, m_linfo[ileaf], ctx.m_pres[ileaf], m_p.m_forces.sf[ileaf].pres);
    set_free_edge_force(obj, m_linfo[ileaf], ctx.m_free_edge[ileaf], m_p.m_forces.sf[ileaf].free_edge, m_commissure.m_n, m_commissure.m_O);
    set_sdf_force(obj, m_linfo, ileaf, std::array<Point, 3>{m_commissure.m_A[0], m_commissure.m_A[1], m_commissure.m_A[2]}, ctx.m_contact[ncontact][ileaf], m_p.m_forces.sf[ileaf].contact_sdf, 0.0);
}
void VirtualSuturerCil::set_aortic_constraint_forces(Object3D& obj, const AorticSutureCxt& ctx, int ileaf, LeafBoundaries& bnd){
    set_elastic_force(obj, m_linfo[ileaf], ctx.m_elast[ileaf], m_p.m_forces.sf[ileaf].elast, m_gen_force_dir, m_regenerate_elastic_force);
    set_bending_force(obj, m_linfo[ileaf], ctx.m_bend[ileaf], m_p.m_forces.sf[ileaf].bend, m_gen_force_dir, m_regenerate_elastic_force, m_regenerate_bending_force);
    if (ctx.m_bend[ileaf].use_clamped_bc) set_clamped_bc(obj, m_p.m_forces.sf[ileaf].bend, bnd.clamp_data);
    set_pressure_force(obj, m_linfo[ileaf], ctx.m_pres[ileaf], m_p.m_forces.sf[ileaf].pres);
    set_free_edge_force(obj, m_linfo[ileaf], ctx.m_free_edge[ileaf], m_p.m_forces.sf[ileaf].free_edge, m_commissure.m_n, m_commissure.m_O);
    set_sdf_force(obj, m_linfo, ileaf, std::array<Point, 3>{m_commissure.m_A[0], m_commissure.m_A[1], m_commissure.m_A[2]}, ctx.m_contact[ileaf], m_p.m_forces.sf[ileaf].contact_sdf, 1.0);
    set_aortic_force(obj, m_linfo[ileaf], ctx.m_aortic[ileaf], m_p.m_forces.sf[ileaf].aorta_sdf, m_aorta);
}
static void _remove_obj_force(Object3D& obj, Force_ID& id){
    if (id >= 0){
        obj.remove_force(id);
        id = -1;
    }
}
void VirtualSuturerCil::set_custom_forces(Object3D& obj, const CustomCxt& ctx, int ileaf, LeafBoundaries& bnd){
    auto _remove_obj_force = [](Object3D& obj, Force_ID& id){
        if (id >= 0){
            obj.remove_force(id);
            id = -1;
        }
    };
    if (ctx.m_elast.second) set_elastic_force(obj, m_linfo[ileaf], ctx.m_elast.first[ileaf], m_p.m_forces.sf[ileaf].elast, m_gen_force_dir, m_regenerate_elastic_force);
    else _remove_obj_force(obj, m_p.m_forces.sf[ileaf].elast);
    
    if (ctx.m_bend.second) {
        set_bending_force(obj, m_linfo[ileaf], ctx.m_bend.first[ileaf], m_p.m_forces.sf[ileaf].bend, m_gen_force_dir, m_regenerate_elastic_force, m_regenerate_bending_force);
        if (ctx.m_bend.first[ileaf].use_clamped_bc) set_clamped_bc(obj, m_p.m_forces.sf[ileaf].bend, bnd.clamp_data);
    } else _remove_obj_force(obj, m_p.m_forces.sf[ileaf].bend);
    if (ctx.m_pres.second) set_pressure_force(obj, m_linfo[ileaf], ctx.m_pres.first[ileaf], m_p.m_forces.sf[ileaf].pres);
    else _remove_obj_force(obj, m_p.m_forces.sf[ileaf].pres);
    if (ctx.m_free_edge.second) set_free_edge_force(obj, m_linfo[ileaf], ctx.m_free_edge.first[ileaf], m_p.m_forces.sf[ileaf].free_edge, m_commissure.m_n, m_commissure.m_O);
    else _remove_obj_force(obj, m_p.m_forces.sf[ileaf].free_edge);
    if (ctx.m_contact.second) set_sdf_force(obj, m_linfo, ileaf, std::array<Point, 3>{m_commissure.m_A[0], m_commissure.m_A[1], m_commissure.m_A[2]}, ctx.m_contact.first[ileaf], m_p.m_forces.sf[ileaf].contact_sdf, 1.0);
    else _remove_obj_force(obj, m_p.m_forces.sf[ileaf].contact_sdf);
    if (ctx.m_aortic.second) set_aortic_force(obj, m_linfo[ileaf], ctx.m_aortic.first[ileaf], m_p.m_forces.sf[ileaf].aorta_sdf, m_aorta);
    else _remove_obj_force(obj, m_p.m_forces.sf[ileaf].aorta_sdf);
}
// void VirtualSuturerCil::apply_free_shell_forces(FreeShellCxt& cxt, std::array<LeafBoundaries, 3>& bnd){
//     for (int nobj = 0; nobj < m_nleafs; ++nobj){
//         set_elastic_force(nobj, cxt.m_elast[nobj]);
//         set_bending_force(nobj, cxt.m_bend[nobj]);
//         set_clamped_bc(nobj, bnd[nobj].clamp_data);
//         set_pressure_force(nobj, cxt.m_pres[nobj]);
//         set_free_edge_force(nobj, cxt.m_free_edge[nobj]);
//     }
//     set_damp_force(cxt.m_damps);
// }
void VirtualSuturerCil::apply_free_membrane_suture(World& w, LeafBoundaries* bnd, int nbnds, const FreeMembraneCxt& ctx,
                                                    SolverData& s, WorldForce_ID& fid, double default_sigma){
    bool is_sutured = (nbnds <= 0);
    do{
        bool lis_sutured = true;
        for (int i = 0; i < nbnds; ++i){
            lis_sutured &= bnd[i].addShift(ctx.suture_dt[i]);
        }
        reinterpret_cast<BridsonCollisionManager*>(w.getCollider())->setAsInitialNonCollidingState();
        is_sutured = lis_sutured;
        int nsub_it = ctx.m_damps.size();
        if (nsub_it < 1) nsub_it = 1;

        for (int n = 0; n < nsub_it; ++n){
            if (ctx.m_damps.empty()) {
                set_damp_force(w, DampForceCtx(0), fid, default_sigma);
            } else {
                set_damp_force(w, ctx.m_damps[n], fid, default_sigma);
            }
            solve_problem(s, ctx.m_solve_ctx, fid);
        }
    } while(!is_sutured);
}
void VirtualSuturerCil::apply_free_shell_suture(World& w, const FreeShellCxt& ctx,
                             SolverData& s, WorldForce_ID& fid, double default_sigma){
    set_damp_force(w, ctx.m_damps, fid, default_sigma);
    solve_problem(s, ctx.m_solve_ctx, fid);
}
void VirtualSuturerCil::apply_contact_suture(World& w, 
                            Object3D* objs[], LeafInfo::ContactSurfType tp[], Force_ID contact_sdf[], DReal contact_dt[], int nobjs, 
                            const ContactSutureCxt& ctx, SolverData& s, WorldForce_ID& fid, double default_sigma, bool zero_solve){
    set_damp_force(w, ctx.m_damps, fid, default_sigma);
    bool is_sutured = false;
    std::array<DReal, 3> lambdas = {0};
    for (int i = nobjs; i < 3; ++i) lambdas[i] = 1;
    
    if (zero_solve){
        for (int i = 0; i < nobjs; ++i) set_sdf_lambda(*(objs[i]), tp[i], contact_sdf[i], lambdas[i]);
        solve_problem(s, ctx.m_solve_ctx, fid);
    }
    //solve_problem(m_p.m_csc.m_solve_ctx);
    do{
        is_sutured = true;
        for (int i = 0; i < nobjs; ++i) {
            if (lambdas[i] >= 1) continue;
            lambdas[i] += contact_dt[i];
            if (lambdas[i] >= 1.0) {
                lambdas[i] = 1.0;
            } else is_sutured &= false;
            set_sdf_lambda(*(objs[i]), tp[i], contact_sdf[i], lambdas[i]);
        }
        solve_problem(s, ctx.m_solve_ctx, fid);
        //solve_problem(m_p.m_csc.m_solve_ctx);
    } while (!is_sutured);
    //solve_problem(s, ctx.m_solve_ctx, fid);
}
void VirtualSuturerCil::apply_aortic_constraint(World& w, const AorticSutureCxt& ctx,
                             SolverData& s, WorldForce_ID& fid, double default_sigma){
    set_damp_force(w, ctx.m_damps, fid, default_sigma);
    solve_problem(s, ctx.m_solve_ctx, fid);
}
void VirtualSuturerCil::apply_custom_suture(World& w, const CustomCxt& ctx,
                        SolverData& s, WorldForce_ID& fid, double default_sigma){
    if (ctx.m_damps.second) set_damp_force(w, ctx.m_damps.first, fid, default_sigma);
    else {
        w.removeForce(fid);
        fid = WorldForce_ID();
    }
    solve_problem(s, ctx.m_solve_ctx, fid);
}

// void VirtualSuturerCil::apply_free_membrane_suture(std::array<LeafBoundaries, 3>& bnd){
//     bool is_sutured = (m_nleafs <= 0);
//     do{
//         bool lis_sutured = true;
//         for (int i = 0; i < m_nleafs; ++i){
//             lis_sutured &= bnd[i].addShift(m_p.m_fmc.suture_dt[i]);
//         }
//         reinterpret_cast<BridsonCollisionManager*>(m_w.getCollider())->setAsInitialNonCollidingState();
//         is_sutured = lis_sutured;
//         int nsub_it = m_p.m_fmc.m_damps.size();
//         if (nsub_it < 1) nsub_it = 1;

//         for (int n = 0; n < nsub_it; ++n){
//             if (m_p.m_fmc.m_damps.empty()) {
//                 set_damp_force(DampForceCtx(0));
//             } else {
//                 set_damp_force(m_p.m_fmc.m_damps[n]);
//             }
//             solve_problem(m_p.m_fmc.m_solve_ctx);
//         }
//     } while(!is_sutured);
// }
void VirtualSuturerCil::solve_collisions(SolverData& s){
    auto pcol = reinterpret_cast<BridsonCollisionManager*>(s.m_nsww.pw->getCollider());
    if (pcol){
        pcol->findCollisions();
        pcol->solveCollisions();
    }
    s.m_nlsp.SetInitialGuess(s.m_nsww.getCurrentX());
}
bool VirtualSuturerCil::solve_problem(SolverData& s, const SolverCtx& ctx, WorldForce_ID damp_id){
    s.setContext(ctx);
    bool slvFlag = s.m_nlsp.Solve();
    solve_collisions(s);
    s.m_nsww.RenderIteration();
    if (damp_id.id >= 0)
        s.m_nsww.pw->force(damp_id).target<WVertDampForce>()->resetCurrentState();
            
    return slvFlag;
}
bool VirtualSuturerCil::solve_problem(SolverCtx ctx){
    return solve_problem(m_solver, ctx, m_p.m_forces.damp);
}
void VirtualSuturerCil::set_free_membrane_forces(){
    for (int i = 0; i < m_nleafs; ++i)
        set_free_membrane_forces(i);
}
void VirtualSuturerCil::set_sdf_lambda(Object3D& obj, LeafInfo::ContactSurfType tp, Force_ID& id, DReal lambda){
    if (id < 0) return;
    auto sdf_force = obj.m_forces[id].target<SDFForce>();
    auto sdf = sdf_force->field<>();
    switch (tp){
        case LeafInfo::PLANE:{
            auto psdf = static_cast<CilindricContactSurface::LeafPlaneContactSDF*>(sdf);
            psdf->set_lambda(lambda);
            break;
        }
        case LeafInfo::BEZIER_CILINDRIC:{
            auto psdf = static_cast<CilindricContactSurface::LeafContactSDF*>(sdf);
            psdf->set_lambda(lambda);
            break;
        }
        default: 
            throw std::runtime_error("Faced unknown contact surface type");
    }
}
void VirtualSuturerCil::set_sdf_lambda(int nobj, DReal lambda){
    set_sdf_lambda(m_objs[m_linfo[nobj].m_id], m_linfo[nobj].m_cst, m_p.m_forces.sf[nobj].contact_sdf, lambda);
}
void VirtualSuturerCil::set_aortic_force(Object3D& obj, const LeafInfo& linfo, const AorticForceCtx& ctx, Force_ID& id, const std::shared_ptr<SignedDistanceField>& aorta){
    if (id >= 0 && ctx.cP == 0){
        obj.remove_force(id);
        id = -1;
        return;
    }
    if (ctx.cP == 0) return;
    DReal Ht = ctx.Ht;
    if (Ht < 0) Ht = ctx.cHt * linfo.m_Ht;
    if (id < 0){
        Force sdf_f = SDFForce(aorta, ctx.cP * linfo.m_P, Ht);
        id = obj.add_force(std::move(sdf_f));
    } else {
        obj.m_forces[id].target<SDFForce>()->set_dist_funcs(ctx.cP * linfo.m_P, Ht);
    }
}
void VirtualSuturerCil::set_aortic_force(int nobj, const AorticForceCtx& ctx){
    set_aortic_force(m_objs[m_linfo[nobj].m_id], m_linfo[nobj], ctx,  m_p.m_forces.sf[nobj].aorta_sdf, m_aorta);
}
void VirtualSuturerCil::set_sdf_force(Object3D& obj, const std::array<LeafInfo, 3>& linfo, int ileaf, const std::array<Point, 3>& _A, const SDFForceCtx& ctx, Force_ID& id, DReal lambda){
    std::array<double, 3> init_leaf_lengths = {linfo[0].m_l, linfo[1].m_l, linfo[2].m_l};
    if (id >= 0 && ctx.cP == 0){
        obj.remove_force(id);
        id = -1;
        return;
    }
    if (ctx.cP == 0) return;
    DReal Ht = ctx.Ht;
    if (Ht < 0) Ht = ctx.cHt * linfo[ileaf].m_Ht;
    if (id < 0){
        CilindricContactSurface ccs(_A, init_leaf_lengths);
        std::shared_ptr<SignedDistanceField> contact_surf;
        switch (linfo[ileaf].m_cst){
            case LeafInfo::PLANE:{
                contact_surf = ccs.getLeafPlaneContactSdf(ileaf);
                static_cast<CilindricContactSurface::LeafPlaneContactSDF*>(contact_surf.get())->set_lambda(lambda);
                break;
            }
            case LeafInfo::BEZIER_CILINDRIC:{
                contact_surf = ccs.getLeafContactSdf(ileaf);
                static_cast<CilindricContactSurface::LeafContactSDF*>(contact_surf.get())->set_lambda(lambda);
                break;
            }
            default: 
                throw std::runtime_error("Faced unknown contact surface type");
        }
        Force sdf_f = SDFForce(contact_surf, ctx.cP * linfo[ileaf].m_P, Ht, Ht*ctx.cshift);
        id = obj.add_force(std::move(sdf_f));
    } else {
        obj.m_forces[id].target<SDFForce>()->set_dist_funcs(ctx.cP * linfo[ileaf].m_P, Ht, Ht*ctx.cshift);
        switch (linfo[ileaf].m_cst){
        case LeafInfo::PLANE:
            obj.m_forces[id].target<SDFForce>()->field<CilindricContactSurface::LeafPlaneContactSDF>()->set_lambda(lambda);
            break;
        case LeafInfo::BEZIER_CILINDRIC:
            obj.m_forces[id].target<SDFForce>()->field<CilindricContactSurface::LeafContactSDF>()->set_lambda(lambda);
            break;
        default: 
            throw std::runtime_error("Faced unknown contact surface type");
        }
    }
}
void VirtualSuturerCil::set_sdf_force(int nobj, const SDFForceCtx& ctx, DReal lambda){
    set_sdf_force(m_objs[m_linfo[nobj].m_id], m_linfo, nobj, std::array<Point, 3>{m_commissure.m_A[0], m_commissure.m_A[1], m_commissure.m_A[2]}, ctx, m_p.m_forces.sf[nobj].contact_sdf, lambda);
}

void VirtualSuturerCil::set_elastic_force(Object3D& obj, const LeafInfo& linfo, const ElasticForceCtx& ctx, Force_ID& id, const std::string& gendir, bool& regen_elast_force){
    if (id >= 0){
        obj.remove_force(id);
        id = -1;
    }
    if (ctx.cE != 0.0){
        Force elastic = SVKirchhoffNoPoissonModel(ctx.cE_Ht*linfo.m_Ht, 0.0, ctx.cE*linfo.m_E/2, gendir, regen_elast_force);
        elastic.target<HyperElasticForceBase>()->prepareJacobianFunction(regen_elast_force);
        regen_elast_force = false;
        id = obj.add_force(std::move(elastic));
    }
}
void VirtualSuturerCil::set_elastic_force(int nobj, const ElasticForceCtx& ctx){
    set_elastic_force(m_objs[m_linfo[nobj].m_id], m_linfo[nobj], ctx, m_p.m_forces.sf[nobj].elast, m_gen_force_dir, m_regenerate_elastic_force);
}
void VirtualSuturerCil::set_bdata_tag(int nobj, VirtualSuturerCil::LeafBoundaries& bnd){
    if (bnd.clamp_data.empty() || m_linfo[nobj].suture_b_tag_name.empty()) return;
    auto& leaf = m_objs[m_linfo[nobj].m_id];
    auto& bdata = bnd.clamp_data;

    auto bdit = leaf.m_mesh.add_property_map<V_ind, Vector>(m_linfo[nobj].suture_b_tag_name).first;
    for (auto v: leaf.m_mesh.vertices()) bdit[v] = CGAL::NULL_VECTOR;
    for (auto dat: bdata) {
        auto e = dat.first;
        auto v = vert_around(leaf.m_mesh, e);
        for (int i = 0; i < 2 && dat.second.btype == BendingForce::BoundaryData::CLAMPED; ++i)
            bdit[v[i]] += Vector(dat.second.normal[0], dat.second.normal[1], dat.second.normal[2]);
    }
    for (auto v: leaf.m_mesh.vertices()) {
        double norm = sqrt(bdit[v].squared_length());
        if (norm > 1e-7) bdit[v] /= norm;
    }
}
void VirtualSuturerCil::set_damp_force(World& w, const DampForceCtx& ctx, WorldForce_ID& id, double default_sigma){
    if (id.id >= 0){
        if (ctx.sigma == 0 || ctx.csigma == 0){
            w.removeForce(id);
            id = WorldForce_ID();
        } else {
            DReal sigma = ctx.sigma;
            if (sigma < 0) sigma = ctx.csigma * default_sigma;
            auto dfc = std::make_shared<WVertDampForce::DampFuncConst>([s = sigma](Object3D& obj, const ObjectID id, V_ind v){ return s; });
            w.force(id).target<WVertDampForce>()->setDampFunc(std::move(dfc));
        } 
    } else if (ctx.sigma != 0 || ctx.csigma != 0){
        DReal sigma = ctx.sigma;
        if (sigma < 0) sigma = ctx.csigma * default_sigma;
        auto dfc = std::make_shared<WVertDampForce::DampFuncConst>([s = sigma](Object3D& obj, const ObjectID id, V_ind v){ return s; });
        WVertDampForce Df;
        Df.setDampFunc(std::move(dfc));
        id = w.addForce(std::move(Df));
        w.force(id).target<WVertDampForce>()->registerWorld(w);
    }
}

void VirtualSuturerCil::set_damp_force(World& w, const DampForceCtx& ctx){
    set_damp_force(w, ctx, m_p.m_forces.damp, m_linfo[0].m_E*m_linfo[0].m_Ht * m_commissure.m_R);
}
void VirtualSuturerCil::set_bending_force(Object3D& obj, const LeafInfo& linfo, const BendForceCtx& ctx, Force_ID& id, 
                                        const std::string& gendir, bool& regen_elast_force, bool& regen_bend_force){
    if (id >= 0){
        obj.remove_force(id);
        id = -1;
    }
    if (ctx.cE != 0.0){
        Force proxy = SVKirchhoffNoPoissonModel(ctx.cE_Ht*linfo.m_Ht, 0.0, ctx.cE*linfo.m_E/2, gendir, regen_elast_force);
        Force bending = BendingForce(proxy.target<HyperElasticForceBase>()->f, regen_bend_force);
        bending.target<BendingForce>()->prepareJacobianFunction(regen_bend_force);
        regen_elast_force = regen_bend_force = false;
        id = obj.add_force(std::move(bending));
    }
}
void VirtualSuturerCil::set_bending_force(int nobj, const BendForceCtx& ctx){
    set_bending_force(m_objs[m_linfo[nobj].m_id], m_linfo[nobj], ctx, m_p.m_forces.sf[nobj].bend, m_gen_force_dir, m_regenerate_elastic_force, m_regenerate_bending_force);
}
void VirtualSuturerCil::set_clamped_bc(Object3D& obj, Force_ID& id, const std::map<E_ind, BendingForce::BoundaryData::Entry>& bdata){
    if (bdata.empty()) return;
    if (id >= 0) obj.m_forces[id].target<BendingForce>()->set_boundary_data(obj, bdata);
}
void VirtualSuturerCil::set_clamped_bc(int nobj, const std::map<E_ind, BendingForce::BoundaryData::Entry>& bdata){
    set_clamped_bc(m_objs[m_linfo[nobj].m_id], m_p.m_forces.sf[nobj].bend, bdata);
}
void VirtualSuturerCil::set_pressure_force(Object3D& obj, const LeafInfo& linfo, const PresForceCtx& ctx, Force_ID& id){
    if (ctx.cP == 0.0){
        if (id >= 0){
            obj.remove_force(id);
            id = -1;
        }
    } else {
        if (id >= 0){
            obj.m_forces[id].target<SimplePressureLoad>()->setPressure(linfo.m_P * ctx.cP);
        } else {
            Force Pr = SimplePressureLoad(linfo.m_P * ctx.cP);
            id = obj.add_force(Pr);
        }
    }
}
void VirtualSuturerCil::set_pressure_force(int nobj, const PresForceCtx& ctx){
    set_pressure_force(m_objs[m_linfo[nobj].m_id], m_linfo[nobj], ctx, m_p.m_forces.sf[nobj].pres);
}
void VirtualSuturerCil::set_free_edge_force(Object3D& obj, const LeafInfo& linfo, const FreeEdgeForceCtx& ctx, Force_ID& id,
                                            Vector f, Point p0){
    if (id >= 0){
        obj.remove_force(id);
        id  = -1;
    }
    if (ctx.cFE_PS != 0.0){
        DReal S = 0;
        for (auto f: obj.m_mesh.faces()){
            auto v = vert_around(obj.m_mesh, f);
            S += sqrt((CGAL::cross_product(obj.m_x0[v[1]] - obj.m_x0[v[0]], obj.m_x0[v[2]] - obj.m_x0[v[0]]) / 2).squared_length());
        }
        std::vector<std::pair<V_ind, VertexLoad<double>::VertParam>> free_edge;
        for (auto v: obj.m_mesh.vertices()){
            if ((obj.m_boundary[v] & (4|1)) && !(obj.m_boundary[v] & 2))
                free_edge.emplace_back(v, VertexLoad<double>::VertParam{1, (obj.m_x[v] - p0)*f});
        }
        double W = ctx.cFE_PS*linfo.m_P * S / free_edge.size();
        for (auto& d: free_edge) d.second.weight *= W;
        std::function<bool(double*, double*, double&)> plane_f, plane_df;
        DReal width = ctx.FE_Ht;
        if (width < 0) width = ctx.cFE_Ht * linfo.m_Ht;
        {
            plane_f = [f, p0, width](double* x, double* F, double& d){
                double dist = (Point(x[0], x[1], x[2]) - p0) * f - d;
                double w_f = (2*dist < width) ? (1 - 2*dist / width) : 0;
                for (int k = 0; k < 3; ++k) F[k] = w_f * f[k];
                return true;
            };
            plane_df = [f, p0, width](double* x, double* J, double& d){
                std::fill(J, J+9, 0.0);
                double dist = (Point(x[0], x[1], x[2]) - p0) * f - d;
                if (2*dist < width){
                    double w_j = -2 / width;
                    for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        J[3*i + j] = w_j * f[i] * f[j];
                }
                return true;
            };
        }
        Force Pre = VertexLoad<double>().insertVertices(free_edge).setForceField(plane_f, plane_df);
        id = obj.add_force(Pre);
    }
}
void VirtualSuturerCil::set_free_edge_force(int nobj, const FreeEdgeForceCtx& ctx){
    set_free_edge_force(m_objs[m_linfo[nobj].m_id], m_linfo[nobj], ctx, m_p.m_forces.sf[nobj].free_edge, m_commissure.m_n, m_commissure.m_O);
}
void VirtualSuturerCil::set_free_membrane_forces(Object3D& obj, const FreeMembraneCxt& ctx, int ileaf){
    set_elastic_force(obj, m_linfo[ileaf], ctx.m_elast[ileaf], m_p.m_forces.sf[ileaf].elast, m_gen_force_dir, m_regenerate_elastic_force);
    set_pressure_force(obj, m_linfo[ileaf], ctx.m_pres[ileaf], m_p.m_forces.sf[ileaf].pres);
    set_free_edge_force(obj, m_linfo[ileaf], ctx.m_free_edge[ileaf], m_p.m_forces.sf[ileaf].free_edge, m_commissure.m_n, m_commissure.m_O);
}
void VirtualSuturerCil::set_free_membrane_forces(int nobj){
    set_free_membrane_forces(m_objs[m_linfo[nobj].m_id], m_p.m_fmc, nobj);
}
std::function<Vector(Point, Vector)> VirtualSuturerCil::getCilinderSuturer(Point C, Vector n){
    n /= sqrt(n.squared_length());
    return [C, n](Point x, Vector dir){
        Vector R = (C - x) - ((C - x)*n)*n;
        Vector b = CGAL::cross_product(dir, R);
        b /= sqrt(b.squared_length());
        return b;
    };
}

void VirtualSuturerCil::initial_shift_and_rotate_leaf(int ileaf){
    int i = ileaf;
    auto& o = m_objs[m_linfo[i].m_id];
    Vector f0 = (m_commissure.m_A[(i+2)%3] - m_commissure.m_A[(i+1)%3])/m_commissure.m_a[i], f1 = m_commissure.m_n;
    Vector f2 = CGAL::cross_product(f1, f0);
    align_to_vectors(&o, f0, f1);
    Point po = m_commissure.m_A[(i+2)%3] + (m_commissure.m_A[(i+1)%3] - m_commissure.m_A[(i+2)%3])/2 
                + f2 * m_p.m_shifts[i].c0 * m_commissure.m_R;
    auto init_x = o.m_mesh.add_property_map<V_ind, Point>("v:init_x").first;
    for (auto v: o.m_mesh.vertices()){
        init_x[v] = o.m_x0[v];
        o.m_x0[v] = o.m_x[v] = po + (o.m_x[v] - CGAL::ORIGIN);
    }
    cilindrize(i, f1, -f2, po);
    while (is_obj_intersect_suture_line(i)){
        for (auto v: o.m_mesh.vertices()){
            init_x[v] = o.m_x0[v];
            o.m_x0[v] = o.m_x[v] += f2 * m_p.m_shifts[i].c1 * m_commissure.m_R;
        }
    }
    o.apply_updaters();
}

void VirtualSuturerCil::initial_shift_and_rotate(){
    for (int i = 0; i < m_nleafs; ++i) initial_shift_and_rotate_leaf(i);
}
void VirtualSuturerCil::cilindrize(int nobj, Vector n, Vector r, Point new_origin){
    if (m_p.m_shifts[nobj].cKcil == 0) return;
    if (abs(m_p.m_shifts[nobj].cKcil) > 2) throw std::runtime_error("Too big initial curvature of cilindric surface");
    auto& obj = m_objs[m_linfo[nobj].m_id];
    DReal Kcil = m_p.m_shifts[nobj].cKcil / m_commissure.m_R;
    auto f2 = -r;
    auto f0 = n;
    auto f1 = CGAL::cross_product(f2, f0);
    auto    x0O = (new_origin - CGAL::ORIGIN) * f0,
            x1O = (new_origin - CGAL::ORIGIN) * f1,
            x2O = (new_origin - CGAL::ORIGIN) * f2;
    for (auto v: obj.m_mesh.vertices()) {
        auto r = obj.m_x[v] - CGAL::ORIGIN;
        auto x0 = r * f0, x1 = r * f1;
        auto dx1 = x1 - x1O;
        auto kdx1 = Kcil * dx1;
        auto kdx1p2 = kdx1 * kdx1;
        auto x0_ = x0 - x0O;
        auto x1_ = x1 + dx1*kdx1p2*(kdx1p2 - 20)/120;
        auto x2_ = x2O + dx1*kdx1/40320*(((kdx1p2-56)*kdx1p2+1680)*kdx1p2-20160);
        obj.m_x[v] = CGAL::ORIGIN + x1_ * f1 + x0_ * f0 + x2_ * f2;
        //obj.m_x0[v] = obj.m_x[v];
    }
}
bool VirtualSuturerCil::is_obj_intersect_suture_line(int nobj){
    return false;
}

void VirtualSuturerCil::align_to_vectors(Object3D* leaf, Vector right, Vector up){
    Vector f0 = right / sqrt(right.squared_length()),
            f1 = up / sqrt(up.squared_length());
    
    Point p[3];
    for (auto v: leaf->m_mesh.vertices()){
        if (leaf->m_boundary[v] == (2 | 1) )
            p[0] = leaf->m_x[v];
        else if (leaf->m_boundary[v] == (2 | 4) )
            p[1] = leaf->m_x[v];
    }
    Vector e0 = p[1] - p[0];
    e0 /= sqrt(e0.squared_length());
    //here we invert mesh orientation
    Vector e1(e0[1], -e0[0], 0);
    for (auto v: leaf->m_mesh.vertices()){
        leaf->m_x[v] = Point((leaf->m_x[v]-CGAL::ORIGIN)*e0, (leaf->m_x[v]-CGAL::ORIGIN)*e1, 0);
    }

    p[2] = *leaf->m_x.begin();
    for (auto v: leaf->m_mesh.vertices()){
        if (leaf->m_boundary[v] == (2 | 1) )
            p[0] = leaf->m_x[v];
        else if (leaf->m_boundary[v] == (2 | 4) )
            p[1] = leaf->m_x[v];
        if (p[2].y() > leaf->m_x[v].y())
            p[2] = leaf->m_x[v];
    }

    Point c = p[0] + (p[1] - p[0])/2;
    for (auto v: leaf->m_mesh.vertices()){
        leaf->m_x[v] = CGAL::ORIGIN + (leaf->m_x[v][0] - c[0])* f0 + (leaf->m_x[v][1] - c[1]) * f1;;
    }
}
std::array<VirtualSuturerCil::LeafBoundaries, 3> VirtualSuturerCil::map_cusps_to_suture_lines(){
    std::array<LeafBoundaries, 3> res;
    for (int i = 0; i < m_nleafs; ++i){
        res[i] = map_cusp_to_aortic_suture(m_objs[m_linfo[i].m_id], m_linfo[i].m_suture_line, m_linfo[i].m_smap, m_linfo[i].m_clamp_suture_b);
    }
    return res;
}
VirtualSuturerCil::LeafBoundaries VirtualSuturerCil::map_cusp_to_aortic_suture(Object3D& leaf, Curve& line, SutureMapCFD& cfd, Field& b){
    LeafBoundaries lbnd;
    lbnd.obj = &leaf;
    std::vector<V_ind>& bnd = lbnd.v_suture_bnd;
    std::vector<E_ind>& ebnd = lbnd.e_suture_bnd;
    auto bdata = &lbnd.clamp_data;
    DReal l_l = get_suture_boundary(leaf, bnd, ebnd);
    DReal l_f = get_free_boundary(leaf, lbnd.v_free_bnd, lbnd.e_free_bnd);
    lbnd.suture_shift.resize(bnd.size());
    QuickLineIntegral qli;
    qli.m_curve = [&line](double t) { auto p = line(t); return Point{p[0], p[1], p[2]}; };
    qli.setup();
    DReal l_s = qli.quick_index.back().l;
    double J0 = l_s / l_l;
    cfd.setup(J0, l_l);

    auto& sline = qli.quick_index;
    lbnd.suture_shift.front() = sline.front().p - leaf.m_x[bnd.front()];
    //leaf.m_x[bnd.front()] = sline.front().p;
    DReal t_prev = 0;
    double d1 = 0, d0 = 0, l = 0;
    auto it = sline.begin();
    for (int i = 1; i < bnd.size()-1; ++i){
        l += sqrt((leaf.m_x0[bnd[i]] - leaf.m_x0[bnd[i - 1]]).squared_length());
        assert(l_l >= l && "Wrong length");
        DReal s = cfd(l / l_l);
        assert(s >= 0 && s <= 1 && "Wrong cfd function");
        if (s >= d1){
            struct CompFunc{
                bool operator()(QuickLineIntegral::QuickIntegrData& a, DReal b) const { return a.l < b; } 
                bool operator()(DReal a, QuickLineIntegral::QuickIntegrData& b) const { return a < b.l; }
            };
            if ((it + 1)->l > s * l_s) ++it;
            else it = std::upper_bound(it+1, sline.end(), s*l_s, CompFunc());

            d1 = it->l / l_s;
            d0 = (it-1)->l / l_s;
            assert(d1 > s && d0 <= s && "Wrong dist relation");
        }
        double ww = (s - d0) / (d1 - d0);
        lbnd.suture_shift[i] = ((it-1)->p - CGAL::ORIGIN) * (1 - ww) + (it->p - CGAL::ORIGIN) * ww - (leaf.m_x0[bnd[i]] - CGAL::ORIGIN);
        //leaf.m_x[bnd[i]] = CGAL::ORIGIN + ((it-1)->p - CGAL::ORIGIN) * (1 - ww) + (it->p - CGAL::ORIGIN) * ww;
        DReal t_cur = (it-1)->t * (1 - ww) + it->t * ww;
        if (b && bdata){
            auto eb = b((t_prev + t_cur)/2);
            std::array<double, 3> ext_normal = {eb[0], eb[1], eb[2]};
            (*bdata)[ebnd[i-1]] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, ext_normal);
        }
        t_prev = t_cur;
    }
    lbnd.suture_shift.back() = sline.back().p - leaf.m_x[bnd.back()];
    //leaf.m_x[bnd.back()] = sline.back().p;
    if (b && bdata){
        auto eb = b((t_prev + 1.0)/2);
        std::array<double, 3> ext_normal = {eb[0], eb[1], eb[2]};
        (*bdata)[ebnd.back()] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::CLAMPED, ext_normal);

        for (auto e: leaf.m_mesh.edges()) {
            auto f = face_around(leaf.m_mesh, e);
            if (f[0].second && f[1].second)
                continue;
            if (bdata->find(e) == bdata->end())
                (*bdata)[e] = BendingForce::BoundaryData::Entry(BendingForce::BoundaryData::FREE);
        }
    }
    return lbnd;
}
void VirtualSuturerCil::set_dirichlet_bc_leaf(int ileaf){
    int i = ileaf;
    auto& obj = m_objs[m_linfo[i].m_id];
    static const int FIXED = 2, FREE = 1 | 4 | 8;
    DirichletCond dc = [](Object3D &obj, const V_ind &v) -> unsigned char {
        if (obj.m_boundary[v] & FIXED) return 0;
        return 7;
    };
    obj.setDirichletCond(dc);
}
void VirtualSuturerCil::set_dirichlet_bc(){
    for (int i = 0; i < m_nleafs; ++i) set_dirichlet_bc_leaf(i);
}

DReal VirtualSuturerCil::get_free_boundary(Object3D& leaf, std::vector<V_ind>& bnd, std::vector<E_ind>& ebnd){
    DReal l_l = 0;
    //compute sequential vertices of boundary on template
    std::map<V_ind, std::vector<E_ind>> bnd_map;
    E_ind bnd_start;
    static const int FIXED = 2, FREE = 1 | 4 | 8;
    for (auto e: leaf.m_mesh.edges()) {
        auto v = vert_around(leaf.m_mesh, e);
        if ((leaf.m_boundary[v[0]] & FREE) && (leaf.m_boundary[v[1]] & FREE)) {
            l_l += sqrt((leaf.m_x0[v[1]] - leaf.m_x0[v[0]]).squared_length());
            if (leaf.m_boundary[v[0]] == (2 | 1) || leaf.m_boundary[v[1]] == (2 | 1))
                bnd_start = e;
            bnd_map[v[0]].push_back(e);
            bnd_map[v[1]].push_back(e);
        }
    }
    bnd.reserve(bnd_map.size());
    ebnd.reserve(bnd_map.size() - 1);
    ebnd.push_back(bnd_start);
    auto bnd_v = vert_around(leaf.m_mesh, bnd_start);
    V_ind v = bnd_v[0];
    bnd.push_back(bnd_v[1]);
    if (bnd_map[bnd_v[0]].size() < 2) {
        v = bnd_v[1];
        bnd[0] = bnd_v[0];
    }
    bnd.push_back(v);
    while (bnd.size() < bnd_map.size()){
        auto& bb = bnd_map[v];
        if (bb.size() != 2)
            throw std::runtime_error("Template has non connected boundary or non standart boundary label: "
                                    + std::to_string(leaf.m_boundary[v]) + " v = "  + std::to_string(v));
        if (bb[0] != bnd_start)
            bnd_start = bb[0];
        else
            bnd_start = bb[1];
        bnd_v = vert_around(leaf.m_mesh, bnd_start);
        ebnd.push_back(bnd_start);
        if (bnd_v[0] != v)
            v = bnd_v[0];
        else
            v = bnd_v[1];
        bnd.push_back(v);
    }
    
    return l_l;
}

DReal VirtualSuturerCil::get_suture_boundary(Object3D& leaf, std::vector<V_ind>& bnd, std::vector<E_ind>& ebnd){
    DReal l_l = 0;
    //compute sequential vertices of boundary on template
    std::map<V_ind, std::vector<E_ind>> bnd_map;
    E_ind bnd_start;
    static const int FIXED = 2;
    for (auto e: leaf.m_mesh.edges()) {
        auto v = vert_around(leaf.m_mesh, e);
        if ((leaf.m_boundary[v[0]] & FIXED) && (leaf.m_boundary[v[1]] & FIXED)) {
            l_l += sqrt((leaf.m_x0[v[1]] - leaf.m_x0[v[0]]).squared_length());
            if (leaf.m_boundary[v[0]] == (2 | 1) || leaf.m_boundary[v[1]] == (2 | 1))
                bnd_start = e;
            bnd_map[v[0]].push_back(e);
            bnd_map[v[1]].push_back(e);
        }
    }
    bnd.reserve(bnd_map.size());
    ebnd.reserve(bnd_map.size() - 1);
    ebnd.push_back(bnd_start);
    auto bnd_v = vert_around(leaf.m_mesh, bnd_start);
    V_ind v = bnd_v[0];
    bnd.push_back(bnd_v[1]);
    if (bnd_map[bnd_v[0]].size() < 2) {
        v = bnd_v[1];
        bnd[0] = bnd_v[0];
    }
    bnd.push_back(v);
    while (bnd.size() < bnd_map.size()){
        auto& bb = bnd_map[v];
        if (bb.size() != 2)
            throw std::runtime_error("Template has non connected boundary or non standart boundary label: "
                                    + std::to_string(leaf.m_boundary[v]) + " v = "  + std::to_string(v));
        if (bb[0] != bnd_start)
            bnd_start = bb[0];
        else
            bnd_start = bb[1];
        bnd_v = vert_around(leaf.m_mesh, bnd_start);
        ebnd.push_back(bnd_start);
        if (bnd_v[0] != v)
            v = bnd_v[0];
        else
            v = bnd_v[1];
        bnd.push_back(v);
    }
    
    return l_l;
}

int virt_sut_Listener(int argc, char* argv[]){
    InterConnector ic;
    /*do something with ic*/

#ifdef USE_MAGNUM_GUI
    GuiApplication app({argc, argv});
    return app.exec();
#else
    return 0;
#endif

    return 0;
}

void VirtualSuturerCil::start_renderer(World& w, const ViewWindowContext& ctx, std::thread* &front_end){
    if (ctx.with_view) {
        front_end = new std::thread(virt_sut_Listener, ctx.argc, ctx.argv);
        auto renderer = std::make_unique<World3d::DefaultRenderer>();
        renderer->set_interconnector(&g_interconnection);
        w.setRenderer(std::move(renderer));
    }
}
void VirtualSuturerCil::stop_renderer(World& w, const ViewWindowContext& ctx, std::thread* &front_end){
    if (ctx.with_view && front_end != nullptr){
        w.clearRenderer();
        front_end->join();
        delete front_end;
        front_end = nullptr;
    }
}
// void VirtualSuturerCil::start_renderer(){
//     if (m_vw_ctx.with_view) {
//         front_end = new std::thread(virt_sut_Listener, m_vw_ctx.argc, m_vw_ctx.argv);
//         auto renderer = std::make_unique<World3d::DefaultRenderer>();
//         renderer->set_interconnector(&g_interconnection);
//         m_w.setRenderer(std::move(renderer));
//     }
// }
// void VirtualSuturerCil::stop_renderer(){
//     if (m_vw_ctx.with_view){
//         front_end->join();
//         delete front_end;
//         front_end = nullptr;
//     }
// }
void VirtualSuturerCil::SaveCxt::save_step(Object3D& obj, int nobj, int step, std::string step_name) const{
    if (save_directory.empty()) return;
    string rprefix = "";
    if (prefix.empty()){
        rprefix = "leaf_" + to_string(nobj);
    } else {
        rprefix = prefix + "_" + to_string(nobj);
    }
    switch (sign_stat){
      case NUMERIC : {
        obj.save(save_directory + "/" + rprefix + "_suture_step_" + to_string(step) + ".vtk");
        break;
    } case NAMEABLE: {
        obj.save(save_directory + "/" + rprefix+ "_suture_step_" + step_name + ".vtk");
        break;
    } default:
        break;
    }
}
