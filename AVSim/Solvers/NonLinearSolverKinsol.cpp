//
// Created by alex on 10.06.2021.
//

#include "NonLinearSolverKinsol.h"

#if SUNDIALS_VERSION_MAJOR >= 6
    #define SUN_CTX(X) X
    #define SUN_CTX_COMMA(X) SUN_CTX(X), 
    #define SUN_TYPED_CTX(X) SUNContext X
    #define SUN_TYPED_CTX_COMMA(X) SUN_TYPED_CTX(X), 
#else
    #define SUN_CTX(X)
    #define SUN_CTX_COMMA(X)
    #define SUN_TYPED_CTX(X)
    #define SUN_TYPED_CTX_COMMA(X)       
#endif

SUNMatrix SUNMatNewEmptyCustom(SUN_TYPED_CTX(ctx)){
    SUNMatrix m;
    SUNMatrixCustom custom;

    m = NULL;
    m = SUNMatNewEmpty(SUN_CTX(ctx));
    if (m == NULL) return(NULL);

    custom = NULL;
    custom = (SUNMatrixCustom) malloc(sizeof *custom);
    if (custom == NULL) { SUNMatFreeEmpty(m); return(NULL); }
    m->content = custom;
    custom->cid = -1;
    custom->content = NULL;

    return m;
}

void         SUNMatDestroyCustom(SUNMatrix m){
    if (m == NULL) return;
    if (m->content != NULL){
        auto custom = static_cast<SUNMatrixCustom>(m->content);
        custom->cid = -1;
        free(m->content); m->content = NULL;
    }
    SUNMatFreeEmpty(m);
}

SUNMatrixCustom_ID SUNMatCustomID(SUNMatrix v){
    if (v == NULL || v->content == NULL) return -1;
    auto custom = static_cast<SUNMatrixCustom>(v->content);
    return custom->cid;
}

int  SUNMatInitOps_SM(SUNMatrix A){
    if (A == NULL) return SUNMAT_MEM_FAIL;

    /* Attach operations */
    A->ops->getid     = SUNMatGetID_SM;
    A->ops->clone     = SUNMatClone_SM;
    A->ops->destroy   = SUNMatDestroy_SM;
    A->ops->zero      = SUNMatZero_SM;
    A->ops->copy      = SUNMatCopy_SM;
    A->ops->scaleadd  = NULL;
    A->ops->scaleaddi = NULL;
    A->ops->matvec    = SUNMatMatvec_SM;
    A->ops->space     = SUNMatSpace_SM;

    return SUNMAT_SUCCESS;
}

SUNMatrix SUNMatNewEmpty_SM(SUN_TYPED_CTX(ctx)) {
    SUNMatrix m;
    SUNMatrixCustom custom;
    SUNMatrixContent_SM content;

    m = NULL;
    m = SUNMatNewEmptyCustom(SUN_CTX(ctx));
    if (m == NULL) return(NULL);

    SUNMatInitOps_SM(m);

    /* Create content */
    content = NULL;
    content = (SUNMatrixContent_SM) malloc(sizeof *content);
    if (content == NULL) { SUNMatFreeEmpty(m); return(NULL); }
    /* Attach content */
    custom = static_cast<SUNMatrixCustom>(m->content);
    custom->content = content;
    custom->cid = SUNMATRIX_SM;
    content->own_data = false;
    content->data = NULL;

    return m;
}

SUNMatrix SUNMat_SM(SUN_TYPED_CTX_COMMA(ctx) SparseMatrix* data) {
    SUNMatrix v = SUNMatNewEmpty_SM(SUN_CTX(ctx));
    auto content = SUNMatContent_SM(v);

    content->data = data;
    content->own_data = false;

    return v;
}

SUNMatrix SUNMat_SM(SUN_TYPED_CTX_COMMA(ctx) int N) {
    SUNMatrix v = SUNMatNewEmpty_SM(SUN_CTX(ctx));
    auto content = SUNMatContent_SM(v);

    content->data = new SparseMatrix(N);
    content->own_data = true;

    return v;
}

SUNMatrix_ID SUNMatGetID_SM(SUNMatrix A) { (void) A; return SUNMATRIX_CUSTOM; }

void SUNMatDestroy_SM(SUNMatrix m) {
    if (m == NULL) return;
    auto content = SUNMatContent_SM(m);
    if (content != NULL){
        if (content->data != NULL){
            if (content->own_data){
                delete content->data;
                content->own_data = false;
            }
            content->data = NULL;
        } else content->own_data = false;
        free(content);
        static_cast<SUNMatrixCustom>(m->content)->content = NULL;
    }
    SUNMatDestroyCustom(m);
}

SUNMatrix SUNMatCloneEmpty_SM(SUNMatrix A) {
    SUNMatrix m;
    SUNMatrixCustom custom;
    SUNMatrixContent_SM content;

    if (A == NULL) return(NULL);

    /* Create matrix */
    m = NULL;
    m = SUNMatNewEmptyCustom(SUN_CTX(A->sunctx));
    if (m == NULL) return(NULL);
    custom = static_cast<SUNMatrixCustom>(m->content);

    /* Attach operations */
    if (SUNMatCopyOps(A, m)) { SUNMatDestroyCustom(m); return(NULL); }

    /* Create content */
    content = NULL;
    content = (SUNMatrixContent_SM) malloc(sizeof *content);
    if (content == NULL) { SUNMatDestroyCustom(m); return(NULL); }

    /* Attach content */
    custom->content = content;

    /* Initialize content */
    content->own_data      = false;
    content->data          = NULL;

    return(m);
}

SUNMatrix SUNMatClone_SM(SUNMatrix A) {
    SUNMatrix v = SUNMatCloneEmpty_SM(A); if (v == NULL) return(NULL);
    auto content = SUNMatContent_SM(A);
    SparseMatrix* wdat = content ? content->data : NULL;
    if (wdat){
        content->data = new SparseMatrix(*wdat);
        content->own_data = true;
    }
    return v;
}

int SUNMatSpace_SM(SUNMatrix A, long *lenrw, long *leniw) {
    auto content = SUNMatContent_SM(A);
    if (content == NULL || content->data == NULL) {
        *lenrw = 0;
        *leniw = (A ? (A->content ? 1 : 0) : 0) + (content ? 1 : 0);
        return SUNMAT_SUCCESS;
    }
    *lenrw = 0;
    *leniw = 2;
    for (int i = 0; i < content->data->rows(); ++i)
        (*lenrw) += static_cast<long>(content->data->operator[](i).size());
    *leniw += content->data->rows() + 2*(*lenrw);

    return SUNMAT_SUCCESS;
}

int SUNMatZero_SM(SUNMatrix A) {
    auto content = SUNMatContent_SM(A);
    if (content && content->data) content->data->clear();

    return SUNMAT_SUCCESS;
}

int SUNMatCopy_SM(SUNMatrix A, SUNMatrix B) {
    auto cA = SUNMatContent_SM(A), cB = SUNMatContent_SM(B);
    if (cA && cB && cA->data && cB->data){
        *(cB->data) = *(cA->data);
        return SUNMAT_SUCCESS;
    }
    return SUNMAT_OPERATION_FAIL;
}

int SUNMatMatvec_SM(SUNMatrix A, N_Vector x, N_Vector y) {
    SparseMatrix* sm;
    auto cA = SUNMatContent_SM(A);
    sm = cA ? cA->data : NULL;
    auto cx = NV_DATA_S(x), cy = NV_DATA_S(y);
    if (!sm || !cx || !cy) return SUNMAT_MEM_FAIL;
    for (int i = 0; i < sm->rows(); ++i){
        cy[i] = 0.0;
        for (auto it: sm->operator[](i))
            cy[i] += it.second * cx[it.first];
    }

    return SUNMAT_SUCCESS;
}

SUNMatrixContent_SM SUNMatContent_SM(SUNMatrix A) {
    if (!A) return NULL;
    auto custom = static_cast<SUNMatrixCustom>(A->content);
    return custom ? static_cast<SUNMatrixContent_SM>(custom->content) : NULL;
}

SUNLinSolCustom_ID  SUNLinSolCustomID(SUNLinearSolver v){
    if (v == NULL || v->content == NULL) return -1;
    auto custom = static_cast<SUNMatrixCustom>(v->content);
    return custom->cid;
}
SUNLinearSolver SUNLinSolNewEmptyCustom(SUN_TYPED_CTX(ctx)){
    SUNLinearSolver S;
    SUNLinSolCustom content;

    if (!(S = SUNLinSolNewEmpty(SUN_CTX(ctx)))) return NULL;
    content = (SUNLinSolCustom) malloc(sizeof *content);
    if (content == NULL) { SUNLinSolFreeEmpty(S); return(NULL); }
    /* Attach content */
    S->content = content;
    content->cid = -1;

    return S;
}

void      SUNLinSolDestroyCustom(SUNLinearSolver m){
    if (m == NULL) return;
    if (m->content != NULL){
        auto custom = static_cast<SUNLinSolCustom>(m->content);
        custom->cid = -1;
        free(m->content); m->content = NULL;
    }
    SUNLinSolFreeEmpty(m);
}

int SUNLinSolInitOps_Com(SUNLinearSolver LS) {
    if (LS == NULL) return SUNLS_MEM_NULL;

    // Attach operations
    LS->ops->gettype    = SUNLinSolGetType_Com;
    LS->ops->getid      = SUNLinSolGetID_Com;
    LS->ops->initialize = SUNLinSolInitialize_Com;
    LS->ops->setup      = SUNLinSolSetup_Com;
    LS->ops->solve      = SUNLinSolSolve_Com;
    LS->ops->numiters   = SUNLinSolNumIters_Com;
    LS->ops->resnorm    = SUNLinSolResNorm_Com;
    LS->ops->free       = SUNLinSolFree_Com;

    return SUNLS_SUCCESS;
}

SUNLinearSolver SUNLinSolNewEmpty_Com(SUN_TYPED_CTX(ctx)) {
    SUNLinearSolver S;
    SUNLinSolCustom custom;
    SUNLinearSolverContent_Com content;

    if (!(S = SUNLinSolNewEmptyCustom(SUN_CTX(ctx)))) return NULL;

    SUNLinSolInitOps_Com(S);

    content = (SUNLinearSolverContent_Com) malloc(sizeof *content);
    if (content == NULL) { SUNLinSolFreeEmpty(S); return(NULL); }
    /* Attach content */
    custom = static_cast<SUNLinSolCustom>(S->content);
    custom->content = content;
    custom->cid = SUNLINSOL_COM;
    content->own_data = false;
    content->data = NULL;
    content->updated_mat = 0;
    content->verbosity = 0;

    return S;
}

SUNLinearSolver SUNLinSol_Com(SUN_TYPED_CTX_COMMA(ctx) LinearSolver *data) {
    SUNLinearSolver S = SUNLinSolNewEmpty_Com(SUN_CTX(ctx));
    auto content = SUNLinSolContent_Com(S);

    content->data = data;
    content->own_data = false;

    return S;
}

SUNLinearSolver SUNLinSol_Com(SUN_TYPED_CTX_COMMA(ctx) std::string name){
    SUNLinearSolver S = SUNLinSolNewEmpty_Com(SUN_CTX(ctx));
    auto content = SUNLinSolContent_Com(S);

    if (name != "")
        content->data = new LinearSolver(name);
    else
        content->data = new LinearSolver();
    content->own_data = true;

    return S;
}

int SUNLinSolFree_Com(SUNLinearSolver S) {
    if (S == NULL) return SUNLS_SUCCESS;
    auto content = SUNLinSolContent_Com(S);
    if (content != NULL){
        if (content->data != NULL){
            if (content->own_data){
                delete content->data;
                content->own_data = false;
            }
            content->data = NULL;
        } else content->own_data = false;
        free(content);
        static_cast<SUNLinSolCustom>(S->content)->content = NULL;
    }
    SUNLinSolDestroyCustom(S);
    return SUNLS_SUCCESS;
}

SUNLinearSolver_Type SUNLinSolGetType_Com(SUNLinearSolver S){ return SUNLINEARSOLVER_MATRIX_ITERATIVE; }
SUNLinearSolver_ID   SUNLinSolGetID_Com(SUNLinearSolver S) { return SUNLINEARSOLVER_CUSTOM; }

int SUNLinSolInitialize_Com(SUNLinearSolver S) {
    /* ensure valid options */
    if (S == NULL) return (SUNLS_MEM_NULL);
    auto cont =  SUNLinSolContent_Com(S);
    if (cont == NULL || cont->data == NULL) return (SUNLS_MEM_NULL);

    return SUNLS_SUCCESS;
}

int SUNLinSolSetup_Com(SUNLinearSolver S, SUNMatrix A) {
    auto cS =  SUNLinSolContent_Com(S);
    auto cA = SUNMatContent_SM(A);
    if (!cS || !cA) return SUNLS_MEM_NULL;
    auto Sd = cS->data;
    auto Ad = cA->data;

    auto verb = cS->verbosity;
    if (verb >= 3 || ((verb == 2))) std::cout << "\tStart precondition" << std::endl;
    Sd->SetMatrix(*Ad, true, false);
    if (verb >= 3 || ((verb == 2))) std::cout << "\tEnd precondition: prec_time = " << Sd->PreconditionerTime() <<  std::endl;
    cS->updated_mat = 1;

    return SUNLS_SUCCESS;
}

int SUNLinSolSolve_Com(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol) {
    auto cS =  SUNLinSolContent_Com(S);
    auto cA = SUNMatContent_SM(A);
    if (!cS || !cA) return SUNLS_MEM_NULL;
    auto Sd = cS->data;
    auto Ad = cA->data;
    auto verb = cS->verbosity;
    if (cS->updated_mat != 1){
        if (verb >= 3 || ((verb == 2))) std::cout << "\tStart precondition" << std::endl;
        Sd->SetMatrix(*Ad, false, true);
        if (verb >= 3 || ((verb == 2))) std::cout << "\tEnd precondition: prec_time = " << Sd->PreconditionerTime() <<  std::endl;
        cS->updated_mat = 0;
    }
    Sd->SetParameterReal("absolute_tolerance", tol);
    auto solve_start = std::chrono::steady_clock::now();
    if (verb >= 3 || ((verb == 2))) std::cout << "\tStart solve with tol = "
                                                             << tol << " | x: " << N_VGetLength(x) << "  b: " << N_VGetLength(b)
                                                             << "  mtx: " << Ad->rows() << " x " << Ad->cols() << std::endl;
    auto flag = Sd->Solve(N_VGetLength(b), NV_DATA_S(b), N_VGetLength(x), NV_DATA_S(x));
    if (verb >= 3 || ((verb == 2))) std::cout << "\tEnd solve: solve_time = " << Sd->IterationsTime() <<  std::endl;
    double t_solve =  std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start).count();
    if (!flag){
        realtype init_resid = std::sqrt(N_VDotProd(b, b));
        if (std::isnan(init_resid)) return SUNLS_PACKAGE_FAIL_UNREC;
        auto resid = Sd->ResidualNorm();
        auto its = Sd->NumIterations();
        if (verb > 0) {
            std::cout << "Solution status: solution_failed Iterations " << its << " Residual " << resid
                      << " (init_resid = " << init_resid << "). Time = " << t_solve << std::endl;
            std::cout << "Reason: " << Sd->GetReason() << std::endl;
        }
        std::string max_its_str;
        bool has_prm = Sd->GetParameter("maximum_iterations", max_its_str);
        bool max_its = (its >= (has_prm ? stoi(max_its_str) : -1));
        bool reduce_resid = (init_resid >= resid);
        bool achieve_tol = (resid <= tol);
        if (reduce_resid && achieve_tol) return SUNLS_PACKAGE_FAIL_UNREC;
        if (reduce_resid && !achieve_tol) return SUNLS_RES_REDUCED;
        if (!reduce_resid && achieve_tol) return SUNLS_PACKAGE_FAIL_UNREC;
        if (!reduce_resid && !achieve_tol) return SUNLS_CONV_FAIL;
    } else {
        auto resid = Sd->ResidualNorm();
        if (verb >= 1)
            std::cout<<"solved_succesful iterations "<<Sd->NumIterations()<<" Residual "<<Sd->ResidualNorm()<<
                     ". Time = " << t_solve << std::endl;
        if (resid <= tol) return SUNLS_SUCCESS;
        else {
            realtype init_resid = std::sqrt(N_VDotProd(b, b));
            if (std::isnan(init_resid)) return SUNLS_PACKAGE_FAIL_UNREC;
            if (init_resid >= resid) return SUNLS_CONV_FAIL;
            else return SUNLS_RES_REDUCED;
        }
    }

    return SUNLS_MEM_FAIL;
}

int SUNLinSolNumIters_Com(SUNLinearSolver LS) {
    auto content = SUNLinSolContent_Com(LS);
    if (!content || !content->data) return -1;
    return static_cast<int>(content->data->NumIterations());
}
realtype SUNLinSolResNorm_Com(SUNLinearSolver LS) {
    auto content = SUNLinSolContent_Com(LS);
    if (!content || !content->data) return -1;
    auto data = content->data;
    return data->ResidualNorm();
}

SUNLinearSolverContent_Com SUNLinSolContent_Com(SUNLinearSolver LS) {
    if (!LS) return NULL;
    auto custom = static_cast<SUNLinSolCustom>(LS->content);
    return custom ? static_cast<SUNLinearSolverContent_Com>(custom->content) : NULL;
}

#if SUNDIALS_VERSION_MAJOR >= 6
#define SUNNLS_CTX m_ctx 
#define SUNNLS_CTX_COMMA m_ctx,
#define COMMA_SUNNLS_CTX ,m_ctx
#else  
#define SUNNLS_CTX  
#define SUNNLS_CTX_COMMA 
#define COMMA_SUNNLS_CTX     
#endif

NonLinearSolverKinsol::NonLinearSolverKinsol() { 
#if SUNDIALS_VERSION_MAJOR >= 6
    m_ctx = sundials::Context(nullptr);    
#endif
    kin = KINCreate(SUNNLS_CTX);
    if (!kin) throw std::runtime_error("Error in initialization of KINCreate"); 
}

#if SUNDIALS_VERSION_MAJOR >= 6
NonLinearSolverKinsol::NonLinearSolverKinsol(sundials::Context&& ctx): m_ctx(std::move(ctx)) {
        kin = KINCreate(m_ctx);
        if (!kin) throw std::runtime_error("Error in initialization of KINCreate");
    }
#endif 

bool NonLinearSolverKinsol::SetParameterReal(std::string param, double val) {
    if (param == "EtaConstValue") return SetEtaConstValue(val);
    if (param == "EtaGamma") return SetEtaParams(val, _ealpha);
    if (param == "EtaAlpha") return SetEtaParams( _egamma, val);
    if (param == "ResMonWmin") return SetResMonParams(val, _wmax);
    if (param == "ResMonWmax") return SetResMonParams(_wmin, val);
    if (param == "ResMonConst") return SetResMonConstValue(val);
    if (param == "MaxNewtonStep") return SetMaxNewtonStep(val);
    if (param == "FuncNormTol") return SetFuncNormTol(val);
    if (param == "ScaledSteptol") return SetScaledStepTol(val);
    if (param == "DampingAA") return SetDampingAA(val);
    return false;
}

bool NonLinearSolverKinsol::SetParameterInt(std::string param, int val) {
    if (param == "MaxIters") return SetNumMaxIters(val);
    if (param == "UseInitSetup") return UseInitSetup(val>0);
    if (param == "UseResMon") return UseResMon(val>0);
    if (param == "MaxSetupCalls") return SetMaxSetupCalls(val);
    if (param == "MaxSubSetupCalls") return SetMaxSubSetupCalls(val);
    if (param == "EtaForm") {
        switch (val) {
            case KIN_ETACHOICE1: return SetEtaForm(CHOISE1);
            case KIN_ETACHOICE2: return SetEtaForm(CHOISE2);
            case KIN_ETACONSTANT: return SetEtaForm(CONSTANT);
            default: return false;
        }
    }
    if (param == "UseMinEps") return UseMinEps(val>0);
    if (param == "MaxBetaFails") return SetMaxBetaFails(val), true;
    if (param == "MAA") return SetMAA(val);
    if (param == "SolveStrategy") {
        switch (val) {
            case KIN_NONE: return SetSolveStrategy(NONE);
            case KIN_LINESEARCH: return SetSolveStrategy(LINESEARCH);
            case KIN_PICARD: return SetSolveStrategy(PICARD);
            default: return false;
        }
    }
    return false;
}

bool NonLinearSolverKinsol::SetParameter(std::string param, std::string val) {
    static std::map<std::string, ParamType> prms = {
            {"MaxIters", INTEGER},
            {"UseInitSetup", BOOLEAN},
            {"UseResMon", BOOLEAN},
            {"MaxSetupCalls", INTEGER},
            {"MaxSubSetupCalls", INTEGER},
            {"EtaForm", ETACHOISEENUM},
            {"EtaConstValue",  REAL},
            {"EtaGamma", REAL},
            {"EtaAlpha", REAL},
            {"ResMonWmin", REAL},
            {"ResMonWmax", REAL},
            {"ResMonConst", REAL},
            {"UseMinEps", BOOLEAN},
            {"MaxNewtonStep", REAL},
            {"MaxBetaFails", INTEGER},
            {"FuncNormTol", REAL},
            {"ScaledSteptol", REAL},
            {"MAA", INTEGER},
            {"DampingAA", REAL},
            {"SolveStrategy", INTEGER}
    };
    auto it = prms.find(param);
    if (it == prms.end()) return false;
    auto type = it->second;
    switch (type) {
        case REAL: return SetParameterReal(param, atof(val.c_str()));
        case INTEGER: return SetParameterInt(param, atoi(val.c_str()));
        case BOOLEAN: {
            std::string tmp = val;
            std::transform(val.begin(), val.end(), tmp.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            if (tmp == "true") return SetParameterInt(param, 1);
            if (tmp == "false") return SetParameterInt(param, 0);
            return SetParameterInt(param, atoi(val.c_str()));
        }
        case ETACHOISEENUM: return SetParameterInt(param, atoi(val.c_str()));
    }
    return false;
}

bool NonLinearSolverKinsol::SetEtaConstValue(double eta) {
    return KINSetEtaConstValue(kin, eta) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetEtaParams(double eta_gamma, double eta_alpha) {
    bool res = KINSetEtaParams(kin, eta_gamma, eta_alpha) == KIN_SUCCESS;
    if (res){
        _egamma = eta_gamma;
        _ealpha = eta_alpha;
    }
    return res;
}

bool NonLinearSolverKinsol::SetResMonParams(double w_min, double w_max) {
    bool res = KINSetResMonParams(kin, _wmin, _wmax) == KIN_SUCCESS;
    if (res){
        _wmin = w_min;
        _wmax = w_max;
    }
    return res;
}

bool NonLinearSolverKinsol::SetResMonConstValue(double w_const) {
    return KINSetResMonConstValue(kin, w_const) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetMaxNewtonStep(double step) {
    return KINSetMaxNewtonStep(kin, step) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetFuncNormTol(double fnorm) {
    return KINSetFuncNormTol(kin, fnorm) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetScaledStepTol(double scl_step_tol) {
    return KINSetScaledStepTol(kin, scl_step_tol) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetDampingAA(double dampAA) {
    if (dampAA > 1 || dampAA < 0) return false;
    return KINSetDampingAA(kin, dampAA) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetNumMaxIters(long maxits) {
    return KINSetNumMaxIters(kin, maxits) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::UseInitSetup(bool use) {
    return KINSetNoInitSetup(kin, !use) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::UseResMon(bool use) {
    return KINSetNoResMon(kin, !use) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetMaxSetupCalls(long max_setup_calls) {
    return KINSetMaxSetupCalls(kin, max_setup_calls) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetMaxSubSetupCalls(long max_sub_setup_calls) {
    return KINSetMaxSubSetupCalls(kin, max_sub_setup_calls) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetEtaForm(NonLinearSolverKinsol::EtaChoice ec) {
    return KINSetEtaForm(kin, ec) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::UseMinEps(bool use) {
    return KINSetNoMinEps(kin, !use) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetMaxBetaFails(int max_fails) {
    return KINSetMaxBetaFails(kin, max_fails) == KIN_SUCCESS;
}

bool NonLinearSolverKinsol::SetMAA(int maa) {
    bool res = KINSetMAA(kin, maa) == KIN_SUCCESS;
    if (res) reinit = true;
    return res;
}

bool NonLinearSolverKinsol::SetSolveStrategy(NonLinearSolverKinsol::SolveStrategy strategy) {
    return _strategy = strategy, true;
}

std::vector<NonLinearSolverKinsol::ParamDiscr> NonLinearSolverKinsol::GetAvaliableSetParameters() {
    static std::vector<ParamDiscr> res = {
            {"MaxIters", "Max. number of nonlinear iterations", "200", INTEGER},
            {"UseInitSetup", "Does make initial matrix setup?", "true", BOOLEAN},
            {"UseResMon", "Does perform residual monitoring?", "true", BOOLEAN},
            {"MaxSetupCalls", "Max. iterations without matrix setup", "10", INTEGER},
            // MaxSetupCalls % MaxSubSetupCalls must be equal 0!
            {"MaxSubSetupCalls", "Max. iterations between checks by the residual monitoring algorithm", "5", INTEGER},
            {"EtaForm", "Form of η coefficient", "1", ETACHOISEENUM},
            {"EtaConstValue", "Constant value of η for constant eta choice", "0.1", REAL},
            {"EtaGamma", "Values of γ for ETECHOICE2", "0.9", REAL},
            {"EtaAlpha", "Values of α for ETECHOICE2", "2.0", REAL},
            {"ResMonWmin", "Value of ω_{min}", "1e-5", REAL},
            {"ResMonWmax", "Value of ω_{max}", "0.9", REAL},
            {"ResMonConst", "Constant value of ω^∗", "0.9", REAL},
            {"UseMinEps", "Do set lower bound on epsilon?", "true", BOOLEAN},
            {"MaxNewtonStep", "Max. scaled length of Newton step", "1000||D_u*u_0||_2", REAL},
            {"MaxBetaFails", "Max. number nonlinear its with β-condition failures", "10", INTEGER},
            {"FuncNormTol", "Function-norm stopping tolerance", "(DBL_EPSILON)^(1/3)", REAL},
            {"ScaledSteptol", "Scaled-step stopping tolerance", "(DBL_EPSILON)^(2/3)", REAL},
            //MAA should be setted before KINInit and should be less, than MaxIters
            {"MAA", "Anderson Acceleration subspace size", "0", INTEGER},
            {"DampingAA", "Anderson Acceleration damping parameter beta (0 < beta ≤ 1.0)", "1", REAL},
            {"SolveStrategy", "Solving strategy: 0 for NONE, 1 for LINESEARCH and 2 for PICARD", "1", INTEGER}
    };
    return res;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetLinearSolver(NonLinearSolverKinsol::Solver s) {
    SUNLinSolFree(LS);
    LS = SUNLinSol_Com(SUNNLS_CTX_COMMA s); if (!LS ) throw std::runtime_error("Error in initialization of SUNLinSol_Com");
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::setVectorInternal(NonLinearSolverBase::Vector v, long vsize, N_Vector *to) {
    if (vsize < 0) { vsize = m_prob.m_dofs; }
    if (vsize > 0) {
        N_VDestroy(*to);
        *to = N_VMake_Serial(vsize, v COMMA_SUNNLS_CTX);
        if (!(*to) ) throw std::runtime_error("Error in initialization of N_V_Serial");
    }
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetMatrix(NonLinearSolverBase::Matrix s) {
    SUNMatDestroy(A);
    A = SUNMat_SM(SUNNLS_CTX_COMMA s); if (!A ) throw std::runtime_error("Error in initialization of SUNMat_SM");
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetResVector(NonLinearSolverBase::Vector v, long vsize) {
    return setVectorInternal(v, vsize, &x);
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetUScaleVector(NonLinearSolverBase::Vector u_scl, long usize) {
    setVectorInternal(u_scl, usize, &u_scale);
    if (!f_scale){
        f_scale = N_VClone(u_scale); if (!f_scale) throw std::runtime_error( "Error in nvclone");
        N_VConst(1.0, f_scale);
    }
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetFScaleVector(NonLinearSolverBase::Vector f_scl, long fsize) {
    return setVectorInternal(f_scl, fsize, &f_scale);
}

std::string NonLinearSolverKinsol::GetReason() const {
    switch (reason_flag) {
        case KIN_SUCCESS:
            return  "Solution success";
        case KIN_INITIAL_GUESS_OK:
            return  "The guess satisfied the system within the tolerance specified";
        case KIN_STEP_LT_STPTOL:
            return  "Solving stopped based on scaled step length. This means that the current iterate may\n"
                    "be an approximate solution of the given nonlinear system, but it is also quite possible\n"
                    "that the algorithm is “stalled” (making insufficient progress) near an invalid solution,\n"
                    "or that the scalar scsteptol is too large (see KINSetScaledStepTol in to\n"
                    "change scsteptol from its default value).";
        case KIN_MEM_NULL:
            return  "The kinsol memory block pointer was NULL.";
        case KIN_ILL_INPUT:
            return  "Initialization of SUNNonlinearSolver is incomplete.";
        case KIN_NO_MALLOC:
            return  "The kinsol memory was not allocated by a call to KINCreate.";
        case KIN_MEM_FAIL:
            return  "A memory allocation failed.";
        case KIN_LINESEARCH_NONCONV:
            return  "The line search algorithm was unable to find an iterate sufficiently distinct from the\n"
                    "current iterate, or could not find an iterate satisfying the sufficient decrease condition.\n"
                    "Failure to satisfy the sufficient decrease condition could mean the current iterate\n"
                    "is “close” to an approximate solution of the given nonlinear system, the difference\n"
                    "approximation of the matrix-vector product J(u)v is inaccurate, or the real scalar\n"
                    "scsteptol is too large.";
        case KIN_MAXITER_REACHED:
            return  "The maximum number of nonlinear iterations has been reached.";
        case KIN_MXNEWT_5X_EXCEEDED:
            return  "Five consecutive steps have been taken that satisfy the inequality ||D_u*p||_{L2} > 0.99\n"
                    "mxnewtstep, where p denotes the current step and mxnewtstep is a scalar upper\n"
                    "bound on the scaled step length. Such a failure may mean that ||D_F*F(u)||_{L2} asymp-\n"
                    "totes from above to a positive value, or the real scalar mxnewtstep is too small. By\n"
                    "default D_u = D_F = I.";
        case KIN_LINESEARCH_BCFAIL:
            return  "The line search algorithm was unable to satisfy the “beta-condition” for MXNBCF +1\n"
                    "nonlinear iterations (not necessarily consecutive), which may indicate the algorithm\n"
                    "is making poor progress.";
        case KIN_LINSOLV_NO_RECOVERY:
            return  "The linear solver's solve function failed recoverably, but the Jacobian data is already current.";
        case KIN_LINIT_FAIL:
            return  "The KINLinearSolver initialization routine (linit) encountered an error.";
        case KIN_LSETUP_FAIL:
            return  "The KINLinearSolver setup routine (lsetup) encountered an error; e.g., the user-supplied routine\n"
                    "pset (used to set up the preconditioner data) encountered an unrecoverable error.";
        case KIN_LSOLVE_FAIL:
            return  "The KINLinearSolver solve routine (lsolve) encountered an error.";
        case KIN_SYSFUNC_FAIL:
            return  "Unrecoverable error in user-supplied functions.";
        case KIN_FIRST_SYSFUNC_ERR:
            return  "User-supplied functions failed recoverably at the first call.";
        case KIN_REPTD_SYSFUNC_ERR:
            return  "User-supplied functions had repeated recoverable errors. No recovery is possible.";
        default:
            return  "Encountered unknown error";
    }
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetConstraints(NonLinearSolverBase::Vector &constraints, long csize) {
    if (csize < 0) csize = m_prob.m_dofs;
    long rsize = csize;
    if (m_prob.m_dofs > 0) rsize = m_prob.m_dofs;
    N_Vector constr;
    if (rsize == csize){
        constr = N_VMake_Serial(rsize, constraints COMMA_SUNNLS_CTX);
    } else {
        constr = N_VNew_Serial(rsize COMMA_SUNNLS_CTX);
        auto pc = N_VGetArrayPointer(constr);
        std::copy(constraints, constraints + std::min(csize, rsize), pc);
        if (csize < rsize) std::fill(pc + csize, pc + rsize, 0.0);
    }
    KINSetConstraints(kin, constr);
    N_VDestroy(constr);
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetVerbosityLevel(int lvl) {
    if (lvl <= 0) KINSetPrintLevel(kin, 0);
    else if (lvl >= 3) KINSetPrintLevel(kin, 3);
    else KINSetPrintLevel(kin, lvl);
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetInfoHandlerFile(std::ostream &out) {
    return SetInfoHandlerFn([&out](const char *module, const char *function, char *msg){
        out << "\n[" << module << "] " << function << "\n  " << msg << std::endl;
    });
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetErrHandlerFile(std::ostream &out) {
    return SetErrHandlerFn([&out](int error_code, const char *module, const char *function, char *msg){
        out << "\n"
               "[" << module << " " << ((error_code == KIN_WARNING) ? "WARNING": "ERROR") << "] " << function << "\n"
                                                                                                                 "  " << msg << "\n" << std::endl;
    });
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetInfoHandlerFn(NonLinearSolverKinsol::InfoHandler ih) {
    ihand = std::move(ih);
    KINSetInfoHandlerFn(kin, infoHandler, this);
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetErrHandlerFn(NonLinearSolverKinsol::ErrorHandler eh) {
    ehand = std::move(eh);
    KINSetErrHandlerFn(kin, errorHandler, this);
    return *this;
}

NonLinearSolverKinsol &NonLinearSolverKinsol::SetInitialGuess(NonLinearSolverBase::ConstVector v, long vsize) {
    if (vsize <= 0) { vsize = m_prob.m_dofs; }
    if (vsize <= 0) { vsize = N_VGetLength(x); }
    if (vsize <= 0) throw std::runtime_error("Can't define size of initial guess vector");
    if (!x || vsize != N_VGetLength(x)){
        N_VDestroy(x); x = N_VNew_Serial(vsize COMMA_SUNNLS_CTX);
    }
    double* px = N_VGetArrayPointer(x);
    if (px != v) std::copy(v, v+vsize, px);
    return *this;
}

bool NonLinearSolverKinsol::Solve() {
    if (reinit) Init();
    N_Vector lu_scale, lf_scale;
    lu_scale = u_scale;
    lf_scale = f_scale ? f_scale : u_scale;
    reason_flag = KINSol(kin, x, _strategy, lu_scale, lf_scale);
    return (reason_flag >= 0);
}

NonLinearSolverKinsol &NonLinearSolverKinsol::Init() {
    if (!x){
        x = N_VNewEmpty_Serial(m_prob.m_dofs COMMA_SUNNLS_CTX);
        if (!x) throw std::runtime_error("Error in N_VNewEmpty_Serial: can't create");
        N_VConst(0.0, x);
    }
    if (!A) {
        A = SUNMat_SM(SUNNLS_CTX_COMMA m_prob.m_dofs);
        if (!A) throw std::runtime_error("Error in SUNMat_SM: can't create");
    }
    if (!u_scale) {
        u_scale = N_VClone(x); if (!u_scale) throw std::runtime_error( "Error in nvclone");
        N_VConst(1.0, u_scale);
    }
    if (!LS){
        LS = SUNLinSol_Com(SUNNLS_CTX);
    }
    int ierr = 0;
#define DO(X, STR) if ((ierr = X) != KIN_SUCCESS) throw std::runtime_error(STR + std::string(": ") + std::to_string(ierr))
    DO(KINInit(kin, assmRHS_interface, x), "Error in KINInit");
    DO(KINSetUserData(kin, this), "Error in KINSetUserData");
    DO(KINSetLinearSolver(kin, LS, A), "Error in KINSetLinearSolver");
    DO(KINSetJacFn(kin, assmMAT_interface), "Error in KINSetJacFn");
#undef DO
    reinit = false;
    return *this;
}

NonLinearSolverKinsol::~NonLinearSolverKinsol() {
    Clear();
}

void NonLinearSolverKinsol::Clear() {
    if (kin) { KINFree(&kin); kin = nullptr; }
    if (LS) { SUNLinSolFree(LS); LS = nullptr; }
    if (A) { SUNMatDestroy(A); A = nullptr; }
    if (x) { N_VDestroy(x); x = nullptr; }
    if (u_scale) { N_VDestroy(u_scale); u_scale = nullptr; }
    if (f_scale) { N_VDestroy(f_scale); u_scale = nullptr; }
}

NonLinearSolverKinsol::NonLinearSolverKinsol(const NLProblem &p, NonLinearSolverKinsol::Solver s,
                                             NonLinearSolverBase::Matrix m, NonLinearSolverBase::Vector x) : NonLinearSolverKinsol()
{
    NonLinearSolverBase::setProblem(p);
    if (s) SetLinearSolver(s);
    if (m) SetMatrix(m);
    if (x) SetResVector(x, p.m_dofs);
}

NonLinearSolverKinsol& NonLinearSolverKinsol::SetNLProblem(const NLProblem& p, Solver s, Matrix m, Vector x){
    NonLinearSolverBase::setProblem(p);
    if (s) SetLinearSolver(s);
    if (m) SetMatrix(m);
    if (x) SetResVector(x, p.m_dofs);
    return *this;
}

int NonLinearSolverKinsol::assmRHS_interface(N_Vector x, N_Vector b, void *this_obj) {
    NonLinearSolverKinsol& sns = *static_cast<NonLinearSolverKinsol*>(this_obj);
    if (!sns.m_prob.m_rhs) return KIN_SYSFUNC_FAIL;
    int stat = sns.m_prob.m_rhs(N_VGetArrayPointer(x), N_VGetArrayPointer(b));
    {
        double *pb = N_VGetArrayPointer(b);
        double bnrm = 0;
        for (int i = 0; i < N_VGetLength(b); ++i) bnrm += pb[i] * pb[i];
        if (!std::isfinite(bnrm)) stat = -10;
    }
    if (stat < 0) return 1;
    return KIN_SUCCESS;
}

int NonLinearSolverKinsol::assmMAT_interface(N_Vector x, N_Vector b, SUNMatrix J, void *this_obj, N_Vector tmp1,
                                             N_Vector tmp2) {
    NonLinearSolverKinsol& sns = *static_cast<NonLinearSolverKinsol*>(this_obj);
    if (!sns.m_prob.m_jac) return KIN_SYSFUNC_FAIL;
    SUNMatContent_SM(J)->data->clear();
    int stat = sns.m_prob.m_jac(N_VGetArrayPointer(x), SUNMatContent_SM(J)->data);
    if (stat < 0) return 1;
    return KIN_SUCCESS;
}

void NonLinearSolverKinsol::errorHandler(int error_code, const char *module, const char *function, char *msg,
                                         void *user_data) {
    NonLinearSolverKinsol& sns = *static_cast<NonLinearSolverKinsol*>(user_data);
    if (sns.ehand) sns.ehand(error_code, module, function, msg);
}

void NonLinearSolverKinsol::infoHandler(const char *module, const char *function, char *msg, void *user_data) {
    NonLinearSolverKinsol& sns = *static_cast<NonLinearSolverKinsol*>(user_data);
    if (sns.ihand) sns.ihand(module, function, msg);
    if (sns.m_monitor) {
        bool iterate_func = false;
        switch (sns._strategy){
            case NONE:
            case LINESEARCH: iterate_func = !strcmp(function, "KINSol"); break;
            case PICARD: iterate_func = !strcmp(function, "KINPicardAA"); break;
        }
        if (iterate_func) sns.m_monitor(N_VGetArrayPointer(sns.x));
    }
}
