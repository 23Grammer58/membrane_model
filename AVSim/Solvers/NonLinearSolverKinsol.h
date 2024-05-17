//
// Created by alex on 10.06.2021.
//

#ifndef AORTIC_VALVE_NONLINEARSOLVERKINSOL_H
#define AORTIC_VALVE_NONLINEARSOLVERKINSOL_H

#include "LinearSolver.h"

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <kinsol/kinsol.h>
#include "NonLinearSolverInterface.h"
#include <chrono>

using SUNMatrixCustom_ID = int;
const int SUNMATRIX_NATIVE = 0;
const int SUNMATRIX_SM = 1;
struct _SUNMatrixCustom{
    SUNMatrixCustom_ID cid;
    void* content;
};
typedef struct _SUNMatrixCustom *SUNMatrixCustom;

#if SUNDIALS_VERSION_MAJOR >= 6
    #define SUN_CTX(X) SUNContext X
    #define SUN_CTX_COMMA(X) SUN_CTX(X), 
#else
    #define SUN_CTX(X) 
    #define SUN_CTX_COMMA(X)      
#endif

SUNMatrixCustom_ID SUNMatCustomID(SUNMatrix v);
SUNMatrix    SUNMatNewEmptyCustom(SUN_CTX(ctx));
void         SUNMatDestroyCustom(SUNMatrix m);

struct _SUNMatrixContent_SM {
    bool own_data;        /* data ownership flag */
    SparseMatrix* data;   /*pointer to INMOST matrix*/
};
typedef struct _SUNMatrixContent_SM *SUNMatrixContent_SM;
int          SUNMatInitOps_SM(SUNMatrix v);
SUNMatrix    SUNMatNewEmpty_SM(SUN_CTX(ctx));
SUNMatrix    SUNMat_SM(SUN_CTX_COMMA(ctx) SparseMatrix* data);
SUNMatrix    SUNMat_SM(SUN_CTX_COMMA(ctx) int N);
void         SUNMatDestroy_SM(SUNMatrix m);
SUNMatrix_ID SUNMatGetID_SM(SUNMatrix A);
SUNMatrix    SUNMatCloneEmpty_SM(SUNMatrix A);
SUNMatrix    SUNMatClone_SM(SUNMatrix A);
int          SUNMatSpace_SM(SUNMatrix A, long int *lenrw, long int *leniw);
int          SUNMatZero_SM(SUNMatrix A);
int          SUNMatCopy_SM(SUNMatrix A, SUNMatrix B);
int          SUNMatMatvec_SM(SUNMatrix A, N_Vector x, N_Vector y);
SUNMatrixContent_SM
             SUNMatContent_SM(SUNMatrix A);

using SUNLinSolCustom_ID = int;
const int SUNLINSOL_NATIVE = 0;
const int SUNLINSOL_COM = 1;
struct _SUNLinSolCustom{
    SUNLinSolCustom_ID cid;
    void* content;
};
typedef struct _SUNLinSolCustom *SUNLinSolCustom;
SUNLinSolCustom_ID  SUNLinSolCustomID(SUNLinearSolver v);
SUNLinearSolver     SUNLinSolNewEmptyCustom(SUN_CTX(ctx));
void                SUNLinSolDestroyCustom(SUNLinearSolver m);

struct _SUNLinearSolverContent_Com {
    bool own_data;                  /* data ownership flag */
    LinearSolver* data;           /*pointer to INMOST solver*/
    int updated_mat;
    int verbosity;
};
typedef struct _SUNLinearSolverContent_Com *SUNLinearSolverContent_Com;

int                  SUNLinSolInitOps_Com(SUNLinearSolver LS);
SUNLinearSolver      SUNLinSolNewEmpty_Com(SUN_CTX(ctx));
SUNLinearSolver      SUNLinSol_Com(SUN_CTX_COMMA(ctx) LinearSolver* data);
SUNLinearSolver      SUNLinSol_Com(SUN_CTX_COMMA(ctx) std::string name = "");
int                  SUNLinSolFree_Com(SUNLinearSolver S);
SUNLinearSolver_Type SUNLinSolGetType_Com(SUNLinearSolver S);
SUNLinearSolver_ID   SUNLinSolGetID_Com(SUNLinearSolver S);
int                  SUNLinSolInitialize_Com(SUNLinearSolver S);
int                  SUNLinSolSetup_Com(SUNLinearSolver S, SUNMatrix A);
int                  SUNLinSolSolve_Com(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
int                  SUNLinSolNumIters_Com(SUNLinearSolver LS);
realtype             SUNLinSolResNorm_Com(SUNLinearSolver LS);
SUNLinearSolverContent_Com
                     SUNLinSolContent_Com(SUNLinearSolver LS);
#undef SUN_CTX_COMMA
#undef SUN_CTX

class NonLinearSolverKinsol: public NonLinearSolverBase{
public:
    using Solver = LinearSolver*;
    using NVContent = N_VectorContent_Serial;
    using SMContent = SUNMatrixContent_SM;
    using LSContent = SUNLinearSolverContent_Com;
    using KinSol = void*;
    using ErrorHandler = std::function<void(int error_code, const char *module, const char *function, char *msg)>;
    using InfoHandler = std::function<void(const char *module, const char *function, char *msg)>;
    enum SolveStrategy{
        NONE = KIN_NONE,
        LINESEARCH = KIN_LINESEARCH,
        PICARD = KIN_PICARD
    };
    enum EtaChoice{
        CHOISE1 = KIN_ETACHOICE1,
        CHOISE2 = KIN_ETACHOICE2,
        CONSTANT = KIN_ETACONSTANT
    };
    enum ParamType{
        REAL,
        INTEGER,
        BOOLEAN,
        ETACHOISEENUM
    };
    struct ParamDiscr{
        std::string name;
        std::string discr;
        std::string default_val;
        ParamType type;
    };

    N_Vector x = nullptr;
    SUNMatrix A = nullptr;
    SUNLinearSolver LS = nullptr;
#if SUNDIALS_VERSION_MAJOR >= 6
    sundials::Context m_ctx;    
#endif 
    KinSol kin = nullptr;

#if SUNDIALS_VERSION_MAJOR >= 6
    NonLinearSolverKinsol(sundials::Context&& ctx);
#endif
    NonLinearSolverKinsol();
    NonLinearSolverKinsol(NonLinearSolverKinsol&&) noexcept = default;
    NonLinearSolverKinsol(const NLProblem& p, Solver s = nullptr, Matrix m = nullptr, Vector x = nullptr);
    
    NonLinearSolverKinsol& SetLinearSolver(Solver s);
    //optional setting of external matrix storage (not in owning)
    NonLinearSolverKinsol& SetMatrix(Matrix s);
    NonLinearSolverKinsol& SetResVector(Vector v, long vsize = -1);
    NonLinearSolverKinsol& SetUScaleVector(Vector u_scl, long usize = -1);
    NonLinearSolverKinsol& SetFScaleVector(Vector f_scl, long fsize = -1);
    void Clear();
    ~NonLinearSolverKinsol();
    /* actually this function is not needed for user */
    NonLinearSolverKinsol& Init();
    bool Solve() override;

    NonLinearSolverKinsol& SetNLProblem(const NLProblem& p, Solver s = nullptr, Matrix m = nullptr, Vector x = nullptr);
    NonLinearSolverKinsol& SetInitialGuess(ConstVector v, long vsize = -1);
    NonLinearSolverKinsol& SetInitialGuess(const std::vector<double>& v) { return SetInitialGuess(v.data(), v.size()); }
    NonLinearSolverKinsol& SetErrHandlerFn(ErrorHandler eh);
    NonLinearSolverKinsol& SetInfoHandlerFn(InfoHandler ih);
    NonLinearSolverKinsol& SetErrHandlerFile(std::ostream& out);
    NonLinearSolverKinsol& SetInfoHandlerFile(std::ostream& out);
    NonLinearSolverKinsol& SetVerbosityLevel(int lvl);
    /*If constraints[i] is
	     0.0 then no constraint is imposed on sol[i]
	     1.0 then sol[i] will be constrained to sol[i] ≥ 0.0
	    −1.0 then sol[i] will be constrained to sol[i] ≤ 0.0
	     2.0 then sol[i] will be constrained to sol[i] > 0.0
	    −2.0 then sol[i] will be constrained to sol[i] < 0.0
     */
    NonLinearSolverKinsol& SetConstraints(Vector& constraints, long csize = -1);

    static std::vector<ParamDiscr> GetAvaliableSetParameters();
    bool SetParameterReal(std::string param, double val) override;
    bool SetParameterInt(std::string param, int val) override;
    bool SetParameter(std::string param, std::string val) override;
    bool SetEtaConstValue(double eta);
    bool SetEtaParams(double eta_gamma, double eta_alpha);
    bool SetResMonParams(double w_min, double w_max);
    bool SetResMonConstValue(double w_const);
    bool SetMaxNewtonStep(double step);
    bool SetFuncNormTol(double fnorm);
    bool SetScaledStepTol(double scl_step_tol);
    bool SetDampingAA(double dampAA);
    bool SetNumMaxIters(long maxits);
    bool UseInitSetup(bool use);
    bool UseResMon(bool use);
    bool SetMaxSetupCalls(long max_setup_calls);
    bool SetMaxSubSetupCalls(long max_sub_setup_calls);
    bool SetEtaForm(EtaChoice ec);
    bool SetSolveStrategy(SolveStrategy strategy);
    bool UseMinEps(bool use);
    bool SetMaxBetaFails(int max_fails);
    bool SetMAA(int maa);

    std::string GetReason() const override;
    int GetReasonFlag() const { return reason_flag; }
    long GetNumFuncEvals() const override      {   long res =  0; KINGetNumFuncEvals(kin, &res);       return res; }
    long GetNumNolinSolvIters() const override {   long res =  0; KINGetNumNonlinSolvIters(kin, &res); return res; }
    long GetNumBetaCondFails() const           {   long res =  0; KINGetNumBetaCondFails(kin, &res);   return res; }
    long GetNumBacktrackOps() const            {   long res =  0; KINGetNumBacktrackOps(kin, &res);    return res; }
    double GetResidualNorm() const override    { double res = -1; KINGetFuncNorm(kin, &res);           return res; }
    double GetStepLength() const               { double res = -1; KINGetStepLength(kin, &res);         return res; }
    long GetNumJacEvals() const override       {   long res =  0; KINGetNumJacEvals(kin, &res);        return res; }
    long GetNumLinFuncEvals() const            {   long res =  0; KINGetNumLinFuncEvals(kin, &res);    return res; }
    long GetNumLinIters() const override       {   long res =  0; KINGetNumLinIters(kin, &res);        return res; }
    long GetNumLinConvFails() const            {   long res =  0; KINGetNumLinConvFails(kin, &res);    return res; }
    long GetLastLinFlag() const                {   long res =  0; KINGetLastLinFlag(kin, &res);        return res; }
    std::string GetLinReturnFlagName() const   {   return KINGetLinReturnFlagName(GetLastLinFlag()); }
    SolveStrategy GetSolveStrategy() const     {   return _strategy; }
    SUNLinearSolver GetLinearSolver()          { return LS; }
    SUNMatrix GetMatrix()                      { return A; }
    N_Vector GetVectorX()                      { return x; }

private:
    NonLinearSolverKinsol&  setVectorInternal(Vector v, long vsize, N_Vector* to);
    static int assmRHS_interface(N_Vector x, N_Vector b, void *this_obj);
    static int assmMAT_interface(N_Vector x, N_Vector b, SUNMatrix J, void *this_obj, N_Vector tmp1, N_Vector tmp2);
    static void errorHandler(int error_code, const char *module, const char *function, char *msg, void *user_data);
    static void infoHandler(const char *module, const char *function, char *msg, void *user_data);

    InfoHandler ihand;
    ErrorHandler ehand;
    N_Vector u_scale = nullptr;
    N_Vector f_scale = nullptr;
    int reason_flag = KIN_MEM_NULL;
    double _egamma = 0.9, _ealpha = 2.0;
    double _wmin = 1e-5, _wmax = 0.9;
    SolveStrategy _strategy = LINESEARCH;
    bool reinit = true;
};

static int test_nlsKinsol(int argc, char* argv[]){
    NLProblem nlp;
    nlp.m_dofs = 2;
    nlp.m_rhs = [](const double *x, double *b) -> int {
        b[0] = x[0] * x[0] + x[0] * x[1] - 3;
        b[1] = x[0] * x[1] + x[1] * x[1] - 6;
        return 0;
    };
    nlp.m_jac = [](const double *x, SparseMatrix *sm) -> int {
        (*sm)[0][0] = 2 * x[0] + x[1];
        (*sm)[0][1] = x[0];
        (*sm)[1][0] = x[1];
        (*sm)[1][1] = x[0] + 2 * x[1];
        return 0;
    };
    NonLinearSolverKinsol nlsp;
    nlsp.setProblem(nlp);
    //LinearSolver ls("eigen");
    //ls.setVerbosityLevel(3);
    //nlsp.SetLinearSolver(&ls);
    nlsp.SetInfoHandlerFn([](const char *module, const char *function, char *msg){
        std::cout << "[" << module << "] " << function << "\n   " << msg << std::endl;
    });
    nlsp.SetVerbosityLevel(1);
    nlsp.SetInitialGuess(std::vector<double>{0.5, 0.5});
#if SUNDIALS_VERSION_MAJOR >= 6
    nlsp.setInterIterationMonitorFunc([&nlp, &ctx = nlsp.m_ctx](const double * x){
        auto y = N_VNew_Serial(2, ctx);
#else
    nlsp.setInterIterationMonitorFunc([&nlp](const double * x){
        auto y = N_VNew_Serial(2);        
#endif
        auto py = N_VGetArrayPointer(y);
        nlp.m_rhs(x, py);
        std::cout << "Monitor fnorm = " << sqrt(N_VDotProd(y, y)) << std::endl;
        N_VDestroy(y);
        return 0;
    });
    bool slvFlag = nlsp.Solve();
    if (true){
        std::cout << nlsp.GetReason() << std::endl;
    }
    return 0;
}


#endif //AORTIC_VALVE_NONLINEARSOLVERKINSOL_H
