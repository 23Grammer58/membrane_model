//
// Created by alex on 08.06.2021.
//

#ifndef AORTIC_VALVE_NONLINEARSOLVERPETSC_H
#define AORTIC_VALVE_NONLINEARSOLVERPETSC_H

#include "NonLinearSolverInterface.h"
#ifdef USE_INMOST
#include <petscsnes.h>

class NonLinearSolverPetsc: public NonLinearSolverBase {
public:
    static PetscErrorCode Initialize(int argc, char* argv[], const char* database = NULL, const char help[] = NULL){
        return PetscInitialize(&argc, &argv, database, help);
    }
    static PetscErrorCode Finalize() { return PetscFinalize(); }
    explicit NonLinearSolverPetsc(std::string name = ""): name{name} { }
    int Init(){
        int ierr = 0;
        ierr = SNESCreate(PETSC_COMM_WORLD,&m_snes); if (ierr != 0) return ierr;
        if (name != ""){
            ierr = SNESSetOptionsPrefix(m_snes, name.c_str());
            if (ierr != 0) return ierr;
        }

        ierr = VecCreate(PETSC_COMM_WORLD, &m_x); if (ierr != 0)  return ierr;
        if (name != ""){
            ierr = VecSetOptionsPrefix(m_x, name.c_str());
            if (ierr != 0) return ierr;
        }

        ierr = MatCreate(PETSC_COMM_WORLD,&m_J); if (ierr != 0) return ierr;
        if (name != ""){
            ierr = MatSetOptionsPrefix(m_J, name.c_str());
            if (ierr != 0) return ierr;
        }

        return ierr;
    }
    int Clear(){
        int ierr = 0;
        if (m_r) {
            ierr = VecDestroy(&m_r);
            if (ierr != 0) return ierr;
        }
        if (m_J) {
            ierr = MatDestroy(&m_J);
            if (ierr != 0) return ierr;
        }
        if (m_x) {
            ierr = VecDestroy(&m_x);
            if (ierr != 0) return ierr;
        }
        if (m_snes) ierr = SNESDestroy(&m_snes);
        return ierr;
    }
    ~NonLinearSolverPetsc(){
        if (m_r || m_J || m_x || m_snes){
            PetscBool flag = PETSC_FALSE;
            PetscInitialized(&flag);
            if (!flag){
                std::cerr << "Can't release petsc resources because petsc is finalized yet" << std::endl;
                std::abort();
            } else Clear();
        }
    }
    int setProblem(const NLProblem& prob) override{
        int ierr = 0;
        m_prob = prob;
        VecSetSizes(m_x, m_prob.m_dofs, m_prob.m_dofs);
        ierr = VecSetFromOptions(m_x); if (ierr != 0) return ierr;
        VecSet(m_x, 0);
        VecDuplicate(m_x,&m_r);
        SNESSetFunction(m_snes,m_r,_funcEval,&m_prob);

        if (prob.m_jac) {
            MatSetSizes(m_J, PETSC_DECIDE, PETSC_DECIDE, m_prob.m_dofs, m_prob.m_dofs);
            ierr = MatSetFromOptions(m_J); if (ierr != 0) return ierr;
            MatSetUp(m_J);
            SNESSetJacobian(m_snes,m_J,m_J,_jacobianFunction,&m_prob);
        }

        KSP ksp = nullptr;
        PC pc = nullptr;
        SNESGetKSP(m_snes, &ksp);
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCILU);
        KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, 2500);

        ierr = SNESSetFromOptions(m_snes); if (ierr != 0) return ierr;

        return 0;
    }
    int setInterIterationMonitorFunc(const std::function<int(ConstVector x)>& monitor) override {
        int ierr = 0;
        if (!m_monitor){
            ierr = SNESMonitorSet(m_snes, _monitor, &m_monitor, 0);
        }
        m_monitor = monitor;
        return ierr;
    }
    int setInitialGuess(const std::vector<double>& x){
        PetscErrorCode ierr = 0;
        double *px = nullptr;
        if ((ierr = VecGetArray(m_x, &px)) < 0) { return ierr;}
        std::copy(x.begin(), x.end(), px);
        ierr = VecRestoreArray(m_x, &px);
        return ierr;
    }
    bool Solve() override{
        PetscErrorCode ierr = SNESSolve(m_snes, NULL, m_x);
        return (ierr == 0);
    }
    std::string GetReason() const override{
        SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
        SNESGetConvergedReason(m_snes, &reason);
        return std::string(SNESConvergedReasons[reason]);
    }

    SNES           m_snes = nullptr;         /* nonlinear solver context */
    Vec            m_x = nullptr;            /* solution vector */
    Vec            m_r = nullptr;            /* residual vectors */
    Mat            m_J = nullptr;            /* jacobian matrix */
    std::string    name = "";
private:
    static PetscErrorCode _funcEval(SNES snes, Vec x, Vec f, void* ctx){
        auto prob = static_cast<NLProblem*>(ctx);
        PetscErrorCode ierr = 0, err = 0;
        const double *px = nullptr;
        double *pf = nullptr;
        if ((ierr = VecGetArrayRead(x, &px)) != 0) { return ierr;}
        if ((ierr = VecGetArray(f, &pf)) != 0) { VecRestoreArrayRead(x, &px); return ierr; }
        memset(pf, 0, prob->m_dofs * sizeof(double ));
        if ((ierr = prob->m_rhs(px, pf)) != 0) { VecRestoreArrayRead(x, &px); VecRestoreArray(f, &pf); return ierr; }
        err = VecRestoreArrayRead(x, &px);
        ierr = VecRestoreArray(f, &pf);
        if (ierr != 0) err = ierr;
        return err;
    }
    static PetscErrorCode _jacobianFunction(SNES snes, Vec x, Mat Amat, Mat Pmat, void *ctx){
        auto prob = static_cast<NLProblem*>(ctx);
        PetscErrorCode ierr = 0;
        const double *px = nullptr;
        if ((ierr = VecGetArrayRead(x, &px)) != 0) { return ierr; }
        SparseMatrix sm;
        sm.resize(prob->m_dofs);
        if ((ierr = prob->m_jac(px, &sm)) != 0) { return ierr; }

        int max = 0;
        {
            std::vector<int> nonzeroes(prob->m_dofs);
            for (int k = 0; k < prob->m_dofs; ++k){
                nonzeroes[k] = sm[k].size();
                if (nonzeroes[k] > max) max = nonzeroes[k];
            }
            ierr = MatXAIJSetPreallocation(Pmat, 1, nonzeroes.data(), PETSC_NULL, PETSC_NULL, PETSC_NULL);
            if (ierr < 0) return ierr;
        }
        if (max > 0) {
            std::vector<int> col_positions; col_positions.reserve(max);
            std::vector<double> col_values; col_values.reserve(max);
            for (int k = 0; k < prob->m_dofs; ++k){
                col_positions.resize(0), col_values.resize(0);
                for (auto it: sm[k]){
                    col_positions.push_back(it.first);
                    col_values.push_back(it.second);
                }
                int cols = static_cast<int>(col_positions.size());
                ierr = MatSetValues(Pmat, 1, &k, cols, col_positions.data(), col_values.data(), INSERT_VALUES);
                if (ierr < 0) return ierr;
            }
        }

        MatAssemblyBegin(Pmat,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Pmat,MAT_FINAL_ASSEMBLY);
        if (Amat != Pmat) {
            MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
        }

        ierr = VecRestoreArrayRead(x, &px);
        return ierr;
    }

    static PetscErrorCode _monitor(SNES snes,PetscInt its,PetscReal fnorm,void *ctx){
        auto monitor = static_cast<std::function<int(ConstVector x)>*>(ctx);
        Vec x = nullptr;
        SNESGetSolution(snes,&x);
        const double* px = nullptr;
        VecGetArrayRead(x, &px);
        monitor->operator()(px);
        VecRestoreArrayRead(x, &px);

        PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);
        return 0;
    }
};

static int test_nlsPetsc(int argc, char* argv[]){
    int ierr = NonLinearSolverPetsc::Initialize(argc, argv);
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
    NonLinearSolverPetsc nlsp;
    nlsp.Init();
    nlsp.setProblem(nlp);
    nlsp.setInitialGuess({0.5, 0.5});
    nlsp.setInterIterationMonitorFunc([](const double *){ return 0; });
    bool slvFlag = nlsp.Solve();
    if (true){
        std::cout << nlsp.GetReason() << std::endl;
    }
    nlsp.Clear();
    NonLinearSolverPetsc::Finalize();
    return 0;
}
#endif

#endif //AORTIC_VALVE_NONLINEARSOLVERPETSC_H
