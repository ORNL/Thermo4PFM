#ifndef Thermo4PFM_included_CALPHADEqConcSolverBinary
#define Thermo4PFM_included_CALPHADEqConcSolverBinary

#include "NewtonSolver.h"
#include "datatypes.h"

namespace Thermo4PFM
{

class CALPHADEqConcSolverBinary
    : public NewtonSolver<2, CALPHADEqConcSolverBinary, JacobianDataType>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// compute equilibrium concentrations cL, cS
    /// conc: initial guess and final solution
    int ComputeConcentration(double* const conc, const double tol,
        const int max_iters, const double alpha = 1.)
    {
        return NewtonSolver::ComputeSolution(conc, tol, max_iters, alpha);
    }

    /// setup model paramater values to be used by solver,
    /// at a given temperature
    void setup(const double RT, const CalphadDataType* const Lmix_L,
        const CalphadDataType* const Lmix_A, const CalphadDataType* const fA,
        const CalphadDataType* const fB);

    /// evaluate RHS of the system of eqautions to solve for
    /// specific to this solver
    void RHS(const double* const x, double* const fvec);

    /// evaluate Jacobian of system of equations
    /// specific to this solver
    void Jacobian(const double* const x, JacobianDataType** const fjac);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

private:
    double RT_;

    ///
    /// energies of 2 species (A and B), in two phase each
    ///
    CalphadDataType fA_[2];
    CalphadDataType fB_[2];

    ///
    /// 4 L coefficients for 2 possible phases (L, A)
    ///
    CalphadDataType Lmix_L_[4];
    CalphadDataType Lmix_A_[4];
};
}
#endif
