#ifndef Thermo4PFM_included_CALPHADTieLineConcSolverTernary
#define Thermo4PFM_included_CALPHADTieLineConcSolverTernary

#include "NewtonSolver.h"

#include "datatypes.h"

namespace Thermo4PFM
{
/// solve for equilibrium compositions along a tie line
/// passing through nominal composition
class CALPHADTieLineConcSolverTernary
    : public NewtonSolver<5, CALPHADTieLineConcSolverTernary, JacobianDataType>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// input x: initial values for cL_0, cL_1, cS_0, cS_1
    /// and phase fraction
    /// output x: ceqL_0, ceqL_1, ceqS_0, ceqS_1, and phi
    int ComputeConcentration(double* const x, const double tol,
        const int max_iters, const double alpha = 1.)
    {
        return NewtonSolver::ComputeSolution(x, tol, max_iters, alpha);
    }

    /// setup model paramater values to be used by solver,
    /// at a given temperature, including nominal composition
    /// c0, c1
    void setup(const double c0, const double c1, const double RTinv,
        const CalphadDataType* const L_AB_L,
        const CalphadDataType* const L_AC_L,
        const CalphadDataType* const L_BC_L,
        const CalphadDataType* const L_AB_S,
        const CalphadDataType* const L_AC_S,
        const CalphadDataType* const L_BC_S,
        const CalphadDataType* const L_ABC_L,
        const CalphadDataType* const L_ABC_S, const CalphadDataType* const fA,
        const CalphadDataType* const fB, const CalphadDataType* const fC);

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
    ///
    /// nominal composition (defining tie line)
    ///
    double conc_[2];

    double RTinv_;
    double RT_;

    ///
    /// L coefficients for 2 possible phases (L and S)
    ///
    CalphadDataType L_AB_L_[4];
    CalphadDataType L_AC_L_[4];
    CalphadDataType L_BC_L_[4];
    CalphadDataType L_ABC_L_[3];

    CalphadDataType L_AB_S_[4];
    CalphadDataType L_AC_S_[4];
    CalphadDataType L_BC_S_[4];
    CalphadDataType L_ABC_S_[3];

    ///
    /// energies of 3 species, in two phase each
    ///
    CalphadDataType fA_[2];
    CalphadDataType fB_[2];
    CalphadDataType fC_[2];
};
}
#endif
