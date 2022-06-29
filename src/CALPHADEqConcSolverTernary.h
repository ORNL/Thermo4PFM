#ifndef Thermo4PFM_included_CALPHADEqConcSolverTernary
#define Thermo4PFM_included_CALPHADEqConcSolverTernary

#include "NewtonSolver.h"
#include "datatypes.h"

namespace Thermo4PFM
{

class CALPHADEqConcSolverTernary
    : public NewtonSolver<4, CALPHADEqConcSolverTernary, JacobianDataType>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// compute equilibrium concentrations cL, cS
    /// conc: initial guess and final solution
    /// cL: conc[0], conc[1]
    /// cS: conc[2], conc[3]
    int ComputeConcentration(double* const conc, const double tol,
        const int max_iters, const double alpha = 1.)
    {
        return NewtonSolver::ComputeSolution(conc, tol, max_iters, alpha);
    }

    /// setup model paramater values to be used by solver,
    /// at a given temperature
    void setup(const double RTinv, const CalphadDataType* const L_AB_L,
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
    double RTinv_;
    double RT_;

    ///
    /// energies of 3 species, in two phase each
    ///
    CalphadDataType fA_[2];
    CalphadDataType fB_[2];
    CalphadDataType fC_[2];

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
};
}
#endif
