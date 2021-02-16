#ifndef included_CALPHADTieLineConcSolverTernary
#define included_CALPHADTieLineConcSolverTernary

#include "NewtonSolver.h"

namespace Thermo4PFM
{
/// solve for equilibrium compositions along a tie line
/// passing through nominal composition
class CALPHADTieLineConcSolverTernary
    : public NewtonSolver<5, CALPHADTieLineConcSolverTernary>
{
public:
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
        const double* const L_AB_L, const double* const L_AC_L,
        const double* const L_BC_L, const double* const L_AB_S,
        const double* const L_AC_S, const double* const L_BC_S,
        const double* const L_ABC_L, const double* const L_ABC_S,
        const double* const fA, const double* const fB, const double* const fC);

    /// evaluate RHS of the system of eqautions to solve for
    /// specific to this solver
    void RHS(const double* const x, double* const fvec);

    /// evaluate Jacobian of system of equations
    /// specific to this solver
    void Jacobian(const double* const x, double** const fjac);

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
    double L_AB_L_[4];
    double L_AC_L_[4];
    double L_BC_L_[4];
    double L_ABC_L_[3];

    double L_AB_S_[4];
    double L_AC_S_[4];
    double L_BC_S_[4];
    double L_ABC_S_[3];

    ///
    /// energies of 3 species, in two phase each
    ///
    double fA_[2];
    double fB_[2];
    double fC_[2];
};
}
#endif
