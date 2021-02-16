#ifndef included_CALPHADEqConcSolverTernary
#define included_CALPHADEqConcSolverTernary

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class CALPHADEqConcSolverTernary
    : public NewtonSolver<4, CALPHADEqConcSolverTernary>
{
public:
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
    void setup(const double RTinv, const double* const L_AB_L,
        const double* const L_AC_L, const double* const L_BC_L,
        const double* const L_AB_S, const double* const L_AC_S,
        const double* const L_BC_S, const double* const L_ABC_L,
        const double* const L_ABC_S, const double* const fA,
        const double* const fB, const double* const fC);

    /// evaluate RHS of the system of eqautions to solve for
    /// specific to this solver
    void RHS(const double* const x, double* const fvec);

    /// evaluate Jacobian of system of equations
    /// specific to this solver
    void Jacobian(const double* const x, double** const fjac);

private:
    double RTinv_;
    double RT_;

    ///
    /// energies of 3 species, in two phase each
    ///
    double fA_[2];
    double fB_[2];
    double fC_[2];

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
};
}
#endif
