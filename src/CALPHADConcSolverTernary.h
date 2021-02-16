#ifndef included_CALPHADConcSolverTernary
#define included_CALPHADConcSolverTernary

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class CALPHADConcSolverTernary
    : public NewtonSolver<4, CALPHADConcSolverTernary>
{
public:
    /// compute "internal" concentrations cL, cS by solving KKK
    /// equations
    /// conc: initial guess and final solution (concentration in each phase)
    int ComputeConcentration(double* const conc, const double tol,
        const int max_iters, const double alpha = 1.)
    {
        return NewtonSolver::ComputeSolution(conc, tol, max_iters, alpha);
    }

    /// setup model paramater values to be used by solver,
    /// including composition "c0" and phase fraction "hphi"
    /// to solve for
    void setup(const double c0, const double c1, const double hphi,
        const double RTinv, const double* const L_AB_L,
        const double* const L_AC_L, const double* const L_BC_L,
        const double* const L_AB_S, const double* const L_AC_S,
        const double* const L_BC_S, const double* const L_ABC_L,
        const double* const L_ABC_S, const double* const fA,
        const double* const fB, const double* const fC);

    /// evaluate RHS of the system of eqautions to solve for
    /// specific to this solver
    void RHS(const double* const c, double* const fvec);

    /// evaluate Jacobian of system of equations
    /// specific to this solver
    void Jacobian(const double* const c, double** const fjac);

private:
    ///
    /// Nominal composition to solve for
    ///
    double c0_[2];

    ///
    /// phase fraction to solve for
    ///
    double hphi_;

    ///
    /// energies of 3 species, in two phase each
    ///
    double fA_[2];
    double fB_[2];
    double fC_[2];

    ///
    /// L coefficients for phase L
    ///
    double L_AB_L_[4];
    double L_AC_L_[4];
    double L_BC_L_[4];
    double L_ABC_L_[3];

    ///
    /// L coefficients for phase S
    ///
    double L_AB_S_[4];
    double L_AC_S_[4];
    double L_BC_S_[4];
    double L_ABC_S_[3];

    double RT_;
};
}
#endif
