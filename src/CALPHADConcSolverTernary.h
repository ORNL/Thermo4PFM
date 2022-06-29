#ifndef Thermo4PFM_included_CALPHADConcSolverTernary
#define Thermo4PFM_included_CALPHADConcSolverTernary

#include "NewtonSolver.h"

#include "datatypes.h"

namespace Thermo4PFM
{

class CALPHADConcSolverTernary
    : public NewtonSolver<4, CALPHADConcSolverTernary, JacobianDataType>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
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
        const double RTinv, const CalphadDataType* const L_AB_L,
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
    void RHS(const double* const c, double* const fvec);

    /// evaluate Jacobian of system of equations
    /// specific to this solver
    void Jacobian(const double* const c, JacobianDataType** const fjac);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

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
    CalphadDataType fA_[2];
    CalphadDataType fB_[2];
    CalphadDataType fC_[2];

    ///
    /// L coefficients for phase L
    ///
    CalphadDataType L_AB_L_[4];
    CalphadDataType L_AC_L_[4];
    CalphadDataType L_BC_L_[4];
    CalphadDataType L_ABC_L_[3];

    ///
    /// L coefficients for phase S
    ///
    CalphadDataType L_AB_S_[4];
    CalphadDataType L_AC_S_[4];
    CalphadDataType L_BC_S_[4];
    CalphadDataType L_ABC_S_[3];

    double RT_;
};
}
#endif
