#ifndef Thermo4PFM_included_ParabolicEqConcSolverBinary
#define Thermo4PFM_included_ParabolicEqConcSolverBinary

#include "NewtonSolver.h"
#include "datatypes.h"

namespace Thermo4PFM
{

class ParabolicEqConcSolverBinary
    : public NewtonSolver<2, ParabolicEqConcSolverBinary, JacobianDataType>
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
    void setup(const double temperature, const double coeffA[][2],
        const double coeffB[][2]);

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
    double temperature_;

    /*
     * Phase A coefficients
     * g(c,T) = (aL1_*T+aL0_)*c*c + (bL1_*T+bL0_)*c + cL1_*T+cL0_
     */
    double aA_[2];
    double bA_[2];
    double cA_[2];

    /*
     * Phase B coefficients
     * g(c,T) = (aL1_*T+aL0_)*c*c + (bL1_*T+bL0_)*c + cL1_*T+cL0_
     */
    double aB_[2];
    double bB_[2];
    double cB_[2];
};
}
#endif
