#ifndef Thermo4PFM_included_ParabolicConcSolverBinary
#define Thermo4PFM_included_ParabolicConcSolverBinary

#include "LinearSolver.h"

namespace Thermo4PFM
{
class ParabolicConcSolverBinary
    : public LinearSolver<2, ParabolicConcSolverBinary, double>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// compute "internal" concentrations cL, cS1, cS2 by solving KKS
    /// equations
    /// conc: solution (concentration in each phase)
    int ComputeConcentration(double* const conc)
    {
        return LinearSolver::ComputeSolution(conc);
    }

    /// setup model paramater values to be used by solver,
    /// including composition "c0" and phase fraction "hphi"
    /// to solve for
    void setup(const double c0, const double hphi0, const double hphi1,
        const double temperature, const double coeffL[][2],
        const double coeffA[][2]);

    /// evaluate RHS of the system of eqautions to solve for
    /// specific to this solver
    void RHS(const double* const x, double* const fvec);

    /// evaluate Jacobian of system of equations
    /// specific to this solver
    void Jacobian(const double* const x, double** const fjac);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

private:
    ///
    /// Nominal composition to solve for
    ///
    double c0_;

    ///
    /// phase fractions to solve for
    ///
    double hphi0_;
    double hphi1_;

    double temperature_;

    /*
     * Phase L coefficients
     * g(c,T) = (aL1_*T+aL0_)*c*c + (bL1_*T+bL0_)*c + cL1_*T+cL0_
     */
    double aL_[2];
    double bL_[2];

    /*
     * Phase A coefficients
     * g(c,T) = (aL1_*T+aL0_)*c*c + (bL1_*T+bL0_)*c + cL1_*T+cL0_
     */
    double aA_[2];
    double bA_[2];
};
}
#endif
