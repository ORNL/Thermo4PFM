#ifndef included_KKSdiluteBinaryConcSolver
#define included_KKSdiluteBinaryConcSolver

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class KKSdiluteBinaryConcSolver
    : public NewtonSolver<2, KKSdiluteBinaryConcSolver>
{
public:
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    /// conc: initial guess and final solution
    /// (concentration in each phase)
    int ComputeConcentration(double* const conc, const double tol,
        const int max_iters, const double alpha = 1.)
    {
        return NewtonSolver::ComputeSolution(conc, tol, max_iters, alpha);
    }

    /// setup model paramater values to be used by solver,
    /// including composition "c0" and phase fraction "hphi"
    /// to solve for
    void setup(
        const double c0, const double hphi, const double fA, const double fB);

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
    /// phase fraction to solve for
    ///
    double hphi_;

    ///
    /// model parameters
    ///
    double fA_;
    double fB_;
};
}

#endif
