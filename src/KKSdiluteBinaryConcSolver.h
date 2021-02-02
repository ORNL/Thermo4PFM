#ifndef included_KKSdiluteBinaryConcSolver
#define included_KKSdiluteBinaryConcSolver

#include "NewtonSolver.h"

namespace Thermo4PFM
{

class KKSdiluteBinaryConcSolver
    : public NewtonSolver<2, KKSdiluteBinaryConcSolver>
{
public:
    // conc: initial guess and final solution (concentration in each phase)
    int ComputeConcentration(double* const conc)
    {
        return NewtonSolver::ComputeSolution(conc);
    }

    void setup(
        const double c0, const double hphi, const double fA, const double fB);

    void RHS(const double* const x, double* const fvec);

    void Jacobian(const double* const x, double** const fjac);

private:
    double fA_;
    double fB_;

    double c0_;
    double hphi_;
};
}

#endif
