#ifndef included_CALPHADEqConcSolverBinary
#define included_CALPHADEqConcSolverBinary

#include "DampedNewtonSolver.h"

#include <cmath>

namespace Thermo4PFM
{

class CALPHADEqConcSolverBinary : public DampedNewtonSolver
{
public:
    CALPHADEqConcSolverBinary() : DampedNewtonSolver(2)
    {
        for (unsigned i = 0; i < 4; i++)
            Lmix_L_[i] = std::nan("");
        for (unsigned i = 0; i < 4; i++)
            Lmix_A_[i] = std::nan("");
    };

    virtual ~CALPHADEqConcSolverBinary(){};

    int ComputeConcentration(double* const conc, const double RTinv,
        const double* const Lmix_L, const double* const Lmix_A,
        const double* const fA, const double* const fB);

protected:
    virtual void RHS(const double* const x, double* const fvec);

    virtual void Jacobian(const double* const x, double** const fjac);

    double RTinv_;
    double RT_;
    double c0_;
    double hphi_;

    // energies of 2 species (A and B), in two phase each
    double fA_[2];
    double fB_[2];

    // 4 L coefficients for 2 possible phases (L, A)
    double Lmix_L_[4];
    double Lmix_A_[4];
};
}
#endif
