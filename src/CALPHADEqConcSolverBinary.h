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
        const double* const Lmix_B, const double* const fA,
        const double* const fB);

protected:
    virtual void RHS(const double* const x, double* const fvec);

    virtual void Jacobian(const double* const x, double** const fjac);

    double RTinv_;
    double RT_;
    double c0_;
    double hphi_;

    // energies of 2 species, in three phase each
    double fA_[3];
    double fB_[3];

    // 4 L coefficients for 3 possible phases (L, A and B)
    double Lmix_L_[4];
    double Lmix_A_[4];
};
}
#endif
