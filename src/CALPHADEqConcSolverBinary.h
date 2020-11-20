#ifndef included_CALPHADEqConcSolverBinary
#define included_CALPHADEqConcSolverBinary

#include "DampedNewtonSolver.h"

#include <math.h>

class CALPHADEqConcentrationSolverBinary : public DampedNewtonSolver
{
public:
    CALPHADEqConcentrationSolverBinary()
    {
        for (unsigned i = 0; i < 3; i++)
            L0_[i] = std::nan("");
        for (unsigned i = 0; i < 3; i++)
            L1_[i] = std::nan("");
        for (unsigned i = 0; i < 3; i++)
            L2_[i] = std::nan("");
        for (unsigned i = 0; i < 3; i++)
            L3_[i] = std::nan("");
    };

    virtual ~CALPHADEqConcentrationSolverBinary(){};

    int ComputeConcentration(double* const conc, const double RTinv,
        const double* const L0, const double* const L1, const double* const L2,
        const double* const L3, const double* const fA, const double* const fB);

protected:
    virtual void RHS(const double* const x, double* const fvec);

    virtual void Jacobian(const double* const x, double** const fjac);

    double RTinv_;
    double RT_;
    double c0_;
    double hphi_;
    double heta_;

    // energies of 2 species, in three phase each
    double fA_[3];
    double fB_[3];

    // L coefficients for 3 possible phases (L, A and B)
    double L0_[3];
    double L1_[3];
    double L2_[3];
    double L3_[3];
};

#endif
