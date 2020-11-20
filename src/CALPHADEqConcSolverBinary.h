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
            d_L0[i] = std::nan("");
        for (unsigned i = 0; i < 3; i++)
            d_L1[i] = std::nan("");
        for (unsigned i = 0; i < 3; i++)
            d_L2[i] = std::nan("");
        for (unsigned i = 0; i < 3; i++)
            d_L3[i] = std::nan("");
    };

    virtual ~CALPHADEqConcentrationSolverBinary(){};

    int ComputeConcentration(double* const conc, const double RTinv,
        const double* const L0, const double* const L1, const double* const L2,
        const double* const L3, const double* const fA, const double* const fB);

protected:
    virtual void RHS(const double* const x, double* const fvec);

    virtual void Jacobian(const double* const x, double** const fjac);

    double d_RTinv;
    double d_RT;
    double d_c0;
    double d_hphi;
    double d_heta;

    // energies of 2 species, in three phase each
    double d_fA[3];
    double d_fB[3];

    // L coefficients for 3 possible phases (L, A and B)
    double d_L0[3];
    double d_L1[3];
    double d_L2[3];
    double d_L3[3];
};

#endif
