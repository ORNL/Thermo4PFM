#include "ParabolicEqConcSolverBinary.h"
#include "Determinant.h"

#include <iostream>

#define EQSCALING 1.e-4

namespace Thermo4PFM
{
//=======================================================================

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

//=======================================================================

// solve for c=(c_L, c_A)
void ParabolicEqConcSolverBinary::RHS(
    const double* const conc, // composition of species A in various phases
    double* const fvec)
{
    const double cA = conc[0];
    const double cB = conc[1];
    const double t  = temperature_;

    double fA = 0.5 * (aA_[1] * t + aA_[0]) * cA * cA
                + (bA_[1] * t + bA_[0]) * cA + (cA_[1] * t + cA_[0]);
    double fB = 0.5 * (aB_[1] * t + aB_[0]) * cB * cB
                + (bB_[1] * t + bB_[0]) * cB + (cB_[1] * t + cB_[0]);

    const double muA = (aA_[1] * t + aA_[0]) * cA + (bA_[1] * t + bA_[0]);
    const double muB = (aB_[1] * t + aB_[0]) * cB + (bB_[1] * t + bB_[0]);

    // common tangent
    fvec[0] = fA + (cB - cA) * muA - fB;

    // equal chemical potentials
    fvec[1] = EQSCALING * (muA - muB);
}

//=======================================================================

void ParabolicEqConcSolverBinary::Jacobian(
    const double* const conc, JacobianDataType** const fjac)
{
    const double cA = conc[0];
    const double cB = conc[1];
    const double t  = temperature_;

    const double muA = (aA_[1] * t + aA_[0]) * cA + (bA_[1] * t + bA_[0]);

    const double dmuAdcA = aA_[1] * t + aA_[0];

    fjac[0][0] = (aA_[1] * t + aA_[0]) * cA + (bA_[1] * t + bA_[0]) - muA
                 + (cB - cA) * dmuAdcA;
    fjac[0][1] = muA - ((aB_[1] * t + aB_[0]) * cB + bB_[1] * t + bB_[0]);

    fjac[1][0] = EQSCALING * dmuAdcA;
    fjac[1][1] = EQSCALING * (-aB_[1] * t - aB_[0]);
}

// set values of internal variables
void ParabolicEqConcSolverBinary::setup(const double temperature,
    const double coeffA[][2], const double coeffB[][2])
{
    temperature_ = temperature;

    aA_[0] = coeffA[0][0];
    aA_[1] = coeffA[0][1];
    bA_[0] = coeffA[1][0];
    bA_[1] = coeffA[1][1];
    cA_[0] = coeffA[2][0];
    cA_[1] = coeffA[2][1];

    aB_[0] = coeffB[0][0];
    aB_[1] = coeffB[0][1];
    bB_[0] = coeffB[1][0];
    bB_[1] = coeffB[1][1];
    cB_[0] = coeffB[2][0];
    cB_[1] = coeffB[2][1];
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
