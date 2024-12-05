#include "ParabolicConcSolverBinaryThreePhase.h"
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

// solve for c=(c_L, c_A, c_B)
void ParabolicConcSolverBinaryThreePhase::RHS(
    const double* const c, double* const fvec)
{
    fvec[0] = -c0_ + hphi0_ * c[0] + hphi1_ * c[1] + hphi2_ * c[2];

    const double t   = temperature_;
    const double muL = (aL_[1] * t + aL_[0]) * c[0] + (bL_[1] * t + bL_[0]);
    const double muA = (aA_[1] * t + aA_[0]) * c[1] + (bA_[1] * t + bA_[0]);
    const double muB = (aB_[1] * t + aB_[0]) * c[2] + (bB_[1] * t + bB_[0]);
    fvec[1]          = EQSCALING * (muA - muL);
    fvec[2]          = EQSCALING * (muB - muL);
}

//=======================================================================

void ParabolicConcSolverBinaryThreePhase::Jacobian(
    const double* const c, double** const fjac)
{
    (void)c;

    const double t = temperature_;
    fjac[0][0]     = hphi0_;
    fjac[0][1]     = hphi1_;
    fjac[0][2]     = hphi2_;

    fjac[1][0] = EQSCALING * (-1. * (aL_[1] * t + aL_[0]));
    fjac[1][1] = EQSCALING * ((aA_[1] * t + aA_[0]));
    fjac[1][2] = 0.;

    fjac[2][0] = EQSCALING * (-1. * (aL_[1] * t + aL_[0]));
    fjac[2][1] = 0.;
    fjac[2][2] = EQSCALING * (aB_[1] * t + aB_[0]);
#if 0
    std::cout << "Matrix:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            std::cout << fjac[i][j] << " ";
        std::cout << std::endl;
    }
    double d = evalDeterminant<3, double>(fjac);
    std::cout << "Determinant: " << d << std::endl;
#endif
}

// set values of internal variables
void ParabolicConcSolverBinaryThreePhase::setup(const double c0,
    const double hphi0, const double hphi1, const double hphi2,
    const double temperature, const double coeffL[][2],
    const double coeffA[][2], const double coeffB[][2])
{
    c0_          = c0;
    hphi0_       = hphi0;
    hphi1_       = hphi1;
    hphi2_       = hphi2;
    temperature_ = temperature;

    aL_[0] = coeffL[0][0];
    aL_[1] = coeffL[0][1];
    bL_[0] = coeffL[1][0];
    bL_[1] = coeffL[1][1];

    aA_[0] = coeffA[0][0];
    aA_[1] = coeffA[0][1];
    bA_[0] = coeffA[1][0];
    bA_[1] = coeffA[1][1];

    aB_[0] = coeffB[0][0];
    aB_[1] = coeffB[0][1];
    bB_[0] = coeffB[1][0];
    bB_[1] = coeffB[1][1];
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
