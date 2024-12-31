#include "ParabolicConcSolverBinary.h"

#define EQSCALING 1.e-4

namespace Thermo4PFM
{
//=======================================================================

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

//=======================================================================

// solve for c=(c_L, c_A, c_B)
void ParabolicConcSolverBinary::RHS(const double* const c, double* const fvec)
{
    const double t   = temperature_;
    const double muL = (aL_[1] * t + aL_[0]) * c[0] + (bL_[1] * t + bL_[0]);
    const double muA = (aA_[1] * t + aA_[0]) * c[1] + (bA_[1] * t + bA_[0]);

    fvec[0] = -c0_ + hphi0_ * c[0] + hphi1_ * c[1];
    fvec[1] = EQSCALING * (muA - muL);
}

//=======================================================================

void ParabolicConcSolverBinary::Jacobian(
    const double* const c, double** const fjac)
{
    (void)c;

    const double t = temperature_;
    fjac[0][0]     = hphi0_;
    fjac[0][1]     = hphi1_;

    fjac[1][0] = EQSCALING * (-1. * (aL_[1] * t + aL_[0]));
    fjac[1][1] = EQSCALING * ((aA_[1] * t + aA_[0]));
}

// set values of internal variables
void ParabolicConcSolverBinary::setup(const double c0, const double hphi0,
    const double hphi1, const double temperature, const double coeffL[][2],
    const double coeffA[][2])
{
    c0_          = c0;
    hphi0_       = hphi0;
    hphi1_       = hphi1;
    temperature_ = temperature;

    aL_[0] = coeffL[0][0];
    aL_[1] = coeffL[0][1];
    bL_[0] = coeffL[1][0];
    bL_[1] = coeffL[1][1];

    aA_[0] = coeffA[0][0];
    aA_[1] = coeffA[0][1];
    bA_[0] = coeffA[1][0];
    bA_[1] = coeffA[1][1];
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
