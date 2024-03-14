#include "KKSdiluteBinaryConcSolver.h"
#include "xlogx.h"

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
void KKSdiluteBinaryConcSolver::RHS(const double* const c, double* const fvec)
{
    fvec[0] = -c0_ + hphi1_ * c[0] + hphi0_ * c[1];
    fvec[1] = xlogx_deriv(c[0]) - xlogx_deriv(1. - c[0]) - xlogx_deriv(c[1])
              + xlogx_deriv(1. - c[1]) - (fA_ - fB_);
}

//=======================================================================

void KKSdiluteBinaryConcSolver::Jacobian(
    const double* const c, JacobianDataType** const fjac)
{
    fjac[0][0] = hphi1_;
    fjac[0][1] = hphi0_;

    fjac[1][0] = xlogx_deriv2(c[0]) + xlogx_deriv2(1. - c[0]);
    fjac[1][1] = -xlogx_deriv2(c[1]) - xlogx_deriv2(1. - c[1]);
}

void KKSdiluteBinaryConcSolver::setup(const double c0, const double hphi0,
    const double hphi1, const double fA, const double fB)
{
    c0_    = c0;
    hphi0_ = hphi0;
    hphi1_ = hphi1;
    fA_    = fA;
    fB_    = fB;
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
