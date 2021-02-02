#include "KKSdiluteBinaryConcSolver.h"
#include "xlogx.h"

#include <cassert>
#include <cmath>
#include <iostream>

namespace Thermo4PFM
{

//=======================================================================

// solve for c=(c_L, c_A)
void KKSdiluteBinaryConcSolver::RHS(const double* const c, double* const fvec)
{
    fvec[0] = -c0_ + (1.0 - hphi_) * c[0] + hphi_ * c[1];
    fvec[1] = xlogx_deriv(c[0]) - xlogx_deriv(1. - c[0]) - xlogx_deriv(c[1])
              + xlogx_deriv(1. - c[1]) - (fA_ - fB_);
}

//=======================================================================

void KKSdiluteBinaryConcSolver::Jacobian(
    const double* const c, double** const fjac)
{
    fjac[0][0] = (1.0 - hphi_);
    fjac[0][1] = hphi_;

    fjac[1][0] = xlogx_deriv2(c[0]) + xlogx_deriv2(1. - c[0]);
    fjac[1][1] = -xlogx_deriv2(c[1]) - xlogx_deriv2(1. - c[1]);
}

/*
 * c0: local composition
 */
void KKSdiluteBinaryConcSolver::setup(
    const double c0, const double hphi, const double fA, const double fB)
{
    // std::cout<<"KKSdiluteBinaryConcSolver::ComputeConcentration()"<<endl;
    c0_   = c0;
    hphi_ = hphi;
    fA_   = fA;
    fB_   = fB;
}

/*
 * conc: initial guess and final solution (concentration in each phase)
 */
}
