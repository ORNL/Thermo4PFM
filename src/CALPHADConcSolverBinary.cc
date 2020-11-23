#include "CALPHADConcSolverBinary.h"
#include "CALPHADFunctions.h"
#include "xlogx.h"

#include <cassert>
#include <cmath>
#include <iostream>

static const double s_smallc     = 1.0e-8;
static const double s_inv_smallc = 1. / s_smallc;

namespace Thermo4PFM
{
//=======================================================================

CALPHADConcentrationSolverBinary::CALPHADConcentrationSolverBinary(
    const bool with_third_phase)
{
    with_third_phase_ = with_third_phase;
    N_                = with_third_phase_ ? 3 : 2;
}

//=======================================================================

void CALPHADConcentrationSolverBinary::computeXi(
    const double* const c, double xi[3]) const
{
    // std::cout<<"CALPHADConcentrationSolverBinary::computeXi()"<<endl;
    // loop over phases
    {
        double omega = CALPHADcomputeFMix_derivBinary(
            Lmix_L_[0], Lmix_L_[1], Lmix_L_[2], Lmix_L_[3], c[0]);

        double eps = fA_[0] - fB_[0];

        xi[0] = RTinv_ * (eps + omega);
        // std::cout << "L2_["<<ii<<"] = " << L2_[ii] << std::endl;
    }

    {
        double omega = CALPHADcomputeFMix_derivBinary(
            Lmix_A_[0], Lmix_A_[1], Lmix_A_[2], Lmix_A_[3], c[1]);

        double eps = fA_[1] - fB_[1];

        xi[1] = RTinv_ * (eps + omega);
        // std::cout << "L2_["<<ii<<"] = " << L2_[ii] << std::endl;
    }

    if (with_third_phase_)
    {
        double omega = CALPHADcomputeFMix_derivBinary(
            Lmix_B_[0], Lmix_B_[1], Lmix_B_[2], Lmix_B_[3], c[2]);

        double eps = fA_[2] - fB_[2];

        xi[2] = RTinv_ * (eps + omega);
        // std::cout << "L2_["<<ii<<"] = " << L2_[ii] << std::endl;
    }
}
//=======================================================================

// solve for c=(c_L, c_A, c_B)
void CALPHADConcentrationSolverBinary::RHS(
    const double* const c, double* const fvec)
{
    double xi[3] = { 0., 0., 0. };

    computeXi(c, xi);

    fvec[0] = -c0_ + (1.0 - hphi_) * c[0] + hphi_ * c[1];
    fvec[1] = xlogx_deriv(c[0]) - xlogx_deriv(1. - c[0]) - xlogx_deriv(c[1])
              + xlogx_deriv(1. - c[1]) + (xi[0] - xi[1]);

    if (with_third_phase_)
    {
        fvec[2] = xlogx_deriv(c[0]) - xlogx_deriv(1. - c[0]) - xlogx_deriv(c[2])
                  + xlogx_deriv(1. - c[2]) + (xi[0] - xi[2]);
    }
}

//=======================================================================

void CALPHADConcentrationSolverBinary::computeDxiDc(
    const double* const c, double dxidc[3]) const
{
    // std::cout<<"CALPHADConcentrationSolverBinary::computeDxiDc()"<<endl;
    // loop over phases
    dxidc[0] = RTinv_
               * CALPHADcomputeFMix_deriv2Binary(
                     Lmix_L_[0], Lmix_L_[1], Lmix_L_[2], Lmix_L_[3], c[0]);

    dxidc[1] = RTinv_
               * CALPHADcomputeFMix_deriv2Binary(
                     Lmix_A_[0], Lmix_A_[1], Lmix_A_[2], Lmix_A_[3], c[1]);

    if (with_third_phase_)
    {
        dxidc[2] = RTinv_
                   * CALPHADcomputeFMix_deriv2Binary(
                         Lmix_B_[0], Lmix_B_[1], Lmix_B_[2], Lmix_B_[3], c[2]);
    }
}

//=======================================================================

void CALPHADConcentrationSolverBinary::Jacobian(
    const double* const c, double** const fjac)
{
    // compute dxidc for up to 3 phases
    double dxidc[3];
    computeDxiDc(c, dxidc);

    fjac[0][0] = (1.0 - hphi_);
    fjac[0][1] = hphi_;
    if (N_ > 2)
    {
        fjac[0][1] -= hphi_ * heta_;
        fjac[0][2] = hphi_ * heta_;
    }

    fjac[1][0] = dxidc[0] + xlogx_deriv2(c[0]) + xlogx_deriv2(1. - c[0]);

    fjac[1][1] = -dxidc[1] - xlogx_deriv2(c[1]) - xlogx_deriv2(1. - c[1]);
}

/*
 ********************************************************************
 * conc: initial guess and final solution (concentration in each phase)
 * c0: local composition
 ********************************************************************
 */
int CALPHADConcentrationSolverBinary::ComputeConcentration(double* const conc,
    const double c0, const double hphi, const double heta, const double RTinv,
    const double* const Lmix_L, const double* const Lmix_A,
    const double* const Lmix_B, const double* const fA, const double* const fB)
{
    // std::cout<<"CALPHADConcentrationSolverBinary::ComputeConcentration()"<<endl;
    c0_    = c0;
    hphi_  = hphi;
    heta_  = heta;
    RTinv_ = RTinv;

    for (int ii = 0; ii < 4; ii++)
        Lmix_L_[ii] = Lmix_L[ii];
    for (int ii = 0; ii < 4; ii++)
        Lmix_A_[ii] = Lmix_A[ii];
    for (int ii = 0; ii < 4; ii++)
        Lmix_B_[ii] = Lmix_B[ii];
    for (int ii = 0; ii < N_; ii++)
        fA_[ii] = fA[ii];
    for (int ii = 0; ii < N_; ii++)
        fB_[ii] = fB[ii];

    return NewtonSolver::ComputeSolution(conc, N_);
}
}
