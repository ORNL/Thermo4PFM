// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "CALPHADConcSolverBinary.h"
#include "CALPHADFunctions.h"
#include "xlogx.h"

#include <cassert>
#include <cmath>
#include <iostream>

static const double s_smallc     = 1.0e-8;
static const double s_inv_smallc = 1. / s_smallc;

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
    for (int ii = 0; ii < N_; ii++)
    {

        double omega = CALPHADcomputeFMix_derivBinary(
            L0_[ii], L1_[ii], L2_[ii], L3_[ii], c[ii]);

        double eps = fA_[ii] - fB_[ii];

        xi[ii] = RTinv_ * (eps + omega);

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

    if (N_ > 2)
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
    for (int ii = 0; ii < N_; ii++)
    {
        dxidc[ii] = RTinv_
                    * CALPHADcomputeFMix_deriv2Binary(
                          L0_[ii], L1_[ii], L2_[ii], L3_[ii], c[ii]);
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
    const double* const L0, const double* const L1, const double* const L2,
    const double* const L3, const double* const fA, const double* const fB)
{
    // std::cout<<"CALPHADConcentrationSolverBinary::ComputeConcentration()"<<endl;
    c0_    = c0;
    hphi_  = hphi;
    heta_  = heta;
    RTinv_ = RTinv;

    for (int ii = 0; ii < N_; ii++)
        L0_[ii] = L0[ii];
    for (int ii = 0; ii < N_; ii++)
        L1_[ii] = L1[ii];
    for (int ii = 0; ii < N_; ii++)
        L2_[ii] = L2[ii];
    for (int ii = 0; ii < N_; ii++)
        L3_[ii] = L3[ii];
    for (int ii = 0; ii < N_; ii++)
        fA_[ii] = fA[ii];
    for (int ii = 0; ii < N_; ii++)
        fB_[ii] = fB[ii];

    return NewtonSolver::ComputeSolution(conc, N_);
}
