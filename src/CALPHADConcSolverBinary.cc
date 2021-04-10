#include "CALPHADConcSolverBinary.h"
#include "CALPHADFunctions.h"
#include "xlogx.h"

namespace Thermo4PFM
{
//=======================================================================

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
void CALPHADConcSolverBinary::computeXi(
    const double* const c, double xi[2]) const
{
    // std::cout<<"CALPHADConcSolverBinary::computeXi()"<<endl;
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
}

//=======================================================================

// solve for c=(c_L, c_A, c_B)
void CALPHADConcSolverBinary::RHS(const double* const c, double* const fvec)
{
    double xi[2] = { 0., 0. };

    computeXi(c, xi);

    fvec[0] = -c0_ + (1.0 - hphi_) * c[0] + hphi_ * c[1];
    fvec[1] = xlogx_deriv(c[0]) - xlogx_deriv(1. - c[0]) - xlogx_deriv(c[1])
              + xlogx_deriv(1. - c[1]) + (xi[0] - xi[1]);
}

//=======================================================================

void CALPHADConcSolverBinary::computeDxiDc(
    const double* const c, double dxidc[3]) const
{
    // std::cout<<"CALPHADConcSolverBinary::computeDxiDc()"<<endl;
    // loop over phases
    dxidc[0] = RTinv_
               * CALPHADcomputeFMix_deriv2Binary(
                     Lmix_L_[0], Lmix_L_[1], Lmix_L_[2], Lmix_L_[3], c[0]);

    dxidc[1] = RTinv_
               * CALPHADcomputeFMix_deriv2Binary(
                     Lmix_A_[0], Lmix_A_[1], Lmix_A_[2], Lmix_A_[3], c[1]);
}

//=======================================================================

void CALPHADConcSolverBinary::Jacobian(
    const double* const c, double** const fjac)
{
    // compute dxidc for 2 phases
    double dxidc[2];
    computeDxiDc(c, dxidc);

    fjac[0][0] = (1.0 - hphi_);
    fjac[0][1] = hphi_;

    fjac[1][0] = dxidc[0] + xlogx_deriv2(c[0]) + xlogx_deriv2(1. - c[0]);

    fjac[1][1] = -dxidc[1] - xlogx_deriv2(c[1]) - xlogx_deriv2(1. - c[1]);
}

// set values of internal variables used to evaluate
// terms in Newton iterations
void CALPHADConcSolverBinary::setup(const double c0, const double hphi,
    const double RTinv, const double* const Lmix_L, const double* const Lmix_A,
    const double* const fA, const double* const fB)
{
    c0_    = c0;
    hphi_  = hphi;
    RTinv_ = RTinv;

    for (int ii = 0; ii < 4; ii++)
        Lmix_L_[ii] = Lmix_L[ii];
    for (int ii = 0; ii < 4; ii++)
        Lmix_A_[ii] = Lmix_A[ii];
    for (int ii = 0; ii < 2; ii++)
        fA_[ii] = fA[ii];
    for (int ii = 0; ii < 2; ii++)
        fB_[ii] = fB[ii];
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
