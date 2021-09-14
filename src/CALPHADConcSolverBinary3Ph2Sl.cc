#include "CALPHADConcSolverBinary3Ph2Sl.h"
#include "CALPHADFunctions.h"
#include "xlogx.h"
#include <iostream>

namespace Thermo4PFM
{
//=======================================================================

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
void CALPHADConcSolverBinary3Ph2Sl::computeXi(
    const double* const c, double xi[3]) const
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
            Lmix_S0_[0], Lmix_S0_[1], Lmix_S0_[2], Lmix_S0_[3], c[1]);

        double eps = fA_[1] - fB_[1];

        xi[1] = RTinv_ * (eps + omega);
        // std::cout << "L2_["<<ii<<"] = " << L2_[ii] << std::endl;
    }

    {
        double omega = CALPHADcomputeFMix_derivBinary(
            Lmix_S1_[0], Lmix_S1_[1], Lmix_S1_[2], Lmix_S1_[3], c[2]);

        double eps = fA_[2] - fB_[2];

        xi[2] = RTinv_ * (eps + omega);
        // std::cout << "L2_["<<ii<<"] = " << L2_[ii] << std::endl;
    }
}

//=======================================================================

// solve for c=(c_L, c_A, c_B)
void CALPHADConcSolverBinary3Ph2Sl::RHS(
    const double* const c, double* const fvec)
{
    // Get the sublattice site fractions, assuming (A)_p (A,B)_q
    double ypp_A[3];
    for (int i = 0; i < 3; ++i)
    {
        ypp_A[i] = (p_[i] + q_[i]) * c[i] - p_[i];
    }

    // Get xi
    double xi[3] = { 0., 0., 0. };
    computeXi(ypp_A, xi);

    fvec[0] = -c0_ + hphi0_ * c[0] + hphi1_ * c[1] + hphi2_ * c[2];

    // We can choose to enforce two of the three chemical potential equilities.
    // Which two are chosen can impact the convergence rate.

    fvec[1] = q_[0] * (xlogx_deriv(ypp_A[0]) - xlogx_deriv(1. - ypp_A[0]))
              - q_[1] * (xlogx_deriv(ypp_A[1]) - xlogx_deriv(1. - ypp_A[1]))
              + (xi[0] - xi[1]);
    fvec[2] = q_[1] * (xlogx_deriv(ypp_A[1]) - xlogx_deriv(1. - ypp_A[1]))
              - q_[2] * (xlogx_deriv(ypp_A[2]) - xlogx_deriv(1. - ypp_A[2]))
              + (xi[1] - xi[2]);

    // Debugging output
    /*
    std::cout << "compositions: ";
    std::cout << c[0] << " " << c[1] << " " << c[2] << std::endl;

    std::cout << "p: " << p_[0] << " " << p_[1] << " " << p_[2] << std::endl;

    std::cout << "xi: " << xi[0] << " " << xi[1] << " " << xi[2] << std::endl;

    std::cout << "chemical potentials: ";
    double temp = (q_[0] * xlogx_deriv(ypp_A[0])
                   - q_[0] * xlogx_deriv(1. - ypp_A[0]) + xi[0]);
    std::cout << temp << " ";
    temp = (q_[1] * xlogx_deriv(ypp_A[1]) - q_[1] * xlogx_deriv(1. - ypp_A[1])
            + xi[1]);
    std::cout << temp << " ";
    temp = (q_[2] * xlogx_deriv(ypp_A[2]) - q_[2] * xlogx_deriv(1. - ypp_A[2])
            + xi[2]);
    std::cout << temp << std::endl;

    std::cout << "composition sum: ";
    std::cout << c0_ << " ";
    temp = hphi0_ * c[0] + hphi1_ * c[1] + hphi2_ * c[2];
    std::cout << temp << std::endl;

    std::cout << "interpolation terms: ";
    std::cout << hphi0_ << " " << hphi1_ << " " << hphi2_ << std::endl;
    */
}

//=======================================================================

void CALPHADConcSolverBinary3Ph2Sl::computeDxiDc(
    const double* const c, double dxidc[3]) const
{
    // std::cout<<"CALPHADConcSolverBinary::computeDxiDc()"<<endl;
    // loop over phases
    dxidc[0] = RTinv_
               * CALPHADcomputeFMix_deriv2Binary(
                   Lmix_L_[0], Lmix_L_[1], Lmix_L_[2], Lmix_L_[3], c[0]);

    dxidc[1] = RTinv_
               * CALPHADcomputeFMix_deriv2Binary(
                   Lmix_S0_[0], Lmix_S0_[1], Lmix_S0_[2], Lmix_S0_[3], c[1]);

    dxidc[2] = RTinv_
               * CALPHADcomputeFMix_deriv2Binary(
                   Lmix_S1_[0], Lmix_S1_[1], Lmix_S1_[2], Lmix_S1_[3], c[2]);
}

//=======================================================================

void CALPHADConcSolverBinary3Ph2Sl::Jacobian(
    const double* const c, JacobianDataType** const fjac)
{
    // Get the sublattice site fractions, assuming (A)_p (A,B)_q
    double ypp_A[3];
    for (int i = 0; i < 3; ++i)
    {
        ypp_A[i] = (p_[i] + q_[i]) * c[i] - p_[i];
    }

    // compute dxidc for 2 phases
    double dxidc[3];
    computeDxiDc(ypp_A, dxidc);

    fjac[0][0] = hphi0_;
    fjac[0][1] = hphi1_;
    fjac[0][2] = hphi2_;

    fjac[1][0] = (p_[0] + q_[0])
                 * (dxidc[0] + q_[0] * xlogx_deriv2(ypp_A[0])
                     + q_[0] * xlogx_deriv2(1. - ypp_A[0]));
    fjac[1][1] = (p_[1] + q_[1])
                 * (-dxidc[1] - q_[1] * xlogx_deriv2(ypp_A[1])
                     - q_[1] * xlogx_deriv2(1. - ypp_A[1]));
    fjac[1][2] = 0.;

    fjac[2][0] = 0.;
    fjac[2][1] = (p_[1] + q_[1])
                 * (dxidc[1] + q_[1] * xlogx_deriv2(ypp_A[1])
                     - q_[1] * xlogx_deriv2(1. - ypp_A[1]));
    fjac[2][2] = (p_[2] + q_[2])
                 * (-dxidc[2] - q_[2] * xlogx_deriv2(ypp_A[2])
                     - q_[2] * xlogx_deriv2(1. - ypp_A[2]));
}

// set values of internal variables used to evaluate
// terms in Newton iterations
void CALPHADConcSolverBinary3Ph2Sl::setup(const double c0, const double hphi0,
    const double hphi1, const double hphi2, const double RTinv,
    const CalphadDataType* const Lmix_L, const CalphadDataType* const Lmix_S0,
    const CalphadDataType* const Lmix_S1, const CalphadDataType* const fA,
    const CalphadDataType* const fB, const int* const p, const int* const q)
{
    c0_    = c0;
    hphi0_ = hphi0;
    hphi1_ = hphi1;
    hphi2_ = hphi2;
    RTinv_ = RTinv;

    for (int ii = 0; ii < 4; ii++)
    {
        Lmix_L_[ii]  = Lmix_L[ii];
        Lmix_S0_[ii] = Lmix_S0[ii];
        Lmix_S1_[ii] = Lmix_S1[ii];
    }
    for (int ii = 0; ii < 3; ii++)
    {
        fA_[ii] = fA[ii];
        fB_[ii] = fB[ii];
        p_[ii]  = static_cast<double>(p[ii]);
        q_[ii]  = static_cast<double>(q[ii]);
    }
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
