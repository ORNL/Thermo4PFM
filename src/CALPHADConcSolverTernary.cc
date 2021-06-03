#include "CALPHADConcSolverTernary.h"
#include "CALPHADFunctions.h"

namespace Thermo4PFM
{
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

// solve for c=(c_L, c_A)
void CALPHADConcSolverTernary::RHS(const double* const c, double* const fvec)
{
    const double* const cL = &c[0];
    const double* const cS = &c[2];

    /*
     * derivatives in liquid
     */
    double derivIdealMixL[2];
    CALPHADcomputeFIdealMix_derivTernary(RT_, cL[0], cL[1], derivIdealMixL);

    double derivFMixL[2];
    CALPHADcomputeFMix_derivTernary(
        L_AB_L_, L_AC_L_, L_BC_L_, L_ABC_L_, cL[0], cL[1], derivFMixL);

    double dfLdciL[2];
    // 1st species
    dfLdciL[0] = fA_[0] - fC_[0] + derivFMixL[0] + derivIdealMixL[0];

    // 2nd species
    dfLdciL[1] = fB_[0] - fC_[0] + derivFMixL[1] + derivIdealMixL[1];

    /*
     * derivatives in solid
     */
    double derivIdealMixS[2];
    CALPHADcomputeFIdealMix_derivTernary(RT_, cS[0], cS[1], derivIdealMixS);

    double derivFMixS[2];
    CALPHADcomputeFMix_derivTernary(
        L_AB_S_, L_AC_S_, L_BC_S_, L_ABC_S_, cS[0], cS[1], derivFMixS);

    double dfSdciS[2];
    // 1st species
    dfSdciS[0] = fA_[1] - fC_[1] + derivFMixS[0] + derivIdealMixS[0];

    // 2nd species
    dfSdciS[1] = fB_[1] - fC_[1] + derivFMixS[1] + derivIdealMixS[1];

    /*
     * system of 4 equations
     */

    // equation for 1st species
    fvec[0] = (1.0 - hphi_) * cL[0] + hphi_ * cS[0] - c0_[0];

    // equation for 2nd species
    fvec[1] = (1.0 - hphi_) * cL[1] + hphi_ * cS[1] - c0_[1];

    // equal chemical potential equation for 1st species
    fvec[2] = dfLdciL[0] - dfSdciS[0];

    // equal chemical potential equation for 2nd species
    fvec[3] = dfLdciL[1] - dfSdciS[1];
}

//=======================================================================

void CALPHADConcSolverTernary::Jacobian(
    const double* const c, JacobianDataType** const fjac)
{
    const double* const cL = &c[0];
    const double* const cS = &c[2];

    double deriv2IdealMixL[4];
    CALPHADcomputeFIdealMix_deriv2Ternary(RT_, cL[0], cL[1], deriv2IdealMixL);

    double deriv2FMixL[4];
    CALPHADcomputeFMix_deriv2Ternary(
        L_AB_L_, L_AC_L_, L_BC_L_, L_ABC_L_, cL[0], cL[1], deriv2FMixL);

    double d2fLdciL2[3]; // include only one cross term (other one equal by
                         // symmetry)
    d2fLdciL2[0] = deriv2FMixL[0] + deriv2IdealMixL[0];
    d2fLdciL2[1] = deriv2FMixL[1] + deriv2IdealMixL[1];
    d2fLdciL2[2] = deriv2FMixL[3] + deriv2IdealMixL[3];

    double deriv2IdealMixS[4];
    CALPHADcomputeFIdealMix_deriv2Ternary(RT_, cS[0], cS[1], deriv2IdealMixS);

    double deriv2FMixS[4];
    CALPHADcomputeFMix_deriv2Ternary(
        L_AB_S_, L_AC_S_, L_BC_S_, L_ABC_S_, cS[0], cS[1], deriv2FMixS);

    double d2fSdciS2[3];
    d2fSdciS2[0] = deriv2FMixS[0] + deriv2IdealMixS[0];
    d2fSdciS2[1] = deriv2FMixS[1] + deriv2IdealMixS[1];
    d2fSdciS2[2] = deriv2FMixS[3] + deriv2IdealMixS[3];

    /*
     * Jacobian:
     * f[i][j]=df[i]/dc[j]
     */
    fjac[0][0] = (1.0 - hphi_);
    fjac[0][1] = 0.;
    fjac[0][2] = hphi_;
    fjac[0][3] = 0.;

    fjac[1][0] = 0.;
    fjac[1][1] = (1.0 - hphi_);
    fjac[1][2] = 0.;
    fjac[1][3] = hphi_;

    fjac[2][0] = d2fLdciL2[0];
    fjac[2][1] = d2fLdciL2[1];
    fjac[2][2] = -d2fSdciS2[0];
    fjac[2][3] = -d2fSdciS2[1];

    fjac[3][0] = d2fLdciL2[1];
    fjac[3][1] = d2fLdciL2[2];
    fjac[3][2] = -d2fSdciS2[1];
    fjac[3][3] = -d2fSdciS2[2];
}

//=======================================================================

void CALPHADConcSolverTernary::setup(const double c0, const double c1,
    const double hphi, const double RTinv, const CalphadDataType* const L_AB_L,
    const CalphadDataType* const L_AC_L, const CalphadDataType* const L_BC_L,
    const CalphadDataType* const L_AB_S, const CalphadDataType* const L_AC_S,
    const CalphadDataType* const L_BC_S, const CalphadDataType* const L_ABC_L,
    const CalphadDataType* const L_ABC_S, const CalphadDataType* const fA,
    const CalphadDataType* const fB, const CalphadDataType* const fC)
{
    c0_[0] = c0;
    c0_[1] = c1;
    hphi_  = hphi;
    RT_    = 1. / RTinv;

    for (int ii = 0; ii < 4; ii++)
        L_AB_L_[ii] = L_AB_L[ii];
    for (int ii = 0; ii < 4; ii++)
        L_AC_L_[ii] = L_AC_L[ii];
    for (int ii = 0; ii < 4; ii++)
        L_BC_L_[ii] = L_BC_L[ii];

    for (int ii = 0; ii < 4; ii++)
        L_AB_S_[ii] = L_AB_S[ii];
    for (int ii = 0; ii < 4; ii++)
        L_AC_S_[ii] = L_AC_S[ii];
    for (int ii = 0; ii < 4; ii++)
        L_BC_S_[ii] = L_BC_S[ii];

    for (int ii = 0; ii < 3; ii++)
        L_ABC_S_[ii] = L_ABC_S[ii];
    for (int ii = 0; ii < 3; ii++)
        L_ABC_L_[ii] = L_ABC_L[ii];

    // loop over phases
    for (int ii = 0; ii < 2; ii++)
        fA_[ii] = fA[ii];
    for (int ii = 0; ii < 2; ii++)
        fB_[ii] = fB[ii];
    for (int ii = 0; ii < 2; ii++)
        fC_[ii] = fC[ii];
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
