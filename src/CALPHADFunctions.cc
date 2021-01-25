#include "CALPHADFunctions.h"
#include "xlogx.h"

//#include <cassert>

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double CALPHADcomputeFMixBinary(const double l0, const double l1,
    const double l2, const double l3, const double conc)
{
    double two_c_minus_one = 2.0 * conc - 1.0;

    double fmix
        = conc * (1.0 - conc)
          * (l0 + l1 * two_c_minus_one + l2 * two_c_minus_one * two_c_minus_one
                + l3 * two_c_minus_one * two_c_minus_one * two_c_minus_one);

    return fmix;
}

double CALPHADcomputeFMix_derivBinary(const double l0, const double l1,
    const double l2, const double l3, const double conc)
{
    double two_c_minus_one = 2.0 * conc - 1.0;

    double fmix_deriv
        = conc * (1.0 - conc)
              * (2.0 * l1 + 4.0 * l2 * two_c_minus_one
                    + 6.0 * l3 * two_c_minus_one * two_c_minus_one)
          - two_c_minus_one
                * (l0 + l1 * two_c_minus_one
                      + l2 * two_c_minus_one * two_c_minus_one
                      + l3 * two_c_minus_one * two_c_minus_one
                            * two_c_minus_one);

    return fmix_deriv;
}

double CALPHADcomputeFMix_deriv2Binary(const double l0, const double l1,
    const double l2, const double l3, const double conc)
{
    const double conc2 = conc * conc;
    double fmix_deriv2
        = -2.
          * (l0 + l1 * (-3. + 6. * conc) + l2 * (5. - 24. * conc + 24. * conc2)
                + l3 * (-7. + 54. * conc - 120. * conc2 + 80. * conc2 * conc));

    return fmix_deriv2;
}

double CALPHADcomputeFIdealMixBinary(const double rt, const double conc)
{
    double fmix = rt * (xlogx(conc) + xlogx(1.0 - conc));

    return fmix;
}

double CALPHADcomputeFIdealMix_derivBinary(const double rt, const double conc)
{
    double fmix_deriv = rt * (xlogx_deriv(conc) - xlogx_deriv(1.0 - conc));

    return fmix_deriv;
}

double CALPHADcomputeFIdealMix_deriv2Binary(const double rt, const double conc)
{
    double fmix_deriv2 = rt * (xlogx_deriv2(conc) + xlogx_deriv2(1.0 - conc));

    return fmix_deriv2;
}

void CALPHADcomputeFIdealMix_deriv2Ternary(
    const double rt, const double cA, const double cB, double* deriv)
{
    deriv[0] = rt * (xlogx_deriv2(cA) + xlogx_deriv2(1.0 - cA - cB));

    deriv[1] = rt * (xlogx_deriv2(1.0 - cA - cB));

    deriv[2] = deriv[1];

    deriv[3] = rt * (xlogx_deriv2(cB) + xlogx_deriv2(1.0 - cA - cB));
}

double CALPHADcomputeGMix_mixDeriv2(const double l0, const double l1,
    const double l2, const double l3, const double c0, const double c1)
{
    const double dc  = (c0 - c1);
    const double dc2 = dc * dc;

    return l0 + 2. * l1 * dc + l2 * (3. * dc2 - 2. * c0 * c1)
           + l3 * (4. * dc2 * dc - 6. * c0 * c1 * dc);
}

double CALPHADcomputeFMixTernary(const double* lAB, const double* lAC,
    const double* lBC, const double* lABC, const double cA, const double cB)
{
    double cC = 1. - cA - cB;

    double fmix
        = cA * cB
          * (lAB[0] + lAB[1] * (cA - cB) + lAB[2] * (cA - cB) * (cA - cB)
                + lAB[3] * (cA - cB) * (cA - cB) * (cA - cB));

    fmix += cA * cC
            * (lAC[0] + lAC[1] * (cA - cC) + lAC[2] * (cA - cC) * (cA - cC)
                  + lAC[3] * (cA - cC) * (cA - cC) * (cA - cC));

    fmix += cB * cC
            * (lBC[0] + lBC[1] * (cB - cC) + lBC[2] * (cB - cC) * (cB - cC)
                  + lBC[3] * (cB - cC) * (cB - cC) * (cB - cC));

    fmix += cA * cB * cC * (cA * lABC[0] + cB * lABC[1] + cC * lABC[2]);

    return fmix;
}

double CALPHADcomputeFIdealMixTernary(
    const double rt, const double conc0, const double conc1)
{
    double fmix
        = rt * (xlogx(conc0) + xlogx(conc1) + xlogx(1.0 - conc0 - conc1));

    return fmix;
}

void CALPHADcomputeFMix_derivTernary(const double* lAB, const double* lAC,
    const double* lBC, const double* lABC, const double cA, const double cB,
    double* deriv)
{
    double cC = 1. - cA - cB;

    //
    // d/dcA
    //

    // AB terms
    deriv[0] = cB
               * (lAB[0] + lAB[1] * (cA - cB) + lAB[2] * (cA - cB) * (cA - cB)
                     + lAB[3] * (cA - cB) * (cA - cB) * (cA - cB));

    deriv[0] += cA * cB
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    // AC terms
    deriv[0] += (1. - 2. * cA - cB)
                * (lAC[0] + lAC[1] * (cA - cC) + lAC[2] * (cA - cC) * (cA - cC)
                      + lAC[3] * (cA - cC) * (cA - cC) * (cA - cC));

    deriv[0]
        += cA * cC
           * (lAC[1] * 2. + // factor 2. because d(cA-cC)/dcA=2
                 lAC[2] * 4. * (cA - cC) + lAC[3] * 6. * (cA - cC) * (cA - cC));

    // BC terms
    deriv[0] -= cB
                * (lBC[0] + lBC[1] * (cB - cC) + lBC[2] * (cB - cC) * (cB - cC)
                      + lBC[3] * (cB - cC) * (cB - cC) * (cB - cC));

    deriv[0] += cB * cC
                * (lBC[1] + lBC[2] * 2. * (cB - cC)
                      + lBC[3] * 3. * (cB - cC) * (cB - cC));

    //
    // d/dcB
    //

    // AB terms
    deriv[1] = cA
               * (lAB[0] + lAB[1] * (cA - cB) + lAB[2] * (cA - cB) * (cA - cB)
                     + lAB[3] * (cA - cB) * (cA - cB) * (cA - cB));

    deriv[1] -= cA * cB
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    // AC terms
    deriv[1] -= cA
                * (lAC[0] + lAC[1] * (cA - cC) + lAC[2] * (cA - cC) * (cA - cC)
                      + lAC[3] * (cA - cC) * (cA - cC) * (cA - cC));

    deriv[1] += cA * cC
                * (lAC[1] + lAC[2] * 2. * (cA - cC)
                      + lAC[3] * 3. * (cA - cC) * (cA - cC));

    // BC terms
    deriv[1] += (1. - cA - 2. * cB)
                * (lBC[0] + lBC[1] * (cB - cC) + lBC[2] * (cB - cC) * (cB - cC)
                      + lBC[3] * (cB - cC) * (cB - cC) * (cB - cC));

    deriv[1] += cB * cC
                * (lBC[1] * 2. + lBC[2] * 4. * (cB - cC)
                      + lBC[3] * 6. * (cB - cC) * (cB - cC));

    // ABC terms
    const double lphi = cA * lABC[0] + cB * lABC[1] + cC * lABC[2];
    deriv[0]
        += cB * cC * lphi - cA * cB * lphi + cA * cB * cC * (lABC[0] - lABC[2]);
    deriv[1]
        += cA * cC * lphi - cA * cB * lphi + cA * cB * cC * (lABC[1] - lABC[2]);
}

// compute the 4 components of the second order derivative
// with respect to cA and CB
void CALPHADcomputeFMix_deriv2Ternary(const double* lAB, const double* lAC,
    const double* lBC, const double* lABC, const double cA, const double cB,
    double* deriv)
{
    // assert(deriv != 0);

    double cC = 1. - cA - cB;

    //
    // d/dcA*dcA
    //

    // AB term
    deriv[0] = cB
               * (lAB[1] + lAB[2] * 2. * (cA - cB)
                     + lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[0] += cB
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[0] += cA * cB * (lAB[2] * 2. + lAB[3] * 6. * (cA - cB));

    // AC terms
    deriv[0] += -2.
                * (lAC[0] + lAC[1] * (cA - cC) + lAC[2] * (cA - cC) * (cA - cC)
                      + lAC[3] * (cA - cC) * (cA - cC) * (cA - cC));

    deriv[0] += (1. - 2. * cA - cB)
                * (lAC[1] * 2. + lAC[2] * 4. * (cA - cC)
                      + lAC[3] * 6. * (cA - cC) * (cA - cC));

    deriv[0] += (1. - 2. * cA - cB) * // d(cA*cC)/dcA=1-cB-2.*cA
                (lAC[1] * 2. + lAC[2] * 4. * (cA - cC)
                    + lAC[3] * 6. * (cA - cC) * (cA - cC));

    deriv[0] += cA * cC * (lAC[2] * 8. + lAC[3] * 24. * (cA - cC));

    // BC terms
    deriv[0] -= cB
                * (lBC[1] + lBC[2] * 2. * (cB - cC)
                      + lBC[3] * 3. * (cB - cC) * (cB - cC));

    deriv[0] -= cB
                * (lBC[1] + lBC[2] * 2. * (cB - cC)
                      + lBC[3] * 3. * (cB - cC) * (cB - cC));

    deriv[0] += cB * cC * (lBC[2] * 2. + lBC[3] * 6. * (cB - cC));

    //
    // d/dcA*dcB (cross term)
    //

    // AB terms
    deriv[1] = (lAB[0] + lAB[1] * (cA - cB) + lAB[2] * (cA - cB) * (cA - cB)
                + lAB[3] * (cA - cB) * (cA - cB) * (cA - cB));

    deriv[1] -= cB
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[1] += cA
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[1] -= cA * cB * (lAB[2] * 2. + lAB[3] * 6. * (cA - cB));

    // AC terms
    deriv[1] += -1.
                * (lAC[0] + lAC[1] * (cA - cC) + lAC[2] * (cA - cC) * (cA - cC)
                      + lAC[3] * (cA - cC) * (cA - cC) * (cA - cC));

    deriv[1] += (1. - 2. * cA - cB)
                * (lAC[1] + lAC[2] * 2. * (cA - cC)
                      + lAC[3] * 3. * (cA - cC) * (cA - cC));

    deriv[1] += -cA
                * (lAC[1] * 2. + lAC[2] * 4. * (cA - cC)
                      + lAC[3] * 6. * (cA - cC) * (cA - cC));

    deriv[1] += cA * cC * (lAC[2] * 4. + lAC[3] * 12. * (cA - cC));

    // BC tersm
    deriv[1] -= (lBC[0] + lBC[1] * (cB - cC) + lBC[2] * (cB - cC) * (cB - cC)
                 + lBC[3] * (cB - cC) * (cB - cC) * (cB - cC));

    deriv[1] -= cB
                * (lBC[1] * 2. + lBC[2] * 4. * (cB - cC)
                      + lBC[3] * 6. * (cB - cC) * (cB - cC));

    deriv[1] += (1. - cA - 2. * cB)
                * (lBC[1] + lBC[2] * 2. * (cB - cC)
                      + lBC[3] * 3. * (cB - cC) * (cB - cC));

    deriv[1] += cB * cC * (lBC[2] * 4. + lBC[3] * 12. * (cB - cC));

    //
    // d/dcB*dcA (cross term)
    //

    // AB terms
    deriv[2] = (lAB[0] + lAB[1] * (cA - cB) + lAB[2] * (cA - cB) * (cA - cB)
                + lAB[3] * (cA - cB) * (cA - cB) * (cA - cB));

    deriv[2] += cA
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[2] -= cB
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[2] -= cA * cB * (lAB[2] * 2. + lAB[3] * 6. * (cA - cB));

    // AC terms
    deriv[2] += -1.
                * (lAC[0] + lAC[1] * (cA - cC) + lAC[2] * (cA - cC) * (cA - cC)
                      + lAC[3] * (cA - cC) * (cA - cC) * (cA - cC));

    deriv[2] += (1. - 2. * cA - cB)
                * (lAC[1] + lAC[2] * 2. * (cA - cC)
                      + lAC[3] * 3. * (cA - cC) * (cA - cC));

    deriv[2] += -cA
                * (lAC[1] * 2. + lAC[2] * 4. * (cA - cC)
                      + lAC[3] * 6. * (cA - cC) * (cA - cC));

    deriv[2] += cA * cC * (lAC[2] * 4. + lAC[3] * 12. * (cA - cC));

    // BC terms
    deriv[2] -= (lBC[0] + lBC[1] * (cB - cC) + lBC[2] * (cB - cC) * (cB - cC)
                 + lBC[3] * (cB - cC) * (cB - cC) * (cB - cC));

    deriv[2] -= cB
                * (lBC[1] * 2. + lBC[2] * 4. * (cB - cC)
                      + lBC[3] * 6. * (cB - cC) * (cB - cC));

    deriv[2] += (1. - 2. * cB - cA)
                * (lBC[1] + lBC[2] * 2. * (cB - cC)
                      + lBC[3] * 3. * (cB - cC) * (cB - cC));

    deriv[2] += cB * cC * (lBC[2] * 4. + lBC[3] * 12. * (cB - cC));

    //
    // d/dcB*dcB
    //

    // AB terms
    deriv[3] = cA
               * (-lAB[1] + -lAB[2] * 2. * (cA - cB)
                     + -lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[3] -= cA
                * (lAB[1] + lAB[2] * 2. * (cA - cB)
                      + lAB[3] * 3. * (cA - cB) * (cA - cB));

    deriv[3] += cA * cB * (lAB[2] * 2. + lAB[3] * 6. * (cA - cB));

    // AC terms
    deriv[3] += -cA
                * (lAC[1] + lAC[2] * 2. * (cA - cC)
                      + lAC[3] * 3. * (cA - cC) * (cA - cC));

    deriv[3] += -cA
                * (lAC[1] + lAC[2] * 2. * (cA - cC)
                      + lAC[3] * 3. * (cA - cC) * (cA - cC));

    deriv[3] += cA * cC * (lAC[2] * 2. + lAC[3] * 6. * (cA - cC));

    // BC terms
    deriv[3] += -2.
                * (lBC[0] + lBC[1] * (cB - cC) + lBC[2] * (cB - cC) * (cB - cC)
                      + lBC[3] * (cB - cC) * (cB - cC) * (cB - cC));

    deriv[3] += (1. - cA - 2. * cB)
                * (lBC[1] * 2. + lBC[2] * 4. * (cB - cC)
                      + lBC[3] * 6. * (cB - cC) * (cB - cC));

    deriv[3] += (1. - cA - 2. * cB)
                * (lBC[1] * 2. + lBC[2] * 4. * (cB - cC)
                      + lBC[3] * 6. * (cB - cC) * (cB - cC));

    deriv[3] += cB * cC * (lBC[2] * 8. + lBC[3] * 24. * (cB - cC));

    // ABC terms
    const double lphi = cA * lABC[0] + cB * lABC[1] + cC * lABC[2];
    deriv[0]
        += -2. * cB * lphi + 2. * (cB * cC - cA * cB) * (lABC[0] - lABC[2]);
    deriv[1] += (cC - cB - cA) * lphi + cB * (cC - cA) * (lABC[1] - lABC[2])
                + cA * (cC - cB) * (lABC[0] - lABC[2]);
    deriv[2] += (cC - cA - cB) * lphi + cA * (cC - cB) * (lABC[0] - lABC[2])
                + cB * (cC - cA) * (lABC[1] - lABC[2]);
    deriv[3]
        += -2. * cA * lphi + 2. * (cA * cC - cA * cB) * (lABC[1] - lABC[2]);
}

void CALPHADcomputeFIdealMix_derivTernary(
    const double rt, const double cA, const double cB, double* deriv)
{
    deriv[0] = rt * (xlogx_deriv(cA) - xlogx_deriv(1.0 - cA - cB));

    deriv[1] = rt * (xlogx_deriv(cB) - xlogx_deriv(1.0 - cA - cB));
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

void readLmixBinary(pt::ptree& db, double LmixPhase[4][MAX_POL_T_INDEX])
{
    {
        int i = 0;
        for (pt::ptree::value_type& v : db.get_child("L0"))
        {
            LmixPhase[0][i] = v.second.get_value<double>();
            i++;
        }
        if (i < MAX_POL_T_INDEX) LmixPhase[0][i] = 0.;
    }

    {
        int i = 0;
        for (pt::ptree::value_type& v : db.get_child("L1"))
        {
            LmixPhase[1][i] = v.second.get_value<double>();
            i++;
        }
        if (i < MAX_POL_T_INDEX) LmixPhase[1][i] = 0.;
    }

    // L2
    {
        auto child = db.get_child_optional("L2");
        if (child)
        {
            int i = 0;
            for (pt::ptree::value_type& v : db.get_child("L2"))
            {
                LmixPhase[2][i] = v.second.get_value<double>();
                i++;
            }
            if (i < MAX_POL_T_INDEX) LmixPhase[2][i] = 0.;
        }
        else
        {
            for (int i = 0; i < MAX_POL_T_INDEX; i++)
                LmixPhase[2][i] = 0.0;
        }
    }

    // L3
    {
        auto child = db.get_child_optional("L3");
        if (child)
        {
            int i = 0;
            for (pt::ptree::value_type& v : db.get_child("L3"))
            {
                LmixPhase[3][i] = v.second.get_value<double>();
                i++;
            }
            if (i < MAX_POL_T_INDEX) LmixPhase[3][i] = 0.;
        }
        else
        {
            for (int i = 0; i < MAX_POL_T_INDEX; i++)
                LmixPhase[3][i] = 0.0;
        }
    }
}
}
