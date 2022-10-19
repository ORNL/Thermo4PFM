#include "CALPHADFunctions.h"
#include "xlogx.h"

//#include <cassert>

namespace pt = boost::property_tree;

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
double CALPHADcomputeFMixBinary(const CalphadDataType l0,
    const CalphadDataType l1, const CalphadDataType l2,
    const CalphadDataType l3, const double conc)
{
    double two_c_minus_one = 2.0 * conc - 1.0;

    double fmix
        = conc * (1.0 - conc)
          * (l0 + l1 * two_c_minus_one + l2 * two_c_minus_one * two_c_minus_one
                + l3 * two_c_minus_one * two_c_minus_one * two_c_minus_one);

    return fmix;
}

double CALPHADcomputeFMix_derivBinary(const CalphadDataType l0,
    const CalphadDataType l1, const CalphadDataType l2,
    const CalphadDataType l3, const double conc)
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

double CALPHADcomputeFMix_deriv2Binary(const CalphadDataType l0,
    const CalphadDataType l1, const CalphadDataType l2,
    const CalphadDataType l3, const double conc)
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

template <typename DataType>
void CALPHADcomputeFIdealMix_deriv2Ternary(
    const DataType rt, const DataType cA, const DataType cB, DataType* deriv)
{
    const double tmp = xlogx_deriv2<DataType>(1.0 - cA - cB);
    deriv[0]         = rt * (xlogx_deriv2<DataType>(cA) + tmp);

    deriv[1] = rt * tmp;

    deriv[2] = deriv[1];

    deriv[3] = rt * (xlogx_deriv2<DataType>(cB) + tmp);
}

double CALPHADcomputeGMix_mixDeriv2(const CalphadDataType l0,
    const CalphadDataType l1, const CalphadDataType l2,
    const CalphadDataType l3, const double c0, const double c1)
{
    const double dc  = (c0 - c1);
    const double dc2 = dc * dc;

    return l0 + 2. * l1 * dc + l2 * (3. * dc2 - 2. * c0 * c1)
           + l3 * (4. * dc2 * dc - 6. * c0 * c1 * dc);
}

#define OMEGAIJ(omega, ci, cj)                                                 \
    (omega[0]                                                                  \
        + (ci - cj)                                                            \
              * (omega[1] + (ci - cj) * (omega[2] + (ci - cj) * omega[3])))

#define DOMEGAIJDXI(omega, ci, cj)                                             \
    (omega[1] + (ci - cj) * (2. * omega[2] + 3. * omega[3] * (ci - cj)))

#define DOMEGAIJDX2(omega, ci, cj) (2. * omega[2] + 6. * omega[3] * (ci - cj))

double CALPHADcomputeFMixTernary(const CalphadDataType lAB[4],
    const CalphadDataType lAC[4], const CalphadDataType lBC[4],
    const CalphadDataType lABC[3], const double cA, const double cB)
{
    double cC = 1. - cA - cB;

    double fmix = cA * cB * OMEGAIJ(lAB, cA, cB)
                  + cA * cC * OMEGAIJ(lAC, cA, cC)
                  + cB * cC * OMEGAIJ(lBC, cB, cC)
                  + cA * cB * cC * (cA * lABC[0] + cB * lABC[1] + cC * lABC[2]);

    return fmix;
}

double CALPHADcomputeFIdealMixTernary(
    const double rt, const double conc0, const double conc1)
{
    double fmix
        = rt * (xlogx(conc0) + xlogx(conc1) + xlogx(1.0 - conc0 - conc1));

    return fmix;
}

void CALPHADcomputeFMix_derivTernary(const CalphadDataType lAB[4],
    const CalphadDataType lAC[4], const CalphadDataType lBC[4],
    const CalphadDataType lABC[3], const double cA, const double cB,
    double deriv[2])
{
    double cC = 1. - cA - cB;

    double omegaAB = OMEGAIJ(lAB, cA, cB);
    double omegaAC = OMEGAIJ(lAC, cA, cC);
    double omegaBC = OMEGAIJ(lBC, cB, cC);

    double domegaABdcA = DOMEGAIJDXI(lAB, cA, cB);
    double domegaACdcA = 2. * DOMEGAIJDXI(lAC, cA, cC);
    double domegaBCdcB = 2. * DOMEGAIJDXI(lBC, cB, cC);

    // 0 -> d/dcA
    // 1 -> d/dcB

    // AB terms
    deriv[0] = cB * (omegaAB + cA * domegaABdcA);
    deriv[1] = cA * (omegaAB - cB * domegaABdcA);

    // AC terms
    deriv[0] += (cC - cA) * omegaAC;
    deriv[0] += cA * cC * domegaACdcA;

    deriv[1] += cA * (cC * 0.5 * domegaACdcA - omegaAC);

    // BC terms
    deriv[1] += (cC - cB) * omegaBC;
    deriv[1] += cB * cC * domegaBCdcB;

    deriv[0] += cB * (cC * 0.5 * domegaBCdcB - omegaBC);

    // ABC terms
    const double lphi = cA * lABC[0] + cB * lABC[1] + cC * lABC[2];
    deriv[0]
        += cB * cC * lphi - cA * cB * lphi + cA * cB * cC * (lABC[0] - lABC[2]);
    deriv[1]
        += cA * cC * lphi - cA * cB * lphi + cA * cB * cC * (lABC[1] - lABC[2]);
}

// compute the 4 components of the second order derivative
// with respect to cA and CB
void CALPHADcomputeFMix_deriv2Ternary(const CalphadDataType lAB[4],
    const CalphadDataType lAC[4], const CalphadDataType lBC[4],
    const CalphadDataType lABC[3], const double cA, const double cB,
    double deriv[4])
{
    // assert(deriv != 0);

    double cC = 1. - cA - cB;

    double omegaAB = OMEGAIJ(lAB, cA, cB);
    double omegaAC = OMEGAIJ(lAC, cA, cC);
    double omegaBC = OMEGAIJ(lBC, cB, cC);

    double domegaABdcA = DOMEGAIJDXI(lAB, cA, cB);
    double domegaACdcA = 2. * DOMEGAIJDXI(lAC, cA, cC);
    double domegaBCdcB = 2. * DOMEGAIJDXI(lBC, cB, cC);

    double d2omegaABdcA2 = DOMEGAIJDX2(lAB, cA, cB);
    double d2omegaACdcA2 = 4. * DOMEGAIJDX2(lAC, cA, cC);
    double d2omegaBCdcB2 = 4. * DOMEGAIJDX2(lBC, cB, cC);

    //
    // d/dcA*dcA
    //

    // AB term
    deriv[0] = cB * (2. * domegaABdcA + cA * d2omegaABdcA2);

    // AC terms
    deriv[0] -= 2. * omegaAC;
    deriv[0] += 2. * (cC - cA) * domegaACdcA;
    deriv[0] += cA * cC * d2omegaACdcA2;

    // BC terms
    deriv[0] += cB * (cC * 0.25 * d2omegaBCdcB2 - domegaBCdcB);

    //
    // d/dcA*dcB (cross term)
    //

    // AB terms
    deriv[1] = omegaAB;
    deriv[1] += (cA - cB) * domegaABdcA;
    deriv[1] -= cA * cB * d2omegaABdcA2;

    // AC terms
    deriv[1] -= omegaAC;
    deriv[1] += (cC - 3. * cA) * 0.5 * domegaACdcA;
    deriv[1] += cA * cC * 0.5 * d2omegaACdcA2;

    // BC terms
    deriv[1] -= omegaBC;
    deriv[1] += (cC - 3. * cB) * 0.5 * domegaBCdcB;
    deriv[1] += cB * cC * 0.5 * d2omegaBCdcB2;

    //
    // d/dcB*dcA (cross term)
    //

    deriv[2] = deriv[1];

    //
    // d/dcB*dcB
    //

    // AB terms
    deriv[3] = cA * (cB * d2omegaABdcA2 - 2. * domegaABdcA);

    // AC terms
    deriv[3] += cA * (cC * 0.25 * d2omegaACdcA2 - domegaACdcA);

    // BC terms
    deriv[3] -= 2. * omegaBC;
    deriv[3] += 2. * (cC - cB) * domegaBCdcB;
    deriv[3] += cB * cC * d2omegaBCdcB2;

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
    const double alpha = xlogx_deriv(1.0 - cA - cB);

    deriv[0] = rt * (xlogx_deriv(cA) - alpha);
    deriv[1] = rt * (xlogx_deriv(cB) - alpha);
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

void readLmixBinary(
    pt::ptree& db, CalphadDataType LmixPhase[4][MAX_POL_T_INDEX])
{
    {
        int i = 0;
        for (pt::ptree::value_type& v : db.get_child("L0"))
        {
            LmixPhase[0][i] = v.second.get_value<CalphadDataType>();
            i++;
        }
        if (i < MAX_POL_T_INDEX) LmixPhase[0][i] = 0.;
    }

    {
        int i = 0;
        for (pt::ptree::value_type& v : db.get_child("L1"))
        {
            LmixPhase[1][i] = v.second.get_value<CalphadDataType>();
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
                LmixPhase[2][i] = v.second.get_value<CalphadDataType>();
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
                LmixPhase[3][i] = v.second.get_value<CalphadDataType>();
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

void readSublatticeStoichiometry(
    boost::property_tree::ptree& db, int sublatticeStoichiometryPhase[2])
{
    if (db.get_child_optional("p"))
    {
        sublatticeStoichiometryPhase[0] = db.get_child("p").get_value<int>();
    }
    else
    {
        sublatticeStoichiometryPhase[0] = 0;
    }

    if (db.get_child_optional("q"))
    {
        sublatticeStoichiometryPhase[1] = db.get_child("q").get_value<int>();
    }
    else
    {
        sublatticeStoichiometryPhase[1] = 1;
    }
}

bool checkSingleSublattice(boost::property_tree::ptree& db)
{
    int p = 0;
    int q = 1;

    if (db.get_child_optional("p"))
    {
        p = db.get_child("p").get_value<int>();
    }

    if (db.get_child_optional("q"))
    {
        q = db.get_child("q").get_value<int>();
    }

    if (p == 0 && q == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool checkSublatticeSpecies(boost::property_tree::ptree& species_db)
{
    auto phase_db = species_db.get_child("PhaseL");
    if (phase_db.get_child_optional("p")) return true;
    if (phase_db.get_child_optional("q")) return true;
    phase_db = species_db.get_child("PhaseA");
    if (phase_db.get_child_optional("p")) return true;
    if (phase_db.get_child_optional("q")) return true;
    if (species_db.get_child_optional("PhaseB"))
    {
        phase_db = species_db.get_child("PhaseB");
        if (phase_db.get_child_optional("p")) return true;
        if (phase_db.get_child_optional("q")) return true;
    }

    return false;
}

// check if database contains sublattice parameters "p" or "q"
bool checkSublattice(boost::property_tree::ptree& db)
{
    boost::property_tree::ptree& speciesA_db = db.get_child("SpeciesA");
    if (checkSublatticeSpecies(speciesA_db)) return true;

    boost::property_tree::ptree& speciesB_db = db.get_child("SpeciesB");
    if (checkSublatticeSpecies(speciesB_db)) return true;

    return false;
}

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
template void CALPHADcomputeFIdealMix_deriv2Ternary<float>(
    const float rt, const float cA, const float cB, float* deriv);
template void CALPHADcomputeFIdealMix_deriv2Ternary<double>(
    const double rt, const double cA, const double cB, double* deriv);
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
