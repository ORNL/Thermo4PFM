#ifndef included_CALPHADFunctions
#define included_CALPHADFunctions

#include <boost/property_tree/ptree.hpp>

namespace Thermo4PFM
{
double CALPHADcomputeFMixBinary(const double l0, const double l1,
    const double l2, const double l3, const double conc);
double CALPHADcomputeFMix_derivBinary(const double l0, const double l1,
    const double l2, const double l3, const double conc);
double CALPHADcomputeFMix_deriv2Binary(const double l0, const double l1,
    const double l2, const double l3, const double conc);
double CALPHADcomputeFIdealMixBinary(const double rt, const double conc);
double CALPHADcomputeFIdealMix_derivBinary(const double rt, const double conc);
double CALPHADcomputeFIdealMix_deriv2Binary(const double rt, const double conc);
double CALPHADcomputeFMixTernary(const double* lAB, const double* lAC,
    const double* lBC, const double* lABC, const double cA, const double cB);
double CALPHADcomputeFIdealMixTernary(
    const double rt, const double conc0, const double conc1);
void CALPHADcomputeFIdealMix_derivTernary(
    const double rt, const double cA, const double cB, double* deriv);
void CALPHADcomputeFIdealMix_deriv2Ternary(
    const double rt, const double cA, const double cB, double* deriv);

void CALPHADcomputeFMix_derivTernary(const double* lAB, const double* lAC,
    const double* lBC, const double* lABC, const double cA, const double cB,
    double* deriv);
void CALPHADcomputeFMix_deriv2Ternary(const double* lAB, const double* lAC,
    const double* lBC, const double* lABC, const double cA, const double cB,
    double* deriv);
void readLmixBinary(
    boost::property_tree::ptree& db, double LmixPhase[4][MAX_POL_T_INDEX]);
}

#endif
