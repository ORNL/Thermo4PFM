#ifndef Thermo4PFM_included_CALPHADFunctions
#define Thermo4PFM_included_CALPHADFunctions

#include "datatypes.h"

#include <boost/property_tree/ptree.hpp>

namespace Thermo4PFM
{
double CALPHADcomputeFMixBinary(const CalphadDataType l0,
    const CalphadDataType l1, const CalphadDataType l2,
    const CalphadDataType l3, const double conc);
double CALPHADcomputeFMix_derivBinary(const CalphadDataType l0,
    const CalphadDataType l1, const CalphadDataType l2,
    const CalphadDataType l3, const double conc);
double CALPHADcomputeFMix_deriv2Binary(const CalphadDataType l0,
    const CalphadDataType l1, const CalphadDataType l2,
    const CalphadDataType l3, const double conc);
double CALPHADcomputeFIdealMixBinary(const double rt, const double conc);
double CALPHADcomputeFIdealMix_derivBinary(const double rt, const double conc);
double CALPHADcomputeFIdealMix_deriv2Binary(const double rt, const double conc);

double CALPHADcomputeFMixTernary(const CalphadDataType* lAB,
    const CalphadDataType* lAC, const CalphadDataType* lBC,
    const CalphadDataType* lABC, const double cA, const double cB);
double CALPHADcomputeFIdealMixTernary(
    const double rt, const double conc0, const double conc1);
void CALPHADcomputeFIdealMix_derivTernary(
    const double rt, const double cA, const double cB, double* deriv);
template <typename DataType>
void CALPHADcomputeFIdealMix_deriv2Ternary(
    const DataType rt, const DataType cA, const DataType cB, DataType* deriv);

void CALPHADcomputeFMix_derivTernary(const CalphadDataType* lAB,
    const CalphadDataType* lAC, const CalphadDataType* lBC,
    const CalphadDataType* lABC, const double cA, const double cB,
    double* deriv);
void CALPHADcomputeFMix_deriv2Ternary(const CalphadDataType* lAB,
    const CalphadDataType* lAC, const CalphadDataType* lBC,
    const CalphadDataType* lABC, const double cA, const double cB,
    double* deriv);
void readLmixBinary(boost::property_tree::ptree& db,
    CalphadDataType LmixPhase[4][MAX_POL_T_INDEX]);
void readSublatticeStoichiometry(
    boost::property_tree::ptree& db, int sublatticeStoichiometryPhase[2]);
bool checkSingleSublattice(boost::property_tree::ptree& db);
bool checkSublattice(boost::property_tree::ptree& db);
}

#endif
