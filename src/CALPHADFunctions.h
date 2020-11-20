#ifndef included_CALPHADFunctions
#define included_CALPHADFunctions

#include <vector>

double CALPHADcomputeFMixBinary(const double l0, const double l1,
    const double l2, const double l3, const double conc);
double CALPHADcomputeFMix_derivBinary(const double l0, const double l1,
    const double l2, const double l3, const double conc);
double CALPHADcomputeFMix_deriv2Binary(const double l0, const double l1,
    const double l2, const double l3, const double conc);
double CALPHADcomputeFIdealMixBinary(const double rt, const double conc);
double CALPHADcomputeFIdealMix_derivBinary(const double rt, const double conc);
double CALPHADcomputeFIdealMix_deriv2Binary(const double rt, const double conc);
void CALPHADcomputeFIdealMix_deriv2(const double rt,
    const std::vector<double>& conc, std::vector<double>& d2fdc2);
double CALPHADcomputeGMix_deriv2(const double l1, const double l2,
    const double l3, const std::vector<double>& conc, const int ic);
double CALPHADcomputeGMix_mixDeriv2(const double l0, const double l1,
    const double l2, const double l3, const std::vector<double>& conc,
    const int ic0, const int ic1);
double CALPHADcomputeFMix_mixDeriv2(const double l0, const double l1,
    const double l2, const double l3, const std::vector<double>& concf,
    const int ic0, const int ic1);

double CALPHADcomputePenalty(const double alpha1, const double p12,
    const double p13, const double alpha2, const double p22, const double p23,
    const double conc);
double CALPHADcomputeDerivPenalty(const double alpha1, const double p12,
    const double p13, const double alpha2, const double p22, const double p23,
    const double conc);
double CALPHADcompute2ndDerivPenalty(const double alpha1, const double p12,
    const double p13, const double alpha2, const double p22, const double p23,
    const double conc);

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
#endif
