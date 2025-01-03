#ifndef Thermo4PFM_included_ParabolicFreeEnergyFunctionsBinary
#define Thermo4PFM_included_ParabolicFreeEnergyFunctionsBinary

#include "InterpolationType.h"
#include "Phases.h"
#include "datatypes.h"
#include "functions.h"

#include <boost/property_tree/ptree.hpp>

#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>

namespace Thermo4PFM
{

class ParabolicFreeEnergyFunctionsBinary
{
public:
    ParabolicFreeEnergyFunctionsBinary(const double Tref,
        const double coeffL[][2], const double coeffA[][2],
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    ~ParabolicFreeEnergyFunctionsBinary() { delete[] fenergy_diag_filename_; };

    double computeFreeEnergy(const double temperature, const double* const conc,
        const PhaseIndex pi, const bool gp = false);
    void computeDerivFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, double*);
    void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2);

    int computePhaseConcentrations(const double temperature, const double* conc,
        const double* const phi, double* x);
    void printEnergyVsComposition(const double temperature, std::ostream& os,
        const double cmin, const double cmax, const int npts = 100);
    double fchem(const double* const phi, const double* const conc,
        const double temperature);
    bool computeCeqT(const double temperature, double* ceq, const int maxits,
        const bool verbose);

private:
    double Tref_;
    /*
     * Phase L coefficients
     * g(c,T) = (aL1_*T+aL0_)*c*c + (bL1_*T+bL0_)*c + cL1_*T+cL0_
     */
    double coeffL_[3][2];

    /*
     * Phase A coefficients
     * g(c,T) = (aL1_*T+aL0_)*c*c + (bL1_*T+bL0_)*c + cL1_*T+cL0_
     */
    double coeffA_[3][2];

    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    void getPhaseCoeffs(const PhaseIndex pi, double& a0, double& a1, double& b0,
        double& b1, double& c0, double& c1);

    char* fenergy_diag_filename_;

    void computePhasesFreeEnergies(const double temperature,
        const double* const hphi, const double conc, double& fl, double& fa);
};
}
#endif
