#ifndef Thermo4PFM_included_QuadraticFreeEnergyFunctionsTernaryThreePhase
#define Thermo4PFM_included_QuadraticFreeEnergyFunctionsTernaryThreePhase

#include "InterpolationType.h"
#include "Phases.h"
#include "datatypes.h"
#include "functions.h"

#include <fstream>
#include <iostream>

namespace Thermo4PFM
{

class QuadraticFreeEnergyFunctionsTernaryThreePhase
{
public:
    QuadraticFreeEnergyFunctionsTernaryThreePhase(const double Al[2],
        const double ceql[2], const double Aa[2], const double ceqa[2],
        const double Ab[2], const double ceqb[2],
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    ~QuadraticFreeEnergyFunctionsTernaryThreePhase(){};

    double computeFreeEnergy(const double temperature, const double* const conc,
        const PhaseIndex pi, const bool gp = false);
    void computeDerivFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, double* deriv);
    void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2);


    int computePhaseConcentrations(const double temperature,
        const double* const conc, const double* const phi, double* x);
    double fchem(const double* const phi, const double* const conc,
        const double temperature);

private:
    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    const double Al_[2];
    const double Aa_[2];
    const double Ab_[2];

    const double ceql_[2];
    const double ceqa_[2];
    const double ceqb_[2];

    void computePhasesFreeEnergies(const double temperature,
        const double* const hphi, const double conc0, const double conc1,
        double& fl, double& fa, double& fb);
};
}
#endif
