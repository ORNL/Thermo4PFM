#ifndef Thermo4PFM_included_QuadraticFreeEnergyFunctionsBinaryThreePhase
#define Thermo4PFM_included_QuadraticFreeEnergyFunctionsBinaryThreePhase

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

class QuadraticFreeEnergyFunctionsBinaryThreePhase
{
public:
    QuadraticFreeEnergyFunctionsBinaryThreePhase(const double Al,
        const double ceql, const double Aa, const double ceqa, const double Ab,
        const double ceqb,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    ~QuadraticFreeEnergyFunctionsBinaryThreePhase()
    {
        delete[] fenergy_diag_filename_;
    };

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

private:
    double Al_;
    double ceql_;

    double Aa_;
    double ceqa_;

    double Ab_;
    double ceqb_;

    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    char* fenergy_diag_filename_;

    void computePhasesFreeEnergies(const double temperature,
        const double* const hphi, const double conc, double& fl, double& fa,
        double& fb);
};
}
#endif
