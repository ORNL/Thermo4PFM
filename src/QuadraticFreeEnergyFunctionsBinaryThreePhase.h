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

    bool computeCeqT(const double temperature, double* ceq,
        const int maxits = 20, const bool verbose = false);

    void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.);

    int computePhaseConcentrations(const double temperature, const double* conc,
        const double* const phi, double* x);
    void energyVsPhiAndC(const double temperature, const double* const ceq,
        const bool found_ceq, const double phi_well_scale,
        const int npts_phi = 51,
        const int npts_c   = 50); // # of compositions to use (>1)
    void printEnergyVsComposition(const double temperature, std::ostream& os,
        const double cmin, const double cmax, const int npts = 100);
    double fchem(const double* const phi, const double* const conc,
        const double temperature);
    void printEnergyVsPhiHeader(const double temperature, const int nphi,
        const int nc, const double cmin, const double cmax, const double slopec,
        std::ostream& os) const;
    void printEnergyVsPhi(const double* const conc, const double temperature,
        const double phi_well_scale, const int npts, const double slopec,
        std::ostream& os);

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
