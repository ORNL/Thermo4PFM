#ifndef Thermo4PFM_included_QuadraticFreeEnergyFunctionsBinary
#define Thermo4PFM_included_QuadraticFreeEnergyFunctionsBinary

#include "InterpolationType.h"
#include "Phases.h"
#include "functions.h"

#include <boost/property_tree/ptree.hpp>

#include <cassert>
#include <fstream>
#include <iostream>
#include <math.h>

namespace Thermo4PFM
{

class QuadraticFreeEnergyFunctionsBinary
{
public:
    QuadraticFreeEnergyFunctionsBinary(const double Tref, const double Al,
        const double ceql, const double m_liquid, const double Aa,
        const double ceqa, const double m_solid,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    ~QuadraticFreeEnergyFunctionsBinary() { delete[] fenergy_diag_filename_; };

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
    char* fenergy_diag_filename_;

    const double Tref_;

    const double Al_;
    const double ceql_;
    const double m_liquid_;

    const double Aa_;
    const double ceqa_;
    const double m_solid_;

    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    void computePhasesFreeEnergies(const double temperature,
        const double* const hphi, const double conc, double& fl, double& fa);

    double computeLiquidConcentration(
        const double temp, const double hphi, const double c) const;
    double computeSolidAConcentration(
        const double temp, const double hphi, const double c) const;
};
}
#endif
