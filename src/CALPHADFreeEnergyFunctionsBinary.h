#ifndef included_CALPHADFreeEnergyFunctionsBinary
#define included_CALPHADFreeEnergyFunctionsBinary

#include "CALPHADFreeEnergyFunctions.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
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

class CALPHADFreeEnergyFunctionsBinary : public CALPHADFreeEnergyFunctions
{
public:
    CALPHADFreeEnergyFunctionsBinary(boost::property_tree::ptree& input_db,
        boost::optional<boost::property_tree::ptree&> newton_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    ~CALPHADFreeEnergyFunctionsBinary(){};

    virtual double computeFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, const bool gp = false);
    virtual void computeDerivFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, double*);
    virtual void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2);

    virtual bool computeCeqT(const double temperature, double* ceq,
        const int maxits = 20, const bool verbose = false);

    void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.);

    int computePhaseConcentrations(const double temperature, const double* conc,
        const double phi, double* x);
    void energyVsPhiAndC(const double temperature, const double* const ceq,
        const bool found_ceq, const double phi_well_scale,
        const int npts_phi = 51,
        const int npts_c   = 50); // # of compositions to use (>1)
    void printEnergyVsComposition(
        const double temperature, const int npts = 100);
    double fchem(
        const double phi, const double* const conc, const double temperature);
    void printEnergyVsPhiHeader(const double temperature, const int nphi,
        const int nc, const double cmin, const double cmax, const double slopec,
        std::ostream& os) const;
    void printEnergyVsPhi(const double* const conc, const double temperature,
        const double phi_well_scale, const int npts, const double slopec,
        std::ostream& os);

private:
    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    void readNewtonparameters(boost::property_tree::ptree& newton_db);

    void computeTdependentParameters(const double temperature, double* Lmix_L,
        double* Lmix_A, double* fA, double* fB);

    std::string fenergy_diag_filename_;

    double newton_tol_;
    double newton_alpha_;
    int newton_maxits_;
    bool newton_verbose_;

    // size 2 for species 0 and 1
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseL_[2];
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseA_[2];

    // size 4 for L0, L1, L2, L3,
    // can contain up to 3 coefficients a,b,c for a+b*T,
    // possibly +c*T*ln(T) if compiled with -DLMIX_WTLOGT
    double LmixPhaseL_[4][MAX_POL_T_INDEX];
    double LmixPhaseA_[4][MAX_POL_T_INDEX];

    double (*fun_ptr_arr_[3])(const double){ linear_interp_func,
        pbg_interp_func, harmonic_interp_func };

    void readParameters(boost::property_tree::ptree& calphad_db);

    // energy of species "is" in phase L,A
    double getFenergyPhaseL(const short is, const double temperature)
    {
        return g_species_phaseL_[is].fenergy(temperature);
    }
    double getFenergyPhaseA(const short is, const double temperature)
    {
        return g_species_phaseA_[is].fenergy(temperature);
    }

    double lmixPhase(
        const unsigned index, const PhaseIndex pi, const double temperature)
    {
        assert(index < 4);

        switch (pi)
        {
            case PhaseIndex::phaseL:
                return LmixPhaseL_[index][0]
                       + LmixPhaseL_[index][1] * temperature
#ifdef LMIX_WTLOGT
                       + LmixPhaseL_[index][2] * temperature * log(temperature)
#endif
                    ;
            case PhaseIndex::phaseA:
                return LmixPhaseA_[index][0]
                       + LmixPhaseA_[index][1] * temperature
#ifdef LMIX_WTLOGT
                       + LmixPhaseA_[index][2] * temperature * log(temperature)
#endif
                    ;
            default:
                return NAN;
        }
    }

    void computePhasesFreeEnergies(const double temperature, const double hphi,
        const double conc, double& fl, double& fa);
};
}
#endif
