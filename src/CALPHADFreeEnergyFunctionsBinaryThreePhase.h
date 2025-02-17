#ifndef Thermo4PFM_included_CALPHADFreeEnergyFunctionsBinaryThreePhase
#define Thermo4PFM_included_CALPHADFreeEnergyFunctionsBinaryThreePhase

#include "CALPHADSpeciesPhaseGibbsEnergy.h"
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

class CALPHADFreeEnergyFunctionsBinaryThreePhase
{
public:
    CALPHADFreeEnergyFunctionsBinaryThreePhase(
        boost::property_tree::ptree& input_db,
        boost::optional<boost::property_tree::ptree&> newton_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    virtual ~CALPHADFreeEnergyFunctionsBinaryThreePhase()
    {
        delete[] fenergy_diag_filename_;
    };

    double computeFreeEnergy(const double temperature, const double* const conc,
        const PhaseIndex pi, const bool gp = false);
    void computeDerivFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, double*);
    void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2);

    void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.);

    virtual int computePhaseConcentrations(const double temperature,
        const double* conc, const double* const phi, double* x);
    void printEnergyVsComposition(const double temperature, std::ostream& os,
        const double cmin, const double cmax, const int npts = 100);
    double fchem(const double* const phi, const double* const conc,
        const double temperature);

    void computeTdependentParameters(const double temperature,
        CalphadDataType* Lmix_L, CalphadDataType* Lmix_A,
        CalphadDataType* Lmix_B, CalphadDataType* fA, CalphadDataType* fB);

protected:
    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    double newton_tol_;
    double newton_alpha_;
    int newton_maxits_;
    bool newton_verbose_;

private:
    void readNewtonparameters(boost::property_tree::ptree& newton_db);

    char* fenergy_diag_filename_;

    // Single species energies in each phase
    // size 2 for species 0 and 1
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseL_[2];
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseA_[2];
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseB_[2];

    // size 4 for L0, L1, L2, L3,
    // can contain up to 3 coefficients a,b,c for a+b*T,
    // possibly +c*T*ln(T) if compiled with -DLMIX_WTLOGT
    CalphadDataType LmixPhaseL_[4][MAX_POL_T_INDEX];
    CalphadDataType LmixPhaseA_[4][MAX_POL_T_INDEX];
    CalphadDataType LmixPhaseB_[4][MAX_POL_T_INDEX];

    void readParameters(boost::property_tree::ptree& calphad_db);

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
    // energy of species "is" in phase L,A
    double getFenergyPhaseL(const short is, const double temperature)
    {
        return g_species_phaseL_[is].fenergy(temperature);
    }
    double getFenergyPhaseA(const short is, const double temperature)
    {
        return g_species_phaseA_[is].fenergy(temperature);
    }
    double getFenergyPhaseB(const short is, const double temperature)
    {
        return g_species_phaseB_[is].fenergy(temperature);
    }

    CalphadDataType lmixPhase(
        const unsigned index, const PhaseIndex pi, const double temperature)
    {
        // assert(index < 4);

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
            case PhaseIndex::phaseB:
                return LmixPhaseB_[index][0]
                       + LmixPhaseB_[index][1] * temperature
#ifdef LMIX_WTLOGT
                       + LmixPhaseB_[index][2] * temperature * log(temperature)
#endif
                    ;
            default:
                return NAN;
        }
    }
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

    void computePhasesFreeEnergies(const double temperature,
        const double* const hphi, const double conc, double& fl, double& fa,
        double& fb);
};
}
#endif
