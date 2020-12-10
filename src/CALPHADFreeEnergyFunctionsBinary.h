#ifndef included_CALPHADFreeEnergyFunctionsBinary
#define included_CALPHADFreeEnergyFunctionsBinary

#include "CALPHADConcSolverBinary.h"
#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADFreeEnergyFunctions.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "InterpolationType.h"
#include "Phases.h"

#include <boost/property_tree/ptree.hpp>

#include <cassert>
#include <fstream>
#include <iostream>

namespace Thermo4PFM
{

class CALPHADFreeEnergyFunctionsBinary : public CALPHADFreeEnergyFunctions
{
public:
    CALPHADFreeEnergyFunctionsBinary(boost::property_tree::ptree& input_db,
        boost::optional<boost::property_tree::ptree&> newton_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        const bool with_third_phase = false);

    ~CALPHADFreeEnergyFunctionsBinary() { delete solver_; };

    virtual double computeFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, const bool gp = false);
    virtual void computeDerivFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, double*);
    virtual void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi,
        std::vector<double>& d2fdc2);

    virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
        const PhaseIndex pi1, double* ceq, const int maxits = 20,
        const bool verbose = false);

    void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.)
    {
        std::ofstream os1("FlC0vsT.dat", std::ios::out);
        os1 << "#Species 0, Phase L" << std::endl;
        g_species_phaseL_[0].plotFofT(os1, T0, T1);

        std::ofstream os2("FlC1vsT.dat", std::ios::out);
        os2 << "#Species 1, Phase L" << std::endl;
        g_species_phaseL_[1].plotFofT(os2, T0, T1);

        std::ofstream os3("FsC0vsT.dat", std::ios::out);
        os3 << "#Species 0, Phase A" << std::endl;
        g_species_phaseA_[0].plotFofT(os3, T0, T1);

        std::ofstream os4("FsC1vsT.dat", std::ios::out);
        os4 << "#Species 1, Phase A" << std::endl;
        g_species_phaseA_[1].plotFofT(os4, T0, T1);
    }

    int computePhaseConcentrations(const double temperature, const double* conc,
        const double phi, const double eta, double* x);
    void energyVsPhiAndC(const double temperature, const double* const ceq,
        const bool found_ceq, const double phi_well_scale,
        const int npts_phi = 51,
        const int npts_c   = 50); // # of compositions to use (>1)
    void printEnergyVsComposition(
        const double temperature, const int npts = 100);
    double fchem(const double phi, const double eta, const double* const conc,
        const double temperature);
    void printEnergyVsPhiHeader(const double temperature, const int nphi,
        const int nc, const double cmin, const double cmax, const double slopec,
        std::ostream& os) const;
    void printEnergyVsPhi(const double* const conc, const double temperature,
        const double phi_well_scale, const int npts, const double slopec,
        std::ostream& os);
    void printEnergyVsEta(const double* const conc, const double temperature,
        const double eta_well_scale, const int npts, const double slopec,
        std::ostream& os);

protected:
    CALPHADConcentrationSolverBinary* solver_;

    double ceq_l_;
    double ceq_a_;
    double ceq_b_;

    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    bool with_third_phase_;

    void readNewtonparameters(boost::property_tree::ptree& newton_db);

    void computeParametersForSolvers(const double temperature, double* Lmix_L,
        double* Lmix_A, double* Lmix_B, double* fA, double* fB,
        const PhaseIndex* pis, const int nphases);

private:
    std::string fenergy_diag_filename_;

    // size 2 for species 0 and 1
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseL_[2];
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseA_[2];
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseB_[2];

    // size 4 for L0, L1, L2, L3,
    // can contain up to 3 coefficients a,b,c for a+b*T,
    // possibly +c*T*ln(T) if compiled with -DLMIX_WTLOGT
    double LmixPhaseL_[4][MAX_POL_T_INDEX];
    double LmixPhaseA_[4][MAX_POL_T_INDEX];
    double LmixPhaseB_[4][MAX_POL_T_INDEX];

    void readParameters(boost::property_tree::ptree& calphad_db);

    void setupSolver(boost::optional<boost::property_tree::ptree&> newton_db);

    // energy of species "is" in phase L,A,B
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
            case PhaseIndex::phaseB:
                return LmixPhaseB_[index][0]
                       + LmixPhaseB_[index][1] * temperature
#ifdef LMIX_WTLOGT
                       + LmixPhaseB_[index][2] * temperature * log(temperature)
#endif
                    ;
            default:
                std::cout << "CALPHADFreeEnergyStrategy::lmixPhase(), "
                             "undefined phase"
                          << "!!!" << std::endl;
                abort();
                return 0.;
        }
    }

    void computePhasesFreeEnergies(const double temperature, const double hphi,
        const double heta, const double conc, double& fl, double& fa,
        double& fb);
};
}
#endif
