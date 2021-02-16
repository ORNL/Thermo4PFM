#ifndef included_CALPHADFreeEnergyFunctionsTernary
#define included_CALPHADFreeEnergyFunctionsTernary

#include "CALPHADConcSolverTernary.h"
#include "CALPHADEqConcSolverTernary.h"
#include "CALPHADFreeEnergyFunctions.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "CALPHADTieLineConcSolverTernary.h"
#include "InterpolationType.h"
#include "Phases.h"
#include "functions.h"

#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <math.h>

namespace Thermo4PFM
{

class CALPHADFreeEnergyFunctionsTernary : public CALPHADFreeEnergyFunctions
{
public:
    CALPHADFreeEnergyFunctionsTernary(boost::property_tree::ptree& input_db,
        boost::optional<boost::property_tree::ptree&> newton_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    ~CALPHADFreeEnergyFunctionsTernary(){};

    virtual double computeFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, const bool gp = false);
    virtual void computeDerivFreeEnergy(const double temperature,
        const double* const conc, const PhaseIndex pi, double* deriv);
    virtual void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2) override;

    virtual bool computeCeqT(const double temperature, double* ceq,
        const int maxits = 20, const bool verbose = false) override;

    /// Compute compositions and phase fractions ate ends of tie line
    /// passing through nominal composition (c0,c1)
    virtual bool computeTieLine(const double temperature, const double c0,
        const double c1, double* ceq, const int maxits = 20,
        const bool verbose = false);

    void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.);

    int computePhaseConcentrations(const double temperature,
        const double* const conc, const double phi, double* x) override;
    void energyVsPhiAndC(const double temperature, const double* const ceq,
        const bool found_ceq, const double phi_well_scale,
        const int npts_phi = 51,
        const int npts_c   = 50) override; // number of compositions to use (>1)
    void printEnergyVsComposition(const double temperature, std::ostream& os,
        const int npts = 100) override;
    double fchem(
        const double phi, const double* const conc, const double temperature);
    void printEnergyVsPhiHeader(const double temperature, const int nphi,
        const int nc0, const int nc1, const double c0min, const double c0max,
        const double c1min, const double c1max, std::ostream& os) const;
    void printEnergyVsPhi(const double* const conc, const double temperature,
        const double phi_well_scale, const int npts, std::ostream& os);

private:
    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    void readNewtonparameters(boost::property_tree::ptree& newton_db);

    void computeTdependentParameters(const double temperature, double* L_AB_L,
        double* L_AC_L, double* L_BC_L, double* L_ABC_L, double* L_AB_S,
        double* L_AC_S, double* L_BC_S, double* L_ABC_S, double* fA, double* fB,
        double* fC);

    std::string fenergy_diag_filename_;

    double newton_tol_;
    double newton_alpha_;
    int newton_maxits_;
    bool newton_verbose_;

    // size 3 for species A, B, C
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseL_[3];
    CALPHADSpeciesPhaseGibbsEnergy g_species_phaseA_[3];

    // size 4 for L0, L1, L2, L3, with 2 coefficient for linear expansion in T
    // a+b*T
    double LmixABPhaseL_[4][2];
    double LmixABPhaseA_[4][2];

    double LmixACPhaseL_[4][2];
    double LmixACPhaseA_[4][2];

    double LmixBCPhaseL_[4][2];
    double LmixBCPhaseA_[4][2];

    double LmixABCPhaseL_[3][2];
    double LmixABCPhaseA_[3][2];

    double (*fun_ptr_arr_[3])(const double){ linear_interp_func,
        pbg_interp_func, harmonic_interp_func };

    void readParameters(boost::property_tree::ptree& calphad_db);

    // energy of species "is" in phase L,A,B
    double getFenergyPhaseL(const short is, const double temperature)
    {
        return g_species_phaseL_[is].fenergy(temperature);
    }
    double getFenergyPhaseA(const short is, const double temperature)
    {
        return g_species_phaseA_[is].fenergy(temperature);
    }

    double lmix0ABPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix0ABPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix0ABPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix1ABPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix1ABPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix1ABPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix2ABPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix2ABPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix2ABPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix3ABPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix3ABPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix3ABPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix0ABPhaseL(const double temperature)
    {
        return LmixABPhaseL_[0][0] + LmixABPhaseL_[0][1] * temperature;
    }

    double lmix1ABPhaseL(const double temperature)
    {
        return LmixABPhaseL_[1][0] + LmixABPhaseL_[1][1] * temperature;
    }

    double lmix2ABPhaseL(const double temperature)
    {
        return LmixABPhaseL_[2][0] + LmixABPhaseL_[2][1] * temperature;
    }

    double lmix3ABPhaseL(const double temperature)
    {
        return LmixABPhaseL_[3][0] + LmixABPhaseL_[3][1] * temperature;
    }

    double lmix0ABPhaseA(const double temperature)
    {
        return LmixABPhaseA_[0][0] + LmixABPhaseA_[0][1] * temperature;
    }

    double lmix1ABPhaseA(const double temperature)
    {
        return LmixABPhaseA_[1][0] + LmixABPhaseA_[1][1] * temperature;
    }

    double lmix2ABPhaseA(const double temperature)
    {
        return LmixABPhaseA_[2][0] + LmixABPhaseA_[2][1] * temperature;
    }

    double lmix3ABPhaseA(const double temperature)
    {
        return LmixABPhaseA_[3][0] + LmixABPhaseA_[3][1] * temperature;
    }

    double lmix0ACPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix0ACPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix0ACPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix1ACPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix1ACPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix1ACPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix2ACPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix2ACPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix2ACPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix3ACPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix3ACPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix3ACPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix0ACPhaseL(const double temperature)
    {
        return LmixACPhaseL_[0][0] + LmixACPhaseL_[0][1] * temperature;
    }

    double lmix1ACPhaseL(const double temperature)
    {
        return LmixACPhaseL_[1][0] + LmixACPhaseL_[1][1] * temperature;
    }

    double lmix2ACPhaseL(const double temperature)
    {
        return LmixACPhaseL_[2][0] + LmixACPhaseL_[2][1] * temperature;
    }

    double lmix3ACPhaseL(const double temperature)
    {
        return LmixACPhaseL_[3][0] + LmixACPhaseL_[3][1] * temperature;
    }

    double lmix0ACPhaseA(const double temperature)
    {
        return LmixACPhaseA_[0][0] + LmixACPhaseA_[0][1] * temperature;
    }

    double lmix1ACPhaseA(const double temperature)
    {
        return LmixACPhaseA_[1][0] + LmixACPhaseA_[1][1] * temperature;
    }

    double lmix2ACPhaseA(const double temperature)
    {
        return LmixACPhaseA_[2][0] + LmixACPhaseA_[2][1] * temperature;
    }

    double lmix3ACPhaseA(const double temperature)
    {
        return LmixACPhaseA_[3][0] + LmixACPhaseA_[3][1] * temperature;
    }

    double lmix0BCPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix0BCPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix0BCPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix1BCPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix1BCPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix1BCPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix2BCPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix2BCPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix2BCPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix3BCPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix3BCPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix3BCPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix0BCPhaseL(const double temperature)
    {
        return LmixBCPhaseL_[0][0] + LmixBCPhaseL_[0][1] * temperature;
    }

    double lmix1BCPhaseL(const double temperature)
    {
        return LmixBCPhaseL_[1][0] + LmixBCPhaseL_[1][1] * temperature;
    }

    double lmix2BCPhaseL(const double temperature)
    {
        return LmixBCPhaseL_[2][0] + LmixBCPhaseL_[2][1] * temperature;
    }

    double lmix3BCPhaseL(const double temperature)
    {
        return LmixBCPhaseL_[3][0] + LmixBCPhaseL_[3][1] * temperature;
    }

    double lmix0BCPhaseA(const double temperature)
    {
        return LmixBCPhaseA_[0][0] + LmixBCPhaseA_[0][1] * temperature;
    }

    double lmix1BCPhaseA(const double temperature)
    {
        return LmixBCPhaseA_[1][0] + LmixBCPhaseA_[1][1] * temperature;
    }

    double lmix2BCPhaseA(const double temperature)
    {
        return LmixBCPhaseA_[2][0] + LmixBCPhaseA_[2][1] * temperature;
    }

    double lmix3BCPhaseA(const double temperature)
    {
        return LmixBCPhaseA_[3][0] + LmixBCPhaseA_[3][1] * temperature;
    }

    // ABC
    double lmix0ABCPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix0ABCPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix0ABCPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix1ABCPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix1ABCPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix1ABCPhaseA(temperature);
            default:
                return NAN;
        }
    }

    double lmix2ABCPhase(const PhaseIndex pi, const double temperature)
    {
        switch (pi)
        {
            case PhaseIndex::phaseL:
                return lmix2ABCPhaseL(temperature);
            case PhaseIndex::phaseA:
                return lmix2ABCPhaseA(temperature);
            default:
                return NAN;
        }
    }

    // ABC liquid
    double lmix0ABCPhaseL(const double temperature)
    {
        assert(LmixABCPhaseL_[0][0] == LmixABCPhaseL_[0][0]);

        return LmixABCPhaseL_[0][0] + LmixABCPhaseL_[0][1] * temperature;
    }

    double lmix1ABCPhaseL(const double temperature)
    {
        return LmixABCPhaseL_[1][0] + LmixABCPhaseL_[1][1] * temperature;
    }

    double lmix2ABCPhaseL(const double temperature)
    {
        return LmixABCPhaseL_[2][0] + LmixABCPhaseL_[2][1] * temperature;
    }

    // ABC solid
    double lmix0ABCPhaseA(const double temperature)
    {
        return LmixABCPhaseA_[0][0] + LmixABCPhaseA_[0][1] * temperature;
    }

    double lmix1ABCPhaseA(const double temperature)
    {
        return LmixABCPhaseA_[1][0] + LmixABCPhaseA_[1][1] * temperature;
    }

    double lmix2ABCPhaseA(const double temperature)
    {
        return LmixABCPhaseA_[2][0] + LmixABCPhaseA_[2][1] * temperature;
    }

    void computePhasesFreeEnergies(const double temperature, const double hphi,
        const double conc0, const double conc1, double& fl, double& fa);
};
}
#endif
