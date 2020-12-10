#ifndef included_KKSFreeEnergyFunctionDiluteBinary
#define included_KKSFreeEnergyFunctionDiluteBinary

#include "FreeEnergyFunctions.h"
#include "InterpolationType.h"
#include "KKSdiluteBinaryConcentrationSolver.h"
#include "Phases.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace Thermo4PFM
{

class KKSFreeEnergyFunctionDiluteBinary : public FreeEnergyFunctions
{
public:
    KKSFreeEnergyFunctionDiluteBinary(boost::property_tree::ptree& conc_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    ~KKSFreeEnergyFunctionDiluteBinary() { delete solver_; };

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

    void preRunDiagnostics(const double T0 = 300., const double T1 = 3000.) {}

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
    void printEnergyVsEta(const double* const conc, const double temperature,
        const int npts, const double slopec, std::ostream& os);

private:
    KKSdiluteBinaryConcentrationSolver* solver_;

    double ceq_l_;
    double ceq_a_;

    EnergyInterpolationType energy_interp_func_type_;
    ConcInterpolationType conc_interp_func_type_;

    void readNewtonparameters(boost::property_tree::ptree& newton_db);

    void setupFB(const double temperature);

    std::string fenergy_diag_filename_;

    double fA_;
    double fB_;

    double Tm_;
    double me_;
    double ke_;

    void readParameters(boost::property_tree::ptree& conc_db);

    void setupSolver(boost::optional<boost::property_tree::ptree&> newton_db);

    void computePhasesFreeEnergies(const double temperature, const double hphi,
        const double conc, double& fl, double& fa);
};
}
#endif
