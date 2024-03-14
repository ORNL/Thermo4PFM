#ifndef Thermo4PFM_included_CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB
#define Thermo4PFM_included_CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB

#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"

namespace Thermo4PFM
{

class CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB
    : public CALPHADFreeEnergyFunctionsBinaryThreePhase
{
public:
    CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB(const double cB,
        boost::property_tree::ptree& input_db,
        boost::optional<boost::property_tree::ptree&> newton_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type);

    int computePhaseConcentrations(const double temperature, const double* conc,
        const double* const phi, double* x) override;

    bool computeCeqT(const double temperature, double* ceq,
        const int maxits = 20, const bool verbose = false);

private:
    // stochiometric composition of phaseB
    const double cB_;
};
}
#endif
