#ifndef included_CALPHADFreeEnergyFunctions
#define included_CALPHADFreeEnergyFunctions

#include "FreeEnergyFunctions.h"
#include "Phases.h"

#include <vector>

namespace Thermo4PFM
{

class CALPHADFreeEnergyFunctions : public FreeEnergyFunctions
{
public:
    CALPHADFreeEnergyFunctions(){};

    virtual ~CALPHADFreeEnergyFunctions(){};

    virtual int computePhaseConcentrations(const double temperature,
        const double* conc, const double phi, const double eta, double* x)
        = 0;

    virtual void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi,
        std::vector<double>& d2fdc2)
        = 0;
};
}

#endif
