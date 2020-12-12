#ifndef included_CALPHADFreeEnergyFunctions
#define included_CALPHADFreeEnergyFunctions

#include "FreeEnergyFunctions.h"
#include "Phases.h"

namespace Thermo4PFM
{

class CALPHADFreeEnergyFunctions : public FreeEnergyFunctions
{
public:
    CALPHADFreeEnergyFunctions(){};

    virtual ~CALPHADFreeEnergyFunctions(){};

    virtual int computePhaseConcentrations(const double temperature,
        const double* conc, const double phi, double* x)
        = 0;

    virtual void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2)
        = 0;
};
}

#endif
