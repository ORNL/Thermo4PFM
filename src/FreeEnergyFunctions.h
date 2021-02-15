#ifndef included_FreeEnergyFunctions
#define included_FreeEnergyFunctions

#include "Phases.h"

#include <string>

namespace Thermo4PFM
{

class FreeEnergyFunctions
{
public:
    FreeEnergyFunctions(){};

    virtual ~FreeEnergyFunctions(){};

    virtual void energyVsPhiAndC(const double temperature,
        const double* const ceq, const bool found_ceq,
        const double phi_well_scale, const int npts_phi = 51,
        const int npts_c = 50)
        = 0;
    virtual void printEnergyVsComposition(
        const double temperature, std::ostream& os, const int npts = 100)
        = 0;

    virtual void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi, double* d2fdc2)
        = 0;

    virtual bool computeCeqT(const double temperature, double* ceq,
        const int maxits, const bool verbose = false)
    {
        (void)temperature;
        (void)ceq;
        (void)maxits;
        (void)verbose;

        return false;
    };
};
}

#endif
