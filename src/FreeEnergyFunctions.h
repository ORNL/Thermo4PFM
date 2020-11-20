#ifndef included_FreeEnergyFunctions
#define included_FreeEnergyFunctions

#include "Phases.h"

#include <string>
#include <vector>

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
        const double temperature, const int npts = 100)
        = 0;

    virtual void computeSecondDerivativeFreeEnergy(const double temp,
        const double* const conc, const PhaseIndex pi,
        std::vector<double>& d2fdc2)
        = 0;

    virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
        const PhaseIndex pi1, double* ceq, const int maxits,
        const bool verbose = false)
    {
        (void)temperature;
        (void)pi0;
        (void)pi1;
        (void)ceq;
        (void)maxits;
        (void)verbose;

        return false;
    };
};
}

#endif
