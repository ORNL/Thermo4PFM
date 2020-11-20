#ifndef included_CALPHADSpeciesPhaseGibbsEnergyExpansion
#define included_CALPHADSpeciesPhaseGibbsEnergyExpansion

#include <string>

namespace Thermo4PFM
{

class CALPHADSpeciesPhaseGibbsEnergyExpansion
{
public:
    CALPHADSpeciesPhaseGibbsEnergyExpansion(const double a, const double b,
        const double c, const double d2, const double d3, const double d4,
        const double d7, const double dm1, const double dm9);

    double value(const double temperature) const;

private:
    /*
     * Expansion coefficient for energy of species as a function of temperature
     */
    const double a_;
    const double b_;
    const double c_;
    const double d2_;
    const double d3_;
    const double d4_;
    const double d7_;
    const double dm1_;
    const double dm9_;
};
}
#endif
