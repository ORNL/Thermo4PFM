#ifndef included_CALPHADSpeciesPhaseGibbsEnergyExpansion
#define included_CALPHADSpeciesPhaseGibbsEnergyExpansion

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

class CALPHADSpeciesPhaseGibbsEnergyExpansion
{
public:
    CALPHADSpeciesPhaseGibbsEnergyExpansion(){};
    void init(const double a, const double b, const double c, const double d2,
        const double d3, const double d4, const double d7, const double dm1,
        const double dm9);

    double value(const double temperature) const;

private:
    /*
     * Expansion coefficient for energy of species as a function of temperature
     */
    double a_;
    double b_;
    double c_;
    double d2_;
    double d3_;
    double d4_;
    double d7_;
    double dm1_;
    double dm9_;
};

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
#endif
