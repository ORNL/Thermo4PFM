#ifndef Thermo4PFM_included_CALPHADSpeciesPhaseGibbsEnergyExpansion
#define Thermo4PFM_included_CALPHADSpeciesPhaseGibbsEnergyExpansion

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif

template <typename ScalarType>
class CALPHADSpeciesPhaseGibbsEnergyExpansion
{
public:
    CALPHADSpeciesPhaseGibbsEnergyExpansion(){};
    void init(const ScalarType a, const ScalarType b, const ScalarType c,
        const ScalarType d2, const ScalarType d3, const ScalarType d4,
        const ScalarType d7, const ScalarType dm1, const ScalarType dm9);

    double value(const double temperature) const;

private:
    /*
     * Expansion coefficient for energy of species as a function of temperature
     */
    ScalarType a_;
    ScalarType b_;
    ScalarType c_;
    ScalarType d2_;
    ScalarType d3_;
    ScalarType d4_;
    ScalarType d7_;
    ScalarType dm1_;
    ScalarType dm9_;
};

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif
}
#endif
