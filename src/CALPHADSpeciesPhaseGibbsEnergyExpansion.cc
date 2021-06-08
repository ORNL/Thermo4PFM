#include "CALPHADSpeciesPhaseGibbsEnergyExpansion.h"

//#include <iostream>
#include <math.h>

namespace Thermo4PFM
{

#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp declare target
#endif
template <typename ScalarType>
void CALPHADSpeciesPhaseGibbsEnergyExpansion<ScalarType>::init(
    const ScalarType a, const ScalarType b, const ScalarType c,
    const ScalarType d2, const ScalarType d3, const ScalarType d4,
    const ScalarType d7, const ScalarType dm1, const ScalarType dm9)
{
    a_   = a;
    b_   = b;
    c_   = c;
    d2_  = d2;
    d3_  = d3;
    d4_  = d4;
    d7_  = d7;
    dm1_ = dm1;
    dm9_ = dm9;
}

template <typename ScalarType>
double CALPHADSpeciesPhaseGibbsEnergyExpansion<ScalarType>::value(
    const double temperature) const
{
    const double t2 = temperature * temperature;
    const double t4 = t2 * t2;
    // std::cout<<"a="<<a_<<std::endl;
    // std::cout<<"b="<<b_<<std::endl;
    // std::cout<<"c="<<c_<<std::endl;
    // std::cout<<"d2="<<d2_<<std::endl;
    // std::cout<<"d3="<<d3_<<std::endl;
    // std::cout<<"d4="<<d4_<<std::endl;
    // std::cout<<"d7="<<d7_<<std::endl;
    // std::cout<<"dm1="<<dm1_<<std::endl;
    return a_ + b_ * temperature + c_ * temperature * log(temperature)
           + d2_ * t2 + d3_ * t2 * temperature + d4_ * t4
           + d7_ * t4 * t2 * temperature + dm1_ / temperature
           + dm9_ / (t4 * t4 * temperature);
}
#ifdef HAVE_OPENMP_OFFLOAD
#pragma omp end declare target
#endif

template class CALPHADSpeciesPhaseGibbsEnergyExpansion<float>;
template class CALPHADSpeciesPhaseGibbsEnergyExpansion<double>;
}
