#ifndef included_CALPHADSpeciesPhaseGibbsEnergy
#define included_CALPHADSpeciesPhaseGibbsEnergy

#include "CALPHADSpeciesPhaseGibbsEnergyExpansion.h"

#include <boost/property_tree/ptree.hpp>

#include <cassert>
#include <string>

namespace Thermo4PFM
{
typedef double coeffsdatatype;

class CALPHADSpeciesPhaseGibbsEnergy
{
private:
    char* name_;
    int nintervals_;
    double* tc_;
    CALPHADSpeciesPhaseGibbsEnergyExpansion<coeffsdatatype>* expansion_;

public:
    CALPHADSpeciesPhaseGibbsEnergy()
    {
        tc_        = nullptr;
        expansion_ = nullptr;
    }

    ~CALPHADSpeciesPhaseGibbsEnergy()
    {
        assert(expansion_ != nullptr);
        assert(tc_ != nullptr);
        delete[] name_;
        delete[] expansion_;
        delete[] tc_;
    }

    std::string name() const
    {
        std::string name(name_);
        return name;
    }

    void initialize(const std::string& name, boost::property_tree::ptree& db);

    double fenergy(const double T); // expect T in Kelvin
    void plotFofT(
        std::ostream& os, const double T0 = 300., const double T1 = 3000.);
};
}
#endif
