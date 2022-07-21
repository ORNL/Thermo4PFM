#ifndef Thermo4PFM_included_CALPHADSpeciesPhaseGibbsEnergy
#define Thermo4PFM_included_CALPHADSpeciesPhaseGibbsEnergy

#include "CALPHADSpeciesPhaseGibbsEnergyExpansion.h"

#include <boost/property_tree/ptree.hpp>

#include <cassert>
#include <string>

#define MAXNINTERVALS 5

namespace Thermo4PFM
{
typedef double coeffsdatatype;

class CALPHADSpeciesPhaseGibbsEnergy
{
private:
    char* name_;
    int nintervals_;

    double tc_[MAXNINTERVALS + 1];
    CALPHADSpeciesPhaseGibbsEnergyExpansion<coeffsdatatype>
        expansion_[MAXNINTERVALS];

public:
    CALPHADSpeciesPhaseGibbsEnergy() {}

    ~CALPHADSpeciesPhaseGibbsEnergy() { delete[] name_; }

    int nintervals() const { return nintervals_; }

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
