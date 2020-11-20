#ifndef included_CALPHADSpeciesPhaseGibbsEnergy
#define included_CALPHADSpeciesPhaseGibbsEnergy

#include "CALPHADSpeciesPhaseGibbsEnergyExpansion.h"

#include <boost/property_tree/ptree.hpp>

#include <string>
#include <vector>

namespace Thermo4PFM
{

class CALPHADSpeciesPhaseGibbsEnergy
{
private:
    std::string name_;
    std::vector<double> tc_;

    std::vector<CALPHADSpeciesPhaseGibbsEnergyExpansion> expansion_;

public:
    CALPHADSpeciesPhaseGibbsEnergy(){};

    std::string name() const { return name_; }

    void initialize(const std::string& name, boost::property_tree::ptree& db);

    double fenergy(const double T); // expect T in Kelvin
    void plotFofT(
        std::ostream& os, const double T0 = 300., const double T1 = 3000.);
};
}
#endif
