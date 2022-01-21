#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    // very small temperature to have a very small contribution of
    // pure mixing term
    // Cannot be 0 as it would trigger nans in f(T) evaluations
    double temperature = 1.e-8;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphad3phases.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    const int npts = 100;
    std::ofstream os("FvsC.dat", std::ios::out);

    std::cout << " Compute energies..." << std::endl;
    cafe.printEnergyVsComposition(temperature, os, npts);
}
