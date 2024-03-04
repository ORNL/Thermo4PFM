#include "CALPHADFreeEnergyFunctionsBinary.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Requires two arguments: database name + temperature"
                  << std::endl;
        return 1;
    }

    std::string databasename(argv[1]);
    double temperature = atof(argv[2]);

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json(databasename, calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    const int npts = 101;
    std::ofstream os("FvsC.dat", std::ios::out);

    std::cout << " Compute energies..." << std::endl;
    cafe.printEnergyVsComposition(temperature, os, npts);
}
