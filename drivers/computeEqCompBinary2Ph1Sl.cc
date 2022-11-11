#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        std::cerr << "ERROR: program needs 2 argument (database + temperature "
                     "+cl_init +cs_init )"
                  << std::endl;
        return 1;
    }

    std::string databasename(argv[1]);
    double temperature = atof(argv[2]);
    double cl          = atof(argv[3]);
    double cs          = atof(argv[4]);

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    std::cout << "Temperature: " << temperature << std::endl;

    std::cout << "Read CALPHAD database..." << std::endl;
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

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    double conc[2] = { cl, cs };
    cafe.computeCeqT(temperature, conc, 25, true);

    std::cout << "Equilibrium compositions: cl = " << conc[0]
              << ", cs=" << conc[1] << std::endl;
}
