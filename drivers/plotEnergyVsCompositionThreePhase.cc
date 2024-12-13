#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFunctions.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    std::string databasename(argv[1]);
    double temperature = atof(argv[2]);
    const int npts     = atof(argv[3]);

    double cmin = 0.;
    double cmax = 1.;
    if (argc > 4)
    {
        cmin = atof(argv[3]);
        cmax = atof(argv[4]);
    }

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    // pick a very small temperature to have a very small contribution of
    // pure mixing term, but
    // cannot be 0 as it would trigger nans in f(T) evaluations
    assert(temperature > 0.);

    std::cout << " Read CALPHAD database " << databasename << std::endl;
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

    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(2);
    ss << temperature;
    std::ofstream os("FvsC" + ss.str() + ".csv", std::ios::out);
    std::cout << " Compute energies..." << std::endl;

    if (Thermo4PFM::checkSublattice(calphad_db))
    {
        Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(calphad_db,
            newton_db, energy_interp_func_type, conc_interp_func_type);

        cafe.printEnergyVsComposition(temperature, os, cmin, cmax, npts);
    }
    else
    {
        Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase cafe(calphad_db,
            newton_db, energy_interp_func_type, conc_interp_func_type);

        cafe.printEnergyVsComposition(temperature, os, cmin, cmax, npts);
    }
}
