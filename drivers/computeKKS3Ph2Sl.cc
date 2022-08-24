#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    if (argc < 7)
    {
        std::cerr << "ERROR: program needs 6 argument (database + temperature "
                     "+ phiL,phiA,phiB + c)"
                  << std::endl;
        return 1;
    }

    std::string databasename(argv[1]);
    double temperature = atof(argv[2]);
    double phiL        = atof(argv[3]);
    double phiA        = atof(argv[4]);
    double phiB        = atof(argv[5]);
    double c           = atof(argv[6]);

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    std::cout << "Temperature: " << temperature << std::endl;
    double nominalc[2] = { 0.6, 0.035 };

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

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    double phi[3]  = { phiL, phiA, phiB };
    double conc[3] = { c, c, c };
    cafe.computePhaseConcentrations(temperature, &c, phi, conc);
    std::cout << "phi = (" << phiL << "," << phiA << "," << phiB << ")"
              << ", c=" << c << std::endl;
    std::cout << "KKS solution: cl = " << conc[0] << ", ca=" << conc[1]
              << " and cb = " << conc[2] << std::endl;
}
