#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "Phases.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "ERROR: program needs 2 argument (database + temperature "
                  << std::endl;
        return 1;
    }

    std::string databasename(argv[1]);
    double temperature = atof(argv[2]);

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

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    double conc[2];
    cafe.computeCeqT(temperature, conc);

    std::cout << "Equilibrium compositions: cl = " << conc[0]
              << ", cs=" << conc[1] << std::endl;

    double deriv;
    cafe.computeDerivFreeEnergy(
        temperature, conc, Thermo4PFM::PhaseIndex::phaseL, &deriv);

    std::cout << "Chemical potential: mu = " << deriv << std::endl;

    double d2fdc2;
    cafe.computeSecondDerivativeFreeEnergy(
        temperature, conc, Thermo4PFM::PhaseIndex::phaseL, &d2fdc2);
    std::cout << "Second derivative for L phase: " << d2fdc2 << std::endl;

    cafe.computeSecondDerivativeFreeEnergy(
        temperature, conc + 1, Thermo4PFM::PhaseIndex::phaseA, &d2fdc2);
    std::cout << "Second derivative for A phase: " << d2fdc2 << std::endl;
}
