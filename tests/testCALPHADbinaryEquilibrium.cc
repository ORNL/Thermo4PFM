#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD binary equilibrium", "[binary equilibrium]")
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    double temperature = 1423.;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphadAuNi.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // choose pair of phases: phaseL, phaseA, phaseB
    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;

    // initial guesses
    double init_guess[2] = { 0.2, 0.1 };

    double lceq[2] = { init_guess[0], init_guess[1] };

    // compute equilibrium concentrations in each phase
    bool found_ceq = cafe.computeCeqT(temperature, &lceq[0]);
    if (lceq[0] > 1.) found_ceq = false;
    if (lceq[0] < 0.) found_ceq = false;
    if (lceq[1] > 1.) found_ceq = false;
    if (lceq[1] < 0.) found_ceq = false;

    std::cout << "Temperature = " << temperature << std::endl;
    if (found_ceq)
    {
        std::cout << "Found equilibrium concentrations: " << lceq[0] << " and "
                  << lceq[1] << "..." << std::endl;
    }
    REQUIRE(found_ceq);

    double expected_result[2] = { 0.283171, 0.108235 };

    CHECK(lceq[0] == Approx(expected_result[0]).margin(1.e-6));

    CHECK(lceq[1] == Approx(expected_result[1]).margin(1.e-6));
}
