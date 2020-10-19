#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD binary KKS", "[binary kks]")
{
    EnergyInterpolationType energy_interp_func_type
        = EnergyInterpolationType::PBG;
    ConcInterpolationType conc_interp_func_type = ConcInterpolationType::PBG;

    double temperature = 1450.;

    std::string conc_avg_func_type = "a";

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

    CALPHADFreeEnergyFunctionsBinary cafe(calphad_db, newton_db,
        energy_interp_func_type, conc_interp_func_type,
        false); // no 3rd phase

    // initial guesses
    double c_init0 = 0.5;
    double c_init1 = 0.5;

    double sol[2] = { c_init0, c_init1 };

    // compute concentrations satisfying KKS equations
    double conc = 0.3;
    double phi  = 0.5;
    cafe.computePhaseConcentrations(temperature, &conc, phi, 0., &sol[0]);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cS = " << sol[1] << std::endl;

    const PhaseIndex pi0 = PhaseIndex::phaseL;
    const PhaseIndex pi1 = PhaseIndex::phaseA;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS);
    std::cout << "   dfS/dcS = " << derivS << std::endl;

    REQUIRE(derivS == Approx(derivL).epsilon(1.e-5));
}
