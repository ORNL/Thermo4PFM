#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
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
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    double temperature = 820.;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json(
            "../thermodynamic_data/calphadAlCuLTheta.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.5;
    double c_init1 = 0.5;

    double sol[2] = { c_init0, c_init1 };

    // compute concentrations satisfying KKS equations
    double conc = 0.7;
    double phi  = 0.5;
    cafe.computePhaseConcentrations(temperature, &conc, &phi, &sol[0]);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cS = " << sol[1] << std::endl;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS);
    std::cout << "   dfS/dcS = " << derivS << std::endl;

    REQUIRE(derivS == Approx(derivL).margin(1.e-5));
}
