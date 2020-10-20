#define CATCH_CONFIG_MAIN

#include "KKSFreeEnergyFunctionDiluteBinary.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

namespace pt = boost::property_tree;

TEST_CASE("Dilute binary KKS", "[dilute binary kks]")
{
    EnergyInterpolationType energy_interp_func_type
        = EnergyInterpolationType::PBG;
    ConcInterpolationType conc_interp_func_type = ConcInterpolationType::LINEAR;

    // read dilute alloy informaton
    pt::ptree conc_db;
    try
    {
        pt::read_json("dilute_binary.json", conc_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }
    pt::write_json(std::clog, conc_db);

    KKSFreeEnergyFunctionDiluteBinary cafe(
        conc_db, energy_interp_func_type, conc_interp_func_type);

    const PhaseIndex pi0 = PhaseIndex::phaseL;
    const PhaseIndex pi1 = PhaseIndex::phaseA;

    // initial guesses
    double c_init0 = 0.5;
    double c_init1 = 0.5;
    double sol[2]  = { c_init0, c_init1 };

    // solve KKS equations
    double conc        = 0.05;
    double phi         = 0.5;
    double temperature = 1710.;
    cafe.computePhaseConcentrations(temperature, &conc, phi, 0., &sol[0]);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cS = " << sol[1] << std::endl;

    // verify concentrations satisfy equal chemical potentials
    std::cout << "Verification:" << std::endl;

    double derivL;
    double derivS;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS);
    std::cout << "   dfS/dcS = " << derivS << std::endl;

    CHECK(derivS == Approx(derivL).margin(1.e-5));
}
