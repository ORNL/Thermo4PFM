#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD binary three phase, two sublattice KKS",
    "[binary three phase, two sublattice kks]")
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
            "../thermodynamic_data/calphadAlCuLFccTheta.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    pt::ptree newton_db;
    newton_db.put("alpha", 0.5);
    newton_db.put("max_its", 100);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.8;
    double c_init1 = 0.6;
    double c_init2 = 0.67;

    double sol[3] = { c_init0, c_init1, c_init2 };

    // compute concentrations satisfying KKS equations
    // I think the system at 820K is best behaved in the conc = 0.67-0.8 range
    double conc = 0.7;

    // This phi doesn't quite give Sum(h(phi)) = 1, but it is close
    double phi[3] = { 0.3, 0.45, 0.5 };

    cafe.computePhaseConcentrations(temperature, &conc, phi, sol);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi[0] << " "
              << phi[1] << " " << phi[2] << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cA = " << sol[1] << std::endl;
    std::cout << "   cB = " << sol[2] << std::endl;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;
    const Thermo4PFM::PhaseIndex pi2 = Thermo4PFM::PhaseIndex::phaseB;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS1;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS1);
    std::cout << "   dfS1/dcS1 = " << derivS1 << std::endl;

    REQUIRE(derivS1 == Approx(derivL).margin(1.e-5));

    double derivS2;
    cafe.computeDerivFreeEnergy(temperature, &sol[2], pi2, &derivS2);
    std::cout << "   dfS2/dcS2 = " << derivS2 << std::endl;

    REQUIRE(derivS2 == Approx(derivL).margin(1.e-5));
}

TEST_CASE("CALPHAD binary three phase, two sublattice KKS #2",
    "[binary three phase, two sublattice kks #2]")
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
            "../thermodynamic_data/calphadAlCuLFccTheta.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    pt::ptree newton_db;
    newton_db.put("alpha", 1.0);
    newton_db.put("max_its", 20000);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.79267;
    double c_init1 = 0.79267;
    double c_init2 = 0.79267;

    double sol[3] = { c_init0, c_init1, c_init2 };

    // compute concentrations satisfying KKS equations
    // I think the system at 820K is best behaved in the conc = 0.67-0.8 range
    double conc = 0.79267;

    // This phi doesn't quite give Sum(h(phi)) = 1, but it is close
    double phi[3] = { 0.938956, 5.88418e-15, 0.0610444 };

    cafe.computePhaseConcentrations(temperature, &conc, phi, sol);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature = " << temperature << std::endl;
    std::cout << "Result for c = " << conc << " and phi = " << phi[0] << " "
              << phi[1] << " " << phi[2] << std::endl;
    std::cout << "   cL = " << sol[0] << std::endl;
    std::cout << "   cA = " << sol[1] << std::endl;
    std::cout << "   cB = " << sol[2] << std::endl;

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;
    const Thermo4PFM::PhaseIndex pi2 = Thermo4PFM::PhaseIndex::phaseB;

    std::cout << "Verification:" << std::endl;

    double derivL;
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
    std::cout << "   dfL/dcL = " << derivL << std::endl;

    double derivS1;
    cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS1);
    std::cout << "   dfS1/dcS1 = " << derivS1 << std::endl;

    // REQUIRE(derivS1 == Approx(derivL).margin(1.e-5));

    double derivS2;
    cafe.computeDerivFreeEnergy(temperature, &sol[2], pi2, &derivS2);
    std::cout << "   dfS2/dcS2 = " << derivS2 << std::endl;

    REQUIRE(derivS2 == Approx(derivL).margin(1.e-5));
}
