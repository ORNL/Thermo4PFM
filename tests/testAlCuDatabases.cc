#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <iostream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("AlCu database check, two phase", "[AlCu database checks, two phase]")
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphadAlCuLFcc.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    pt::ptree newton_db;
    newton_db.put("alpha", 0.1);
    newton_db.put("max_its", 10000);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;

    double temperature = 0.0;
    double conc        = 1.0;
    double f0          = 0.0;
    double f1          = 0.0;

    // Check the melting temperature
    temperature = 1000.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    REQUIRE(f0 <= f1);

    temperature = 934.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    REQUIRE(f0 <= f1);

    temperature = 933.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    REQUIRE(f0 >= f1);

    temperature = 900.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    REQUIRE(f0 >= f1);

    // Check at the eutectic (even though only two phases)
    temperature = 1000.0;
    conc        = 1.0 - 0.175;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    REQUIRE(f0 <= f1);

    temperature = 822.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    REQUIRE(f0 <= f1);

    // Check in the two-phase liquid-fcc region
    temperature          = 900.0;
    double init_guess[2] = { 0.8, 0.99 };
    double lceq[2]       = { init_guess[0], init_guess[1] };
    bool found_ceq       = cafe.computeCeqT(temperature, &lceq[0], 1000, true);
    REQUIRE(found_ceq);
    double expected_result[2] = { 0.94, 0.99 };
    CHECK(lceq[0] == Approx(expected_result[0]).margin(0.01));
    CHECK(lceq[1] == Approx(expected_result[1]).margin(0.01));

    temperature = 822.0;
    lceq[0]     = 0.8;
    lceq[1]     = 0.99;
    found_ceq   = cafe.computeCeqT(temperature, &lceq[0], 1000, true);
    REQUIRE(found_ceq);
    expected_result[0] = 1.0 - 0.175;
    expected_result[1] = 1.0 - 0.025;
    CHECK(lceq[0] == Approx(expected_result[0]).margin(0.01));
    CHECK(lceq[1] == Approx(expected_result[1]).margin(0.01));

    // Check a point on the curves
    temperature = 812.0;
    conc        = 0.8743844133668341;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    CHECK(f0 == Approx(-36490.44098336194).margin(0.0001));
    CHECK(f1 == Approx(-36449.059817785725).margin(0.0001));
}

TEST_CASE(
    "AlCu database check, three phase", "[AlCu database checks, three phase]")
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

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
    newton_db.put("alpha", 0.1);
    newton_db.put("max_its", 10000);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;
    const Thermo4PFM::PhaseIndex pi2 = Thermo4PFM::PhaseIndex::phaseB;

    double temperature = 0.0;
    double conc        = 1.0;
    double f0          = 0.0;
    double f1          = 0.0;
    double f2          = 0.0;

    // Check the Al melting temperature
    temperature = 1000.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    f2          = cafe.computeFreeEnergy(temperature, &conc, pi2);
    REQUIRE(f0 <= f1);
    REQUIRE(f0 <= f2);

    temperature = 934.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    f2          = cafe.computeFreeEnergy(temperature, &conc, pi2);
    REQUIRE(f0 <= f1);
    REQUIRE(f0 <= f2);

    temperature = 933.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    f2          = cafe.computeFreeEnergy(temperature, &conc, pi2);
    REQUIRE(f0 >= f1);
    REQUIRE(f1 <= f2);

    temperature = 900.0;
    conc        = 1.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    f2          = cafe.computeFreeEnergy(temperature, &conc, pi2);
    REQUIRE(f0 >= f1);
    REQUIRE(f1 <= f2);

    // Check the liquid above the eutectic
    temperature = 1000.0;
    conc        = 1.0 - 0.175;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    f2          = cafe.computeFreeEnergy(temperature, &conc, pi2);
    REQUIRE(f0 <= f1);
    REQUIRE(f0 <= f2);

    temperature = 822.0;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    f2          = cafe.computeFreeEnergy(temperature, &conc, pi2);
    REQUIRE(f0 < f1);
    REQUIRE(f0 < f2);

    // Check a point on the curves
    temperature = 812.0;
    conc        = 0.8743844133668341;
    f0          = cafe.computeFreeEnergy(temperature, &conc, pi0);
    f1          = cafe.computeFreeEnergy(temperature, &conc, pi1);
    f2          = cafe.computeFreeEnergy(temperature, &conc, pi2);
    CHECK(f0 == Approx(-36490.44098336194).margin(0.0001));
    CHECK(f1 == Approx(-36449.059817785725).margin(0.0001));
    CHECK(f2 == Approx(-34045.5790319836).margin(0.0001));
}
