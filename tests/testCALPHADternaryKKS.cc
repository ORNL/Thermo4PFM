#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsTernary.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("KKS ternary", "[KKS ternary]")
{
    EnergyInterpolationType energy_interp_func_type
        = EnergyInterpolationType::PBG;
    ConcInterpolationType conc_interp_func_type = ConcInterpolationType::PBG;

    double temperature = 2923.;

    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphadMoNbTa.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    boost::optional<pt::ptree&> newton_db;

    CALPHADFreeEnergyFunctionsTernary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    // Mo energy at melting point
    {
        double temperature0 = 2896.;
        double conc[2]      = { 1., 0. };

        double el
            = cafe.computeFreeEnergy(temperature0, conc, PhaseIndex::phaseL);
        std::cout << "Mo, Energy in liquid phase: " << el << std::endl;
        double es
            = cafe.computeFreeEnergy(temperature0, conc, PhaseIndex::phaseA);
        std::cout << "Mo, Energy in solid phase: " << es << std::endl;

        REQUIRE(el == Approx(es).margin(1.e-1));
    }

    // Ta energy at melting point
    {
        double temperature0 = 3290.;
        double conc[2]      = { 0., 0. };

        double el
            = cafe.computeFreeEnergy(temperature0, conc, PhaseIndex::phaseL);
        std::cout << "Ta, Energy in liquid phase: " << el << std::endl;
        double es
            = cafe.computeFreeEnergy(temperature0, conc, PhaseIndex::phaseA);
        std::cout << "Ta, Energy in solid phase: " << es << std::endl;

        REQUIRE(el == Approx(es).margin(1.e-1));
    }

    // initial guesses
    double sol[4] = { 0.33, 0.38, 0.32, 0.33 };

    // compute concentrations satisfying KKS equations
    double conc[2] = { 0.33, 0.33 };
    double phi     = 0.5;
    cafe.computePhaseConcentrations(temperature, &conc[0], phi, 0., &sol[0]);

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Temperature: " << temperature << std::endl;
    std::cout << "Result for c = " << conc[0] << "," << conc[1]
              << " and phi = " << phi << std::endl;
    std::cout << "   cL = " << sol[0] << ", " << sol[1] << std::endl;
    std::cout << "   cS = " << sol[2] << ", " << sol[3] << std::endl;

    const PhaseIndex pi0 = PhaseIndex::phaseL;
    const PhaseIndex pi1 = PhaseIndex::phaseA;

    std::cout << "Verification:" << std::endl;

    double derivL[2];
    cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL[0]);

    double derivS[2];
    cafe.computeDerivFreeEnergy(temperature, &sol[2], pi1, &derivS[0]);

    std::cout << "   dfL/dcL0 = " << derivL[0] << std::endl;
    std::cout << "   dfS/dcS0 = " << derivS[0] << std::endl;
    std::cout << "   dfL/dcL1 = " << derivL[1] << std::endl;
    std::cout << "   dfS/dcS1 = " << derivS[1] << std::endl;

    const double tol = 1.e-5;
    CHECK(derivS[0] == Approx(derivL[0]).margin(tol));
    CHECK(derivS[1] == Approx(derivL[1]).margin(tol));
}
