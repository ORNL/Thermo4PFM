#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsTernary.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <array>
#include <fstream>
#include <string>

namespace pt = boost::property_tree;

TEST_CASE("CALPHADternaryEquilibrium", "[CALPHADternaryEquilibrium]")
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

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

    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    double temperature[2] = { 2923., 3073. };

    // initial guesses
    double init_guess[2][5]
        = { { 0.33, 0.38, 0.32, 0.33, 0.8 }, { 0.13, 0.21, 0.12, 0.15, 0.8 } };

    double nominalc[2][2] = { { 0.33, 0.33 }, { 0.16, 0.225 } };

    double expected_cl[2][2] = { { 0.344909, 0.380528 }, { 0.164344, 0.2321 } };
    double expected_cs[2][2]
        = { { 0.329373, 0.327874 }, { 0.13675, 0.186999 } };
    double expected_fs[2] = { 0.959621, 0.157429 };

    int maxits = 20;

    for (int itest = 0; itest < 2; itest++)
    {
        std::cout << "----------------------------------------" << std::endl;
        double lceq[5] = { init_guess[itest][0], init_guess[itest][1], // liquid
            init_guess[itest][2], init_guess[itest][3], // solid
            init_guess[itest][4] };

        bool found_ceq = cafe.computeCeqT(temperature[itest],
            nominalc[itest][0], nominalc[itest][1], &lceq[0], maxits);
        if (lceq[0] != lceq[0]) found_ceq = false;
        if (lceq[0] > 1.) found_ceq = false;
        if (lceq[0] < 0.) found_ceq = false;
        if (lceq[1] > 1.) found_ceq = false;
        if (lceq[1] < 0.) found_ceq = false;

        std::cout << "Temperature = " << temperature[itest] << std::endl;
        std::cout << "Nominal composition " << nominalc[itest][0] << ","
                  << nominalc[itest][1] << std::endl;
        REQUIRE(found_ceq);
        std::cout << "Found equilibrium concentrations: " << std::endl;
        std::cout << "Liquid: " << lceq[0] << "," << lceq[1] << std::endl;
        std::cout << "Solid:  " << lceq[2] << "," << lceq[3] << std::endl;
        std::cout << "Solid fraction: " << lceq[4] << std::endl;

        // test values
        const double tol = 1.e-6;
        CHECK(lceq[0] == Approx(expected_cl[itest][0]).margin(tol));
        CHECK(lceq[1] == Approx(expected_cl[itest][1]).margin(tol));
        CHECK(lceq[2] == Approx(expected_cs[itest][0]).margin(tol));
        CHECK(lceq[3] == Approx(expected_cs[itest][1]).margin(tol));
        CHECK(lceq[4] == Approx(expected_fs[itest]).margin(tol));
    }
}
