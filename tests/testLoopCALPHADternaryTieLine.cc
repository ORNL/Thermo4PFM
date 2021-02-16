#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

#include <omp.h>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD ternary equilibrium", "[ternary equilibrium]")
{
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    const double Tmin     = 2990;
    const double Tmax     = 3020.;
    const int nTintervals = 10;
    const double deltaT   = (Tmax - Tmin) / (double)nTintervals;

    std::cout << " Read CALPHAD database..." << std::endl;
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

    // choose pair of phases: phaseL, phaseA, phaseB
    const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;
    const Thermo4PFM::PhaseIndex pi1 = Thermo4PFM::PhaseIndex::phaseA;

    // initial guesses
    const double init_guess[5] = { 0.33, 0.38, 0.32, 0.33, 0.8 };

    // nominal composition
    double nominalc[2] = { 0.2, 0.3 };

    std::vector<double> cel(2 * (nTintervals + 1));
    std::vector<double> ces(2 * (nTintervals + 1));

    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            const double temperature = Tmin + i * deltaT;

            double lceq[5] = { init_guess[0], init_guess[1], // liquid
                init_guess[2], init_guess[3], // solid
                init_guess[4] };

            // compute equilibrium concentrations in each phase
            // at ends of tie line
            bool found_ceq = cafe.computeTieLine(
                temperature, nominalc[0], nominalc[1], &lceq[0]);
            if (lceq[0] > 1.) found_ceq = false;
            if (lceq[0] < 0.) found_ceq = false;
            if (lceq[1] > 1.) found_ceq = false;
            if (lceq[1] < 0.) found_ceq = false;

            std::cout << "Temperature = " << temperature << std::endl;
            if (found_ceq)
            {
                std::cout << "Found solution: " << lceq[0] << ", " << lceq[1]
                          << ", " << lceq[2] << ", " << lceq[3] << ", "
                          << lceq[4] << "..." << std::endl;
            }
            CHECK(found_ceq);

            cel[2 * i]     = lceq[0];
            cel[2 * i + 1] = lceq[1];
            ces[2 * i]     = lceq[2];
            ces[2 * i + 1] = lceq[3];
        }
    }

// parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        double lceq[5] = { init_guess[0], init_guess[1], // liquid
            init_guess[2], init_guess[3], // solid
            init_guess[4] };

        // compute equilibrium concentrations in each phase
        bool found_ceq = cafe.computeTieLine(
            temperature, nominalc[0], nominalc[1], &lceq[0]);

        CHECK(found_ceq);
        CHECK(lceq[0] == Approx(cel[2 * i]).margin(1.e-6));
        CHECK(lceq[1] == Approx(cel[2 * i + 1]).margin(1.e-6));
        CHECK(lceq[2] == Approx(ces[2 * i]).margin(1.e-6));
        CHECK(lceq[3] == Approx(ces[2 * i + 1]).margin(1.e-6));
    }
}
