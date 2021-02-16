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

TEST_CASE("CALPHAD ternary kks in a loop", "[ternary kks loop]")
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

    std::cout << "Read CALPHAD database..." << std::endl;
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
    const double init_guess[4] = { 0.33, 0.38, 0.32, 0.33 };
    double nominalc[2]         = { 0.33, 0.33 };

    std::vector<double> cl(2 * (nTintervals + 1));
    std::vector<double> cs(2 * (nTintervals + 1));

    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            const double temperature = Tmin + i * deltaT;

            double phi = 0.5;

            double conc[4] = { init_guess[0], init_guess[1], // liquid
                init_guess[2], init_guess[3] }; // solid

            // compute equilibrium concentrations in each phase
            cafe.computePhaseConcentrations(temperature, nominalc, phi, conc);

            std::cout << "Temperature = " << temperature << std::endl;
            std::cout << "Concentrations: cl = (" << conc[0] << "." << conc[1]
                      << ")"
                      << " and cs = (" << conc[2] << "," << conc[3] << ")"
                      << std::endl;

            cl[2 * i]     = conc[0];
            cl[2 * i + 1] = conc[1];
            cs[2 * i]     = conc[2];
            cs[2 * i + 1] = conc[3];
        }
    }

// parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;
        double phi               = 0.5;

        double conc[4] = { init_guess[0], init_guess[1], // liquid
            init_guess[2], init_guess[3] }; // solid

        // compute equilibrium concentrations in each phase
        cafe.computePhaseConcentrations(temperature, nominalc, phi, &conc[0]);

        CHECK(conc[0] == Approx(cl[2 * i]).margin(1.e-6));
        CHECK(conc[1] == Approx(cl[2 * i + 1]).margin(1.e-6));
        CHECK(conc[2] == Approx(cs[2 * i]).margin(1.e-6));
        CHECK(conc[3] == Approx(cs[2 * i + 1]).margin(1.e-6));
    }
}
