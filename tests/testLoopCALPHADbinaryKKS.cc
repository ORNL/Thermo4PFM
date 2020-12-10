#define CATCH_CONFIG_MAIN

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "InterpolationType.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

#include <omp.h>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD binary kks in a loop", "[binary kks loop]")
{
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    const double Tmin     = 1300.;
    const double Tmax     = 1500.;
    const int nTintervals = 10;
    const double deltaT   = (Tmax - Tmin) / (double)nTintervals;

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

    // initial guesses
    const double init_guess[2] = { 0.2, 0.1 };

    std::vector<double> cl(nTintervals + 1);
    std::vector<double> cs(nTintervals + 1);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(
        calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            const double temperature = Tmin + i * deltaT;

            double phi = 0.5;
            double conc[2];

            // compute concentrations in each phase
            cafe.computePhaseConcentrations(temperature, init_guess, phi, conc);

            std::cout << "Temperature = " << temperature << std::endl;
            std::cout << "Concentrations: cl = " << conc[0]
                      << " and cs = " << conc[1] << "..." << std::endl;

            cl[i] = conc[0];
            cs[i] = conc[1];
        }
    }

// parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        double phi = 0.5;
        double conc[2];

        // compute concentrations in each phase
        cafe.computePhaseConcentrations(temperature, init_guess, phi, conc);

        REQUIRE(conc[0] == Approx(cl[i]).margin(1.e-1));
        REQUIRE(conc[1] == Approx(cs[i]).margin(1.e-1));
    }
}
