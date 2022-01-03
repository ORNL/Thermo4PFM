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

TEST_CASE("CALPHAD binary equilibrium", "[binary equilibrium]")
{
#ifdef _OPENMP
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;
#endif

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

    std::vector<double> cel(nTintervals + 1);
    std::vector<double> ces(nTintervals + 1);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary* cafe
        = new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(calphad_db,
            newton_db, energy_interp_func_type, conc_interp_func_type);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            const double temperature = Tmin + i * deltaT;

            double lceq[2] = { init_guess[0], init_guess[1] };

            // compute equilibrium concentrations in each phase
            bool found_ceq = cafe->computeCeqT(temperature, &lceq[0]);
            if (lceq[0] > 1.) found_ceq = false;
            if (lceq[0] < 0.) found_ceq = false;
            if (lceq[1] > 1.) found_ceq = false;
            if (lceq[1] < 0.) found_ceq = false;

            std::cout << "Temperature = " << temperature << std::endl;
            if (found_ceq)
            {
                std::cout << "Found equilibrium concentrations: " << lceq[0]
                          << " and " << lceq[1] << "..." << std::endl;
            }
            CHECK(found_ceq);

            cel[i] = lceq[0];
            ces[i] = lceq[1];
        }
    }

// parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        double lceq[2] = { init_guess[0], init_guess[1] };

        // compute equilibrium concentrations in each phase
        bool found_ceq = cafe->computeCeqT(temperature, &lceq[0]);

        CHECK(found_ceq);
        CHECK(lceq[0] == Approx(cel[i]).margin(1.e-6));
        CHECK(lceq[1] == Approx(ces[i]).margin(1.e-6));
    }

#ifdef HAVE_OPENMP_OFFLOAD
    double* xdev = new double[2 * (nTintervals + 1)];

// clang-format off
#pragma omp target map(to : cafe [0:1]) \
                   map(to : init_guess[:2]) \
                   map(from : xdev[:2 * (nTintervals + 1)])
// clang-format on
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        xdev[2 * i]     = init_guess[0];
        xdev[2 * i + 1] = init_guess[1];

        // compute concentrations in each phase
        cafe->computeCeqT(temperature, &xdev[2 * i]);
    }

    for (int i = 0; i < nTintervals; i++)
    {
        std::cout << "Device: x=" << xdev[2 * i] << "," << xdev[2 * i + 1]
                  << std::endl;
        CHECK(xdev[2 * i] == Approx(cel[i]).margin(1.e-6));
        CHECK(xdev[2 * i + 1] == Approx(ces[i]).margin(1.e-6));
    }
    delete[] xdev;
#endif

    delete cafe;
}
