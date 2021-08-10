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

    // nominal concentration for which to solve KKS eqs.
    const double nominalc = 0.2;

    // initial guesses
    const double init_guess[2] = { 0.5, 0.5 };

    std::vector<double> cl(nTintervals + 1);
    std::vector<double> cs(nTintervals + 1);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary* cafe
        = new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(calphad_db,
            newton_db, energy_interp_func_type, conc_interp_func_type);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            const double temperature = Tmin + i * deltaT;

            double phi = 0.5;
            double conc[2];

            // compute concentrations in each phase
            cafe->computePhaseConcentrations(
                temperature, &nominalc, &phi, conc);

            std::cout << "Temperature = " << temperature << std::endl;
            std::cout << "Concentrations: cl = " << conc[0]
                      << " and cs = " << conc[1] << "..." << std::endl;

            cl[i] = conc[0];
            cs[i] = conc[1];
        }
    }

    std::vector<short> nitscpu(nTintervals + 1);

// parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        double phi     = 0.5;
        double conc[2] = { init_guess[0], init_guess[1] };

        // compute concentrations in each phase
        nitscpu[i] = cafe->computePhaseConcentrations(
            temperature, &nominalc, &phi, conc);

        CHECK(conc[0] == Approx(cl[i]).margin(1.e-6));
        CHECK(conc[1] == Approx(cs[i]).margin(1.e-6));
    }

    short* nits  = new short[nTintervals + 1];
    double* xdev = new double[2 * (nTintervals + 1)];

// clang-format off
#pragma omp target map(to : cafe [0:1]) \
                   map(to : init_guess[:2]) \
                   map(from : xdev[:2 * (nTintervals + 1)]) \
                   map(from : nits[:nTintervals + 1])
// clang-format on
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        double phi = 0.5;
        double c0  = nominalc;

        xdev[2 * i]     = init_guess[0];
        xdev[2 * i + 1] = init_guess[1];

        // compute concentrations in each phase
        nits[i] = cafe->computePhaseConcentrations(
            temperature, &c0, &phi, &xdev[2 * i]);
    }

    for (int i = 0; i < nTintervals + 1; i++)
    {
        std::cout << "Number of Newton iterations: " << nits[i] << std::endl;
        std::cout << "Device: x=" << xdev[2 * i] << "," << xdev[2 * i + 1]
                  << std::endl;
        CHECK(xdev[2 * i] == Approx(cl[i]).margin(1.e-6));
        CHECK(xdev[2 * i + 1] == Approx(cs[i]).margin(1.e-6));
        CHECK(nits[i] == nitscpu[i]);
    }
    delete[] nits;
    delete[] xdev;
}
