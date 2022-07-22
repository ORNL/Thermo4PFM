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

TEST_CASE("CALPHAD binary energy", "[binary energy]")
{
#ifdef _OPENMP
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;
#endif

    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::PBG;

    const double Tmin     = 300.;
    const double Tmax     = 2999.;
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

    std::vector<double> e(nTintervals + 1);

    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary* cafe
        = new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(calphad_db,
            newton_db, energy_interp_func_type, conc_interp_func_type);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            const double temperature = Tmin + i * deltaT;

            // compute energy for species 0
            e[i] = cafe->getFenergyPhaseL(0, temperature);

            std::cout << "Temperature = " << temperature
                      << ", energy = " << e[i] << std::endl;
        }
    }

// parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        // compute energy for species 0
        double etest = cafe->getFenergyPhaseL(0, temperature);

        CHECK(etest == Approx(e[i]).margin(1.e-6));
    }

#ifdef HAVE_OPENMP_OFFLOAD
    double* edev = new double[nTintervals + 1];

// clang-format off
#pragma omp target map(to : cafe [0:1]) \
                   map(from : edev[:(nTintervals + 1)])
// clang-format on
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        // compute concentrations in each phase
        edev[i] = cafe->getFenergyPhaseL(0, temperature);
    }

    for (int i = 0; i < nTintervals + 1; i++)
    {
        std::cout << "Device: e=" << edev[i] << std::endl;
        CHECK(edev[i] == Approx(e[i]).margin(1.e-6));
    }
    delete[] edev;
#endif

    delete cafe;
}
