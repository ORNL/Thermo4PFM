#define CATCH_CONFIG_MAIN

#include "KKSFreeEnergyFunctionDiluteBinary.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

namespace pt = boost::property_tree;

TEST_CASE("Loop dilute binary KKS", "[loop dilute binary kks]")
{
    Thermo4PFM::EnergyInterpolationType energy_interp_func_type
        = Thermo4PFM::EnergyInterpolationType::PBG;
    Thermo4PFM::ConcInterpolationType conc_interp_func_type
        = Thermo4PFM::ConcInterpolationType::LINEAR;

    const double Tmin     = 1600.;
    const double Tmax     = 1720.;
    const int nTintervals = 10;
    const double deltaT   = (Tmax - Tmin) / (double)nTintervals;

    // read dilute alloy informaton
    pt::ptree conc_db;
    try
    {
        pt::read_json("dilute_binary.json", conc_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }
    pt::write_json(std::clog, conc_db);

    Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary* cafe
        = new Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary(
            conc_db, energy_interp_func_type, conc_interp_func_type);

    // initial guesses
    double c_init0 = 0.5;
    double c_init1 = 0.5;

    // solve KKS equations
    double nominalc = 0.05;
    double phi      = 0.5;

    std::vector<double> cl(nTintervals + 1);
    std::vector<double> cs(nTintervals + 1);

    {
        double sol[2];

        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            double temperature = Tmin + i * deltaT;
            ;

            sol[0] = c_init0;
            sol[1] = c_init1;
            cafe->computePhaseConcentrations(
                temperature, &nominalc, &phi, &sol[0]);

            std::cout << "-------------------------------" << std::endl;
            std::cout << "Temperature = " << temperature << std::endl;
            std::cout << "   cL = " << sol[0] << ",   cS = " << sol[1]
                      << std::endl;

            cl[i] = sol[0];
            cs[i] = sol[1];
        }
    }

    std::vector<short> nitscpu(nTintervals + 1);

    // parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        double conc[2] = { c_init0, c_init1 };

        // compute concentrations in each phase
        nitscpu[i] = cafe->computePhaseConcentrations(
            temperature, &nominalc, &phi, conc);

        CHECK(conc[0] == Approx(cl[i]).margin(1.e-6));
        CHECK(conc[1] == Approx(cs[i]).margin(1.e-6));
    }

    short* nits          = new short[nTintervals + 1];
    double* xdev         = new double[2 * (nTintervals + 1)];
    double init_guess[2] = { c_init0, c_init1 };

    std::cout << "=========================\n";
    std::cout << "GPU offload testing:" << std::endl;

// clang-format off
#pragma omp target map(to:cafe[0:1]) \
                   map(to: init_guess) \
                   map(from: xdev[:2*(nTintervals+1)]) \
                   map(from: nits[:nTintervals+1])
// clang-format on
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double temperature = Tmin + i * deltaT;

        double phi      = 0.5;
        double nominalc = 0.05;

        xdev[2 * i]     = init_guess[0];
        xdev[2 * i + 1] = init_guess[1];

        // compute concentrations in each phase
        nits[i] = cafe->computePhaseConcentrations(
            temperature, &nominalc, &phi, &xdev[2 * i]);
    }

    // verify results
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

    delete cafe;
}
