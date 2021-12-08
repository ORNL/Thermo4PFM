#define CATCH_CONFIG_MAIN

#include "CALPHADSpeciesPhaseGibbsEnergy.h"

#include "catch.hpp"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

namespace pt = boost::property_tree;

TEST_CASE("Loop CALPHADSpeciesPhaseGibbsEnergy",
    "[loop CALPHADSpeciesPhaseGibbsEnergy]")
{
    std::cout << std::setprecision(12);

    // Ni melting temperature
    double temperature = 1728.;

    // read calphad database
    std::string calphad_filename("../thermodynamic_data/calphadAuNi.json");
    pt::ptree calphad_db;
    std::cout << "Read " << calphad_filename << std::endl;
    try
    {
        pt::read_json(calphad_filename, calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    // Ni is species "B" in calphadAuNi.json
    std::cout << "Read species B data" << std::endl;
    pt::ptree spB;
    try
    {
        spB = calphad_db.get_child("SpeciesB");
    }
    catch (std::exception& e)
    {
        std::cerr << "Read species: exception caught: " << e.what()
                  << std::endl;
    }

    pt::ptree liqB;
    try
    {
        liqB = spB.get_child("PhaseL");
    }
    catch (std::exception& e)
    {
        std::cerr << "Read phase: exception caught: " << e.what() << std::endl;
    }

    // liquid
    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy energyNiLiquid;
    energyNiLiquid.initialize("NiL", liqB);

    const double Tmin     = temperature - 1.;
    const double Tmax     = temperature + 1.;
    const int nTintervals = 10;
    const double deltaT   = (Tmax - Tmin) / (double)nTintervals;

    // initial guesses
    double c_init0 = 0.5;
    double c_init1 = 0.5;

    // solve KKS equations
    double nominalc = 0.05;
    double phi      = 0.5;

    std::vector<double> energy(nTintervals + 1);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            double t  = Tmin + i * deltaT;
            energy[i] = energyNiLiquid.fenergy(t);

            std::cout << "-------------------------------" << std::endl;
            std::cout << "Temperature = " << t << std::endl;
            std::cout << "energy      = " << energy[i] << std::endl;
        }
    }

    // parallel loop
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double t = Tmin + i * deltaT;
        double e       = energyNiLiquid.fenergy(t);

        CHECK(e == Approx(energy[i]).margin(1.e-6));
    }

    double* xdev = new double[nTintervals + 1];

    std::cout << "=========================\n";
    std::cout << "GPU offload testing:" << std::endl;

// clang-format off
#pragma omp target map(to:energyNiLiquid) \
                   map(from: xdev[:(nTintervals+1)])
// clang-format on
#pragma omp parallel for
    for (int i = 0; i < nTintervals + 1; i++)
    {
        const double t = Tmin + i * deltaT;
        double e       = energyNiLiquid.fenergy(t);

        xdev[i] = e;
    }

    // verify results
    for (int i = 0; i < nTintervals + 1; i++)
    {
        std::cout << "Device: x=" << xdev[i] << "," << std::endl;
        CHECK(xdev[i] == Approx(energy[i]).margin(1.e-6));
    }

    delete[] xdev;
}
