
#include "CALPHADSpeciesPhaseGibbsEnergy.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <string>

#include <omp.h>

namespace pt = boost::property_tree;

int main(int argc, char* argv[])
{
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;

    const double Tmin     = 1300.;
    const double Tmax     = 1500.;
    const int nTintervals = 10;
    const double deltaT   = (Tmax - Tmin) / (double)nTintervals;

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("calphadAuNi.json", calphad_db);
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

    std::vector<double> el(nTintervals + 1);

    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy energyNiLiquid;
    energyNiLiquid.initialize("NiL", liqB);

    {
        // serial loop
        for (int i = 0; i < nTintervals + 1; i++)
        {
            const double temperature = Tmin + i * deltaT;

            double energy = energyNiLiquid.fenergy(temperature);

            std::cout << "Temperature = " << temperature << std::endl;
            std::cout << "Energy liquid = " << energy << "..." << std::endl;
            el[i] = energy;
        }
    }

    double* vel = new double[nTintervals];
#pragma omp target map(to : Tmin, deltaT) map(tofrom : vel)
    {
#pragma omp parallel for
        for (int i = 0; i < nTintervals + 1; i++)
        {
            double temperature = Tmin + i * deltaT;

            vel[i] = energyNiLiquid.fenergy(temperature);
        }
    }

    for (int i = 0; i < nTintervals + 1; i++)
        std::cout << "Energy liquid = " << vel[i] << std::endl;

    delete[] vel;
}
