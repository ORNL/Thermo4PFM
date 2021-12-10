#define CATCH_CONFIG_MAIN

#include "CALPHADSpeciesPhaseGibbsEnergyExpansion.h"

#include "catch.hpp"

#include <iostream>
#include <omp.h>

TEST_CASE("CALPHAD Gibbs energy expansion", "[gibbs energy expansion]")
{
    std::cout << std::setprecision(12);
    double temperature = 100.;

    double ehost = -1.;
    {
        Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergyExpansion<double> cspgee;

        cspgee.init(1., 1., 1., 1., 1., 1., 1., 1., 1.);

        ehost = cspgee.value(temperature);
        std::cout << "Expansion value =" << ehost << std::endl;
    }

#ifdef HAVE_OPENMP_OFFLOAD
    double edev = -1.;

    {
#pragma omp target map(tofrom : edev) map(to : temperature)
        {
            Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergyExpansion<double>
                cspgee_dev;
            cspgee_dev.init(1., 1., 1., 1., 1., 1., 1., 1., 1.);
            edev = cspgee_dev.value(temperature);
        }

        // check result
        std::cout << "Expansion value on device =" << edev << std::endl;
        CHECK(std::abs(edev - ehost) < 1.e-6);
    }

    // do it again with an object constructed on the host
    edev = -1.;
    {
        Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergyExpansion<double> cspgee;
        cspgee.init(1., 1., 1., 1., 1., 1., 1., 1., 1.);

#pragma omp target map(to : cspgee) map(tofrom : edev) map(to : temperature)
        {
            edev = cspgee.value(temperature);
        }
        // check result
        std::cout << "Expansion value on device =" << edev << std::endl;
        CHECK(std::abs(edev - ehost) < 1.e-6);
    }
#endif
}
