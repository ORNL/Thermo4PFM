#define CATCH_CONFIG_MAIN

#include "InterpolationType.h"
#include "functions.h"

#include "catch.hpp"

TEST_CASE("Interpolation functions", "[interpolation]")
{
    double (*fun_ptr_arr_[3])(const double){ Thermo4PFM::linear_interp_func,
        Thermo4PFM::pbg_interp_func, Thermo4PFM::harmonic_interp_func };

    // make sure all the functions equal 0.5 for phi=0.5
    for (int i = 0; i < 3; i++)
    {
        const double hphi = fun_ptr_arr_[i](0.5);
        CHECK(hphi == Approx(0.5).margin(1.e-8));
    }

    const double phi = 0.25;

    double hphi = interp_func(Thermo4PFM::EnergyInterpolationType::LINEAR, phi);
    CHECK(hphi == Approx(phi).margin(1.e-8));

    hphi = interp_func(Thermo4PFM::EnergyInterpolationType::PBG, phi);
    CHECK(hphi
          == Approx(phi * phi * phi * (10. - 15. * phi + 6. * phi * phi))
                 .margin(1.e-8));

    hphi = interp_func(Thermo4PFM::EnergyInterpolationType::HARMONIC, phi);
    CHECK(hphi == Approx(phi * phi * (3. - 2. * phi)).margin(1.e-8));
}
