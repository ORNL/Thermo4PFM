#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "CALPHADFunctions.h"

#include <cmath>
#include <iomanip>
#include <iostream>

TEST_CASE("CALPHAD functions", "[calphad functions]")
{
    std::cout << "Test CALPHAD functions." << std::endl;

    const double epsilon = 1.e-8;
    const double tol     = 1.e-6;

    std::cout << std::setprecision(12);
    std::cerr << std::setprecision(12);

    std::cout << "============" << std::endl;
    std::cout << "Binaries..." << std::endl;
    {
        double l0 = 2.3;
        double l1 = 5.1;
        double l2 = 3.2;
        double l3 = -2.5;
        double c  = 0.33;

        double f0 = Thermo4PFM::CALPHADcomputeFMixBinary(l0, l1, l2, l3, c);
        double f1
            = Thermo4PFM::CALPHADcomputeFMixBinary(l0, l1, l2, l3, c + epsilon);

        double fd = (f1 - f0) / epsilon;
        std::cout << "Numerical derivative   = " << fd << std::endl;
        double ad
            = Thermo4PFM::CALPHADcomputeFMix_derivBinary(l0, l1, l2, l3, c);
        std::cout << "Analytical derivative = " << ad << std::endl;

        REQUIRE(fd == Approx(ad).epsilon(tol));
    }

    std::cout << "======================" << std::endl;
    std::cout << "FIdealMix Ternaries..." << std::endl;
    {
        double rt = 10.;
        double cA = 0.33;
        double cB = 0.21;
        double f0 = Thermo4PFM::CALPHADcomputeFIdealMixTernary(rt, cA, cB);
        double f1
            = Thermo4PFM::CALPHADcomputeFIdealMixTernary(rt, cA + epsilon, cB);
        double f2
            = Thermo4PFM::CALPHADcomputeFIdealMixTernary(rt, cA, cB + epsilon);
        double deriv[2];
        Thermo4PFM::CALPHADcomputeFIdealMix_derivTernary(rt, cA, cB, deriv);
        std::cout << "Analytical derivative = " << deriv[0] << "," << deriv[1]
                  << std::endl;

        double fd1 = (f1 - f0) / epsilon;
        std::cout << "Numerical derivative   = " << fd1 << std::endl;
        REQUIRE(fd1 == Approx(deriv[0]).epsilon(tol));

        double fd2 = (f2 - f0) / epsilon;
        std::cout << "Numerical derivative   = " << fd2 << std::endl;
        REQUIRE(fd2 == Approx(deriv[1]).epsilon(tol));
    }

    std::cout << "=================" << std::endl;
    std::cout << "FMix Ternaries..." << std::endl;
    {
        CalphadDataType lAB[4]  = { 1.2, 3.5, 6.1, -1.2 };
        CalphadDataType lAC[4]  = { 3.2, 4.5, -3.1, 7.1 };
        CalphadDataType lBC[4]  = { 0.2, 0.7, 4.1, -8.2 };
        CalphadDataType lABC[3] = { 3.7, 6.3, -1.1 };

        double cA = 0.33;
        double cB = 0.21;

        double f0 = Thermo4PFM::CALPHADcomputeFMixTernary(
            lAB, lAC, lBC, lABC, cA, cB);
        double f1 = Thermo4PFM::CALPHADcomputeFMixTernary(
            lAB, lAC, lBC, lABC, cA + epsilon, cB);
        double f2 = Thermo4PFM::CALPHADcomputeFMixTernary(
            lAB, lAC, lBC, lABC, cA, cB + epsilon);

        double deriv[2];
        Thermo4PFM::CALPHADcomputeFMix_derivTernary(
            lAB, lAC, lBC, lABC, cA, cB, deriv);
        std::cout << "Analytical derivative = " << deriv[0] << "," << deriv[1]
                  << std::endl;

        double fd1 = (f1 - f0) / epsilon;
        std::cout << "Numerical derivative   = " << fd1 << std::endl;

        REQUIRE(fd1 == Approx(deriv[0]).epsilon(tol));

        double fd2 = (f2 - f0) / epsilon;
        std::cout << "Numerical derivative   = " << fd2 << std::endl;

        REQUIRE(fd2 == Approx(deriv[1]).epsilon(tol));

        double deriv2[4];
        Thermo4PFM::CALPHADcomputeFMix_deriv2Ternary(
            lAB, lAC, lBC, lABC, cA, cB, deriv2);

        double fderiv0[2];
        Thermo4PFM::CALPHADcomputeFMix_derivTernary(
            lAB, lAC, lBC, lABC, cA + epsilon, cB, fderiv0);
        double fderiv1[2];
        Thermo4PFM::CALPHADcomputeFMix_derivTernary(
            lAB, lAC, lBC, lABC, cA, cB + epsilon, fderiv1);

        double fd = (fderiv0[0] - deriv[0]) / epsilon;
        std::cout << "FD=" << fd << ", exact=" << deriv2[0] << std::endl;

        REQUIRE(fd == Approx(deriv2[0]).epsilon(tol));

        fd = (fderiv1[1] - deriv[1]) / epsilon;
        std::cout << "FD=" << fd << ", exact=" << deriv2[3] << std::endl;

        REQUIRE(fd == Approx(deriv2[3]).epsilon(tol));

        // cross derivatives
        fd = (fderiv0[1] - deriv[1]) / epsilon;
        std::cout << "FD=" << fd << ", exact=" << deriv2[1] << std::endl;

        REQUIRE(fd == Approx(deriv2[1]).epsilon(tol));

        fd = (fderiv1[0] - deriv[0]) / epsilon;
        std::cout << "FD=" << fd << ", exact=" << deriv2[2] << std::endl;

        REQUIRE(fd == Approx(deriv2[2]).epsilon(tol));
    }
}
