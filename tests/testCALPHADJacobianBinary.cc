#define CATCH_CONFIG_MAIN

#include "CALPHADConcSolverBinary.h"
#include "CALPHADConcSolverBinaryThreePhase.h"
#include "CALPHADEqConcSolverBinary.h"

#include "InterpolationType.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "catch.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

namespace pt = boost::property_tree;

TEST_CASE("CALPHAD Jacobian binary", "[CALPHAD Jacobian binary]")
{
    // for finite differences
    double epsilon = 1.e-8;
    double tol     = 1.e-6;
    double RT      = 0.1;
    double RTinv   = 1. / RT;

    {
        double c0   = 0.1;
        double hphi = 0.3;

        // energies of 2 species, in two phase each
        CalphadDataType fA[2] = { 2.3, 4.5 };
        CalphadDataType fB[2] = { 0.5, 2.6 };

        // L coefficients for 2 possible phases (L and S)
        CalphadDataType Lmix_L[4] = { 2., 3., 4., 5. };
        CalphadDataType Lmix_S[4] = { 2., 2., 3., 4. };

        double fvec1[2];
        double fvec2[2];
        JacobianDataType* fjac[2];
        for (int i = 0; i < 2; i++)
            fjac[i] = new JacobianDataType[2];

        double x[2] = { 0.1, 0.4 };

        std::cout << "Test CALPHADConcSolverBinary...\n";
        {
            Thermo4PFM::CALPHADConcSolverBinary solver;
            solver.setup(c0, hphi, RTinv, Lmix_L, Lmix_S, fA, fB);

            solver.RHS(x, fvec1);

            solver.Jacobian(x, fjac);

            std::cout << std::setprecision(12);

            // loop over variables
            for (int j = 0; j < 2; j++)
            {
                std::cout << "----------------------------" << std::endl;
                std::cout << "Test variations of variable " << j << std::endl;
                x[j] += epsilon;
                solver.RHS(x, fvec2);

                double fd[2];
                // loop over equations
                for (int i = 0; i < 2; i++)
                    fd[i] = (fvec2[i] - fvec1[i]) / epsilon;

                for (int i = 0; i < 2; i++)
                {
                    std::cout << "Equation " << i << ", FD=" << fd[i]
                              << ", fjac=" << fjac[i][j] << std::endl;
                    CHECK(fd[i] == Approx(fjac[i][j]).margin(tol));
                }

                x[j] -= epsilon;
            }
        }

        std::cout << "Test CALPHADEqConcSolverBinary...\n";
        {
            Thermo4PFM::CALPHADEqConcSolverBinary solver;
            solver.setup(RTinv, Lmix_L, Lmix_S, fA, fB);

            solver.RHS(x, fvec1);

            solver.Jacobian(x, fjac);

            std::cout << std::setprecision(12);

            // loop over variables
            for (int j = 0; j < 2; j++)
            {
                std::cout << "----------------------------" << std::endl;
                std::cout << "Test variations of variable " << j << std::endl;
                x[j] += epsilon;
                solver.RHS(x, fvec2);

                double fd[2];
                // loop over equations
                for (int i = 0; i < 2; i++)
                    fd[i] = (fvec2[i] - fvec1[i]) / epsilon;

                for (int i = 0; i < 2; i++)
                {
                    std::cout << "Equation " << i << ", FD=" << fd[i]
                              << ", fjac=" << fjac[i][j] << std::endl;
                    CHECK(fd[i] == Approx(fjac[i][j]).margin(tol));
                }

                x[j] -= epsilon;
            }
        }

        for (int i = 0; i < 2; i++)
            delete[] fjac[i];
    }

    // CALPHADConcSolverBinaryThreePhase
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Test CALPHADConcSolverBinaryThreePhase...\n";
    {
        double c0 = 0.3;

        // energies of 2 species, in two phase each
        CalphadDataType fA[2] = { 2.3, 4.5 };
        CalphadDataType fB[2] = { 0.5, 2.6 };

        // L coefficients for 2 possible phases (L and S)
        CalphadDataType Lmix_L[4]  = { 2., 3., 4., 5. };
        CalphadDataType Lmix_S0[4] = { 2., 2., 3., 4. };
        CalphadDataType Lmix_S1[4] = { 2., 2., 3., 4. };

        double fvec1[3];
        double fvec2[3];
        JacobianDataType* fjac[3];
        for (int i = 0; i < 3; i++)
            fjac[i] = new JacobianDataType[3];

        // test all 3 cases in class CALPHADConcSolverBinaryThreePhase
        // depending on value of hphi1 and hphi2
        for (int i = 0; i < 3; i++)
        {
            double hphi0 = 0.53 - i * 0.1;
            double hphi1 = 0.43 - i * 0.1;

            std::cout << "----------------------------" << std::endl;
            std::cout << "hphi0 = " << hphi0 << ", hphi1 = " << hphi1
                      << ", hphi2 = " << 1. - hphi0 - hphi1 << std::endl;

            // concentrations in 3 phases we are solving for
            double x[3] = { 0.1, 0.4, 0.3 };

            Thermo4PFM::CALPHADConcSolverBinaryThreePhase solver;
            solver.setup(c0, hphi0, hphi1, 1. - hphi0 - hphi1, RT, Lmix_L,
                Lmix_S0, Lmix_S1, fA, fB);

            solver.RHS(x, fvec1);

            solver.Jacobian(x, fjac);

            std::cout << std::setprecision(12);

            // loop over variables
            for (int j = 0; j < 3; j++)
            {
                std::cout << "----------------------------" << std::endl;
                std::cout << "Test variations of variable " << j << std::endl;
                x[j] += epsilon;
                solver.RHS(x, fvec2);

                double fd[3];
                // loop over equations
                for (int i = 0; i < 3; i++)
                    fd[i] = (fvec2[i] - fvec1[i]) / epsilon;

                for (int i = 0; i < 3; i++)
                {
                    std::cout << "Equation " << i << ", FD=" << fd[i]
                              << ", fjac=" << fjac[i][j] << std::endl;
                    CHECK(fd[i] == Approx(fjac[i][j]).margin(tol));
                }

                x[j] -= epsilon;
            }
        }

        for (int i = 0; i < 3; i++)
            delete[] fjac[i];
    }
}
