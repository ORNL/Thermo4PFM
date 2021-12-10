#include "CALPHADConcSolverBinary.h"
#include "CALPHADFunctions.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "PhysicalConstants.h"

#include <chrono>

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iomanip>
#include <iostream>
#include <string>

#include <omp.h>

namespace pt = boost::property_tree;

typedef std::chrono::high_resolution_clock Clock;

int main(int argc, char* argv[])
{
    const int N = 10000000;

#ifdef _OPENMP
    std::cout << "Compiled by an OpenMP-compliant implementation.\n";
    std::cout << "Run with " << omp_get_max_threads() << " threads"
              << std::endl;
#endif

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

    double temperature = 1450.;

    CalphadDataType LmixPhaseL[4][MAX_POL_T_INDEX];
    CalphadDataType LmixPhaseA[4][MAX_POL_T_INDEX];

    {
        std::string dbnamemixL("LmixPhaseL");
        pt::ptree Lmix0_db = calphad_db.get_child(dbnamemixL);
        Thermo4PFM::readLmixBinary(Lmix0_db, LmixPhaseL);
    }
    {
        std::string dbnamemixA("LmixPhaseA");
        pt::ptree Lmix1_db = calphad_db.get_child(dbnamemixA);
        Thermo4PFM::readLmixBinary(Lmix1_db, LmixPhaseA);
    }

    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy g_species_phaseL[2];
    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy g_species_phaseA[2];

    {
        pt::ptree& species0_db = calphad_db.get_child("SpeciesA");
        std::string dbnameL("PhaseL");
        std::string dbnameA("PhaseA");

        g_species_phaseL[0].initialize("L0", species0_db.get_child(dbnameL));
        g_species_phaseA[0].initialize("A0", species0_db.get_child(dbnameA));

        pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
        g_species_phaseL[1].initialize("L1", speciesB_db.get_child(dbnameL));
        g_species_phaseA[1].initialize("A1", speciesB_db.get_child(dbnameA));
    }

    CalphadDataType fA[2];
    fA[0] = g_species_phaseL[0].fenergy(temperature);
    fA[1] = g_species_phaseA[0].fenergy(temperature);
    // std::cout<<"fA[0]="<<fA[0]<<", fA[1]="<<fA[1]<<std::endl;

    CalphadDataType fB[2];
    fB[0] = g_species_phaseL[1].fenergy(temperature);
    fB[1] = g_species_phaseA[1].fenergy(temperature);
    // std::cout<<"fB[0]="<<fB[0]<<", fB[1]="<<fB[1]<<std::endl;

    CalphadDataType Lmix_L[4];
    for (int i = 0; i < 4; i++)
        Lmix_L[i] = LmixPhaseL[i][0] + temperature * LmixPhaseL[i][1];
    // for(int i=0;i<4;i++)std::cout<<"Lmix_L["<<i<<"]="<<Lmix_L[i]<<std::endl;

    CalphadDataType Lmix_A[4];
    for (int i = 0; i < 4; i++)
        Lmix_A[i] = LmixPhaseA[i][0] + temperature * LmixPhaseA[i][1];
    // for(int i=0;i<4;i++)std::cout<<"Lmix_A["<<i<<"]="<<Lmix_A[i]<<std::endl;

    const double RTinv
        = 1.0 / (GASCONSTANT_R_JPKPMOL * temperature);

    double sol[2] = { 0.5, 0.5 };

    double deviation = 1.e-4;

    double* xhost = new double[2 * N];
    for (int i = 0; i < 2 * N; i++)
    {
        xhost[i] = -1.;
    }

    // Host solve
    {
        short* nits = new short[N];
        auto t1     = Clock::now();

#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
#ifdef _OPENMP
            if (!omp_is_initial_device()) abort();
#endif
            xhost[2 * i]     = sol[0];
            xhost[2 * i + 1] = sol[1];
            double hphi      = 0.5 + (i % 100) * deviation;
            double c0        = 0.3;
            Thermo4PFM::CALPHADConcSolverBinary solver;
            solver.setup(c0, hphi, RTinv, Lmix_L, Lmix_A, fA, fB);
            nits[i] = solver.ComputeConcentration(&xhost[2 * i], 1.e-8, 50);
        }
        auto t2 = Clock::now();
        long int usec
            = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                  .count();
        std::cout << "Host time/us/solve:   " << (double)usec / (double)N
                  << std::endl;

        std::cout << std::setprecision(12);

        int n = N > 20 ? 20 : N;

        for (int i = 0; i < n; i++)
        {
            std::cout << "Host: x=" << xhost[2 * i] << "," << xhost[2 * i + 1]
                      << std::endl;
            std::cout << "nits=" << nits[i] << std::endl;
        }
        delete[] nits;
    }
#ifdef HAVE_OPENMP_OFFLOAD
    // Device solve
    {
        double* xdev = new double[2 * N];
        for (int i = 0; i < 2 * N; i++)
        {
            xdev[i] = -1;
        }

        short* nits = new short[N];

        auto t1 = Clock::now();

// clang-format off
#pragma omp target map(to : sol) \
                   map(tofrom : xdev) \
                   map(to : fA, fB, Lmix_L, Lmix_A) \
                   map(to : RTinv) \
                   map(from : nits)                                                     \
// clang-format on
        {
#pragma omp teams distribute parallel for
            for (int i = 0; i < N; i++)
            {
                // if( omp_is_initial_device() ) abort();
                xdev[2 * i]     = sol[0];
                xdev[2 * i + 1] = sol[1];

                double hphi = 0.5 + (i % 100) * deviation;
                double c0   = 0.3;
                class Thermo4PFM::CALPHADConcSolverBinary solver;
                solver.setup(c0, hphi, RTinv, Lmix_L, Lmix_A, fA, fB);
                nits[i] = solver.ComputeConcentration(&xdev[2 * i], 1.e-8, 50);
            }
        }

        auto t2 = Clock::now();
        long int usec
            = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                  .count();
        std::cout << "Device time/us/solve: " << (double)usec / (double)N
                  << std::endl;

        double tol = 1.e-8;
        for (int i = 0; i < N; i++)
        {
            if (std::abs(xdev[2 * i] - xhost[2 * i]) > tol
                || std::abs(xdev[2 * i + 1] - xhost[2 * i + 1]) > tol || N < 20)
            {
                std::cout << "Device: x=" << xdev[2 * i] << ","
                          << xdev[2 * i + 1] << std::endl;
                std::cout << "Difference: " << xdev[2 * i] - xhost[2 * i]
                          << ", " << xdev[2 * i + 1] - xhost[2 * i + 1]
                          << std::endl;
                std::cout << "nits[" << i << "]=" << nits[i] << std::endl;
            }
        }
        delete[] nits;
        delete[] xdev;
    }
#endif
    delete[] xhost;
}
