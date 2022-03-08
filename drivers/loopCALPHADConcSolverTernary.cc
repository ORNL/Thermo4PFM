#include "CALPHADConcSolverTernary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
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
#include <sys/time.h>
#include <time.h>

#include <omp.h>

namespace pt = boost::property_tree;

typedef std::chrono::high_resolution_clock Clock;

double gtod(void)
{
    struct timeval tv;
    gettimeofday(&tv, (struct timezone*)nullptr);
    return 1.e6*tv.tv_sec + tv.tv_usec;
}

int main(int argc, char* argv[])
{
    const int N = 100000;

#ifdef _OPENMP
    std::cout << "Compiled by an OpenMP-compliant implementation.\n";
    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;
#endif

    std::cout << " Read CALPHAD database..." << std::endl;
    pt::ptree calphad_db;
    try
    {
        pt::read_json("../thermodynamic_data/calphadMoNbTa.json", calphad_db);
    }
    catch (std::exception& e)
    {
        std::cerr << "exception caught: " << e.what() << std::endl;
    }

    double temperature = 2923.;

    CalphadDataType LmixABPhaseL[4][2];
    CalphadDataType LmixABPhaseA[4][2];

    CalphadDataType LmixACPhaseL[4][2];
    CalphadDataType LmixACPhaseA[4][2];

    CalphadDataType LmixBCPhaseL[4][2];
    CalphadDataType LmixBCPhaseA[4][2];

    CalphadDataType LmixABCPhaseL[3][2];
    CalphadDataType LmixABCPhaseA[3][2];

    {
        std::string dbnamemixL("LmixABCPhaseL");
        if (calphad_db.get_child_optional(dbnamemixL))
        {
            pt::ptree& Lmix0_db = calphad_db.get_child(dbnamemixL);
            Thermo4PFM::readLmixTernaryParameters(Lmix0_db, LmixABCPhaseL);
        }
        else
        {
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 2; i++)
                {
                    LmixABCPhaseL[j][i] = 0.;
                }
        }
    }
    {
        std::string dbnamemixL("LmixABCPhaseA");
        if (calphad_db.get_child_optional(dbnamemixL))
        {
            pt::ptree& Lmix0_db = calphad_db.get_child(dbnamemixL);
            Thermo4PFM::readLmixTernaryParameters(Lmix0_db, LmixABCPhaseA);
        }
        else
        {
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 2; i++)
                {
                    LmixABCPhaseA[j][i] = 0.;
                }
        }
    }

    {
        std::string dbnamemixL("LmixABPhaseL");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixL);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixABPhaseL);
    }
    {
        std::string dbnamemixA("LmixABPhaseA");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixA);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixABPhaseA);
    }
    {
        std::string dbnamemixL("LmixACPhaseL");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixL);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixACPhaseL);
    }
    {
        std::string dbnamemixA("LmixACPhaseA");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixA);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixACPhaseA);
    }
    {
        std::string dbnamemixL("LmixBCPhaseL");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixL);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixBCPhaseL);
    }
    {
        std::string dbnamemixA("LmixBCPhaseA");
        pt::ptree Lmix_db = calphad_db.get_child(dbnamemixA);
        Thermo4PFM::readLmixBinary(Lmix_db, LmixBCPhaseA);
    }

    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy g_species_phaseL[3];
    Thermo4PFM::CALPHADSpeciesPhaseGibbsEnergy g_species_phaseA[3];

    {
        std::string dbnameL("PhaseL");
        std::string dbnameA("PhaseA");

        pt::ptree& speciesA_db = calphad_db.get_child("SpeciesA");
        g_species_phaseL[0].initialize("L0", speciesA_db.get_child(dbnameL));
        g_species_phaseA[0].initialize("A0", speciesA_db.get_child(dbnameA));

        pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
        g_species_phaseL[1].initialize("L1", speciesB_db.get_child(dbnameL));
        g_species_phaseA[1].initialize("A1", speciesB_db.get_child(dbnameA));

        pt::ptree& speciesC_db = calphad_db.get_child("SpeciesC");
        g_species_phaseL[2].initialize("L2", speciesC_db.get_child(dbnameL));
        g_species_phaseA[2].initialize("A2", speciesC_db.get_child(dbnameA));
    }

    CalphadDataType fA[2];
    fA[0] = g_species_phaseL[0].fenergy(temperature);
    fA[1] = g_species_phaseA[0].fenergy(temperature);
    // std::cout<<"fA[0]="<<fA[0]<<", fA[1]="<<fA[1]<<std::endl;

    CalphadDataType fB[2];
    fB[0] = g_species_phaseL[1].fenergy(temperature);
    fB[1] = g_species_phaseA[1].fenergy(temperature);
    // std::cout<<"fB[0]="<<fB[0]<<", fB[1]="<<fB[1]<<std::endl;

    CalphadDataType fC[2];
    fC[0] = g_species_phaseL[2].fenergy(temperature);
    fC[1] = g_species_phaseA[2].fenergy(temperature);

    CalphadDataType L_AB_L[4];
    for (int i = 0; i < 4; i++)
        L_AB_L[i] = LmixABPhaseL[i][0] + temperature * LmixABPhaseL[i][1];

    CalphadDataType L_AB_S[4];
    for (int i = 0; i < 4; i++)
        L_AB_S[i] = LmixABPhaseA[i][0] + temperature * LmixABPhaseA[i][1];

    CalphadDataType L_AC_L[4];
    for (int i = 0; i < 4; i++)
        L_AC_L[i] = LmixACPhaseL[i][0] + temperature * LmixACPhaseL[i][1];

    CalphadDataType L_AC_S[4];
    for (int i = 0; i < 4; i++)
        L_AC_S[i] = LmixACPhaseA[i][0] + temperature * LmixACPhaseA[i][1];

    CalphadDataType L_BC_L[4];
    for (int i = 0; i < 4; i++)
        L_BC_L[i] = LmixBCPhaseL[i][0] + temperature * LmixBCPhaseL[i][1];

    CalphadDataType L_BC_S[4];
    for (int i = 0; i < 4; i++)
        L_BC_S[i] = LmixBCPhaseA[i][0] + temperature * LmixBCPhaseA[i][1];

    CalphadDataType L_ABC_L[3];
    for (int i = 0; i < 3; i++)
        L_ABC_L[i] = LmixABCPhaseL[i][0] + temperature * LmixABCPhaseL[i][1];

    CalphadDataType L_ABC_S[3];
    for (int i = 0; i < 3; i++)
        L_ABC_S[i] = LmixABCPhaseA[i][0] + temperature * LmixABCPhaseA[i][1];

    const double RTinv
        = 1.0 / (GASCONSTANT_R_JPKPMOL * temperature);

    double sol[4] = { 0.33, 0.38, 0.32, 0.33 };

    const double deviation = 1.e-4;

    double* xhost = new double[4 * N];
    for (int i = 0; i < 4 * N; i++)
    {
        xhost[i] = -1.;
    }

    // Host solve
    {
        short* nits = new short[N];
        auto t1     = gtod(); //Clock::now();

#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            xhost[4 * i]     = sol[0];
            xhost[4 * i + 1] = sol[1];
            xhost[4 * i + 2] = sol[2];
            xhost[4 * i + 3] = sol[3];
            double hphi      = 0.5 + (i % 100) * deviation;
            double c0        = 0.33;
            double c1        = 0.33;
            Thermo4PFM::CALPHADConcSolverTernary solver;
            solver.setup(c0, c1, hphi, RTinv, L_AB_L, L_AC_L, L_BC_L, L_AB_S,
                L_AC_S, L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);
            nits[i] = solver.ComputeConcentration(&xhost[4 * i], 1.e-8, 10);
        }
        auto t2 = gtod(); //Clock::now();
        long int usec
            = t2-t1; //std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                  //.count();
        std::cout << "Host time/us/solve:   " << (double)usec / (double)N
                  << std::endl;

        std::cout << std::setprecision(12);
        int n = N > 20 ? 20 : N;
        for (int i = 0; i < n; i++)
        {
            std::cout << "Host: x=" << xhost[4 * i] << "," << xhost[4 * i + 1]
                      << "," << xhost[4 * i + 2] << "," << xhost[4 * i + 3]
                      << std::endl;
            std::cout << "nits=" << nits[i] << std::endl;
        }
        delete[] nits;
    }
    delete[] xhost;

    double* xdev = new double[4 * N];
    short* nits  = new short[N];
    for (int i = 0; i < N; i++)nits[i]=-1;

        // warm-up GPU with an empty target region
#pragma omp target
{
}

    // Device solve
    for(int rep = 0; rep < 10; rep++)
    {
        for (int i = 0; i < 4 * N; i++)
        {
            xdev[i] = -1;
        }

        auto t1 = gtod(); //Clock::now();

// clang-format off
#pragma omp target map(from : xdev[:4*N]) \
                   map(from : nits[:N])
{
#pragma omp teams distribute parallel for
            // clang-format on
            for (int i = 0; i < N; i++)
            {
                // if( omp_is_initial_device() ) abort();
                xdev[4 * i]     = sol[0];
                xdev[4 * i + 1] = sol[1];
                xdev[4 * i + 2] = sol[2];
                xdev[4 * i + 3] = sol[3];

                double hphi = 0.5 + (i % 100) * deviation;
                double c0   = 0.33;
                double c1   = 0.33;
                Thermo4PFM::CALPHADConcSolverTernary solver;
                solver.setup(c0, c1, hphi, RTinv, L_AB_L, L_AC_L, L_BC_L,
                    L_AB_S, L_AC_S, L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);
                nits[i] = solver.ComputeConcentration(&xdev[4 * i], 1.e-8, 10);
            }
        }

        auto t2 = gtod(); // Clock::now();

        long int usec
            = t2-t1; //std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                  //.count();
        std::cout << "Repetition "<<rep<<": Device time/us/solve = " << (double)usec / (double)N
                  << std::endl;

        // print out some results
        std::cout << std::setprecision(12);
        int n = N > 20 ? 20 : N;
        for (int i = 0; i < n; i++)
        {
            std::cout << "Dev: x=" << xdev[4 * i] << "," << xdev[4 * i + 1]
                      << "," << xdev[4 * i + 2] << "," << xdev[4 * i + 3]
                      << std::endl;
            std::cout << "nits=" << nits[i] << std::endl;
        }

        // verify results
        double tol = 0.03;
        int count  = 0;
        for (int i = 0; i < N; i++)
        {
            if ((xdev[4 * i + 1] != xdev[4 * i + 1])
                || std::abs(xdev[4 * i] - 0.33) > tol
                || std::abs(xdev[4 * i + 1] - 0.33) > tol
                || std::abs(xdev[4 * i + 2] - 0.33) > tol
                || std::abs(xdev[4 * i + 3] - 0.33) > tol)
            {
                std::cout << "Device: x=" << xdev[4 * i] << ","
                          << xdev[4 * i + 1] << "," << xdev[4 * i + 2] << ","
                          << xdev[4 * i + 3] << std::endl;
                std::cout << "Difference: " << xdev[4 * i] - 0.33 << ", "
                          << xdev[4 * i + 1] - 0.33 << ", "
                          << xdev[4 * i + 2] - 0.33 << ", "
                          << xdev[4 * i + 3] - 0.33 << std::endl;
                std::cout << "nits[" << i << "]=" << nits[i] << std::endl;
                count++;
            }
            if (count > 20) break;
        }
    }

    delete[] xdev;
    delete[] nits;
}
