#include "CALPHADFunctions.h"
#include "CALPHADConcSolverTernary.h"
#include "PhysicalConstants.h"
#include "CALPHADSpeciesPhaseGibbsEnergy.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"

#include <chrono>

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <iostream>

#include <omp.h>

namespace pt = boost::property_tree;

typedef std::chrono::high_resolution_clock Clock;


int main(int argc, char *argv[])
{
    const int N = 10000000;

#ifdef _OPENMP
    printf("Compiled by an OpenMP-compliant implementation.\n");
# endif

    std::cout << "Run test with " << omp_get_max_threads() << " threads"
              << std::endl;

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

    double LmixABPhaseL[4][2];
    double LmixABPhaseA[4][2];

    double LmixACPhaseL[4][2];
    double LmixACPhaseA[4][2];

    double LmixBCPhaseL[4][2];
    double LmixBCPhaseA[4][2];

    double LmixABCPhaseL[3][2];
    double LmixABCPhaseA[3][2];

    {
        std::string dbnamemixL("LmixABCPhaseL");
        if (calphad_db.get_child_optional(dbnamemixL))
        {
            pt::ptree& Lmix0_db = calphad_db.get_child(dbnamemixL);
            Thermo4PFM::readLmixTernaryParameters(Lmix0_db, LmixABCPhaseL);
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

    double fA[2];
    fA[0]     = g_species_phaseL[0].fenergy(temperature);
    fA[1]     = g_species_phaseA[0].fenergy(temperature);
    //std::cout<<"fA[0]="<<fA[0]<<", fA[1]="<<fA[1]<<std::endl;

    double fB[2];
    fB[0]     = g_species_phaseL[1].fenergy(temperature);
    fB[1]     = g_species_phaseA[1].fenergy(temperature);
    //std::cout<<"fB[0]="<<fB[0]<<", fB[1]="<<fB[1]<<std::endl;

    double fC[2];
    fC[0]     = g_species_phaseL[2].fenergy(temperature);
    fC[1]     = g_species_phaseA[2].fenergy(temperature);

    double L_AB_L[4];
    for(int i=0;i<4;i++)L_AB_L[i]=LmixABPhaseL[i][0]+temperature*LmixABPhaseL[i][1];

    double L_AB_S[4];
    for(int i=0;i<4;i++)L_AB_S[i]=LmixABPhaseA[i][0]+temperature*LmixABPhaseA[i][1];

    double L_AC_L[4];
    for(int i=0;i<4;i++)L_AC_L[i]=LmixACPhaseL[i][0]+temperature*LmixACPhaseL[i][1];

    double L_AC_S[4];
    for(int i=0;i<4;i++)L_AC_S[i]=LmixACPhaseA[i][0]+temperature*LmixACPhaseA[i][1];

    double L_BC_L[4];
    for(int i=0;i<4;i++)L_BC_L[i]=LmixBCPhaseL[i][0]+temperature*LmixBCPhaseL[i][1];

    double L_BC_S[4];
    for(int i=0;i<4;i++)L_BC_S[i]=LmixBCPhaseA[i][0]+temperature*LmixBCPhaseA[i][1];

    double L_ABC_L[3];
    for(int i=0;i<3;i++)L_ABC_L[i]=LmixABCPhaseL[i][0]+temperature*LmixABCPhaseL[i][1];

    double L_ABC_S[3];
    for(int i=0;i<3;i++)L_ABC_S[i]=LmixABCPhaseA[i][0]+temperature*LmixABCPhaseA[i][1];


    const double RTinv = 1.0 / (Thermo4PFM::gas_constant_R_JpKpmol * temperature);

    double sol[4]={ 0.33, 0.38, 0.32, 0.33 };

    double deviation = 0.01/(double)N;

    double xhost[4*N];
    for(int i=0;i<4*N;i++)
    {
        xhost[i]=-1.;
    }


// Host solve
{
    short nits[N];
    auto t1 = Clock::now();

#pragma omp parallel for
for(int i=0;i<N;i++)
{
    if(! omp_is_initial_device() ) abort();

    xhost[4*i]=sol[0];
    xhost[4*i+1]=sol[1];
    xhost[4*i+2]=sol[2];
    xhost[4*i+3]=sol[3];
    double hphi = 0.5+i*deviation;
    double c0 = 0.33;
    double c1 = 0.33;
    Thermo4PFM::CALPHADConcSolverTernary solver;
    solver.setup(c0, c1, hphi, RTinv, L_AB_L, L_AC_L, L_BC_L,
        L_AB_S, L_AC_S, L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);
    nits[i] = solver.ComputeConcentration(&xhost[4*i], 1.e-8, 50);
}
    auto t2 = Clock::now();
    long int usec=
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout<<"Host time/us/solve:   "<<(double)usec/(double)N<<std::endl;


if(N<20)
for(int i=0;i<N;i++)
{
    std::cout<<"Host: x="<<xhost[4*i]<<","<<xhost[4*i+1]<<","
                         <<xhost[4*i+2]<<","<<xhost[4*i+3]<<std::endl;
    std::cout<<"nits="<<nits[i]<<std::endl;
}

}

// Device solve
{
    double xdev[4*N];
    for(int i=0;i<4*N;i++)
    {
        xdev[i]=-1;
    }

    short nits[N];

    auto t1 = Clock::now();

# pragma omp target \
    map (to: sol ) map ( tofrom: xdev ) \
    map (to: fA, fB, fC) \
    map (to: L_AB_L, L_AC_L, L_BC_L, L_AB_S, L_AC_S, L_BC_S, L_ABC_L, L_ABC_S) \
    map (to: RTinv), map ( from: nits)
{
#pragma omp teams distribute parallel for
for(int i=0;i<N;i++)
{
    //if( omp_is_initial_device() ) abort();
    xdev[4*i]=sol[0];
    xdev[4*i+1]=sol[1];
    xdev[4*i+2]=sol[2];
    xdev[4*i+3]=sol[3];

    double hphi = 0.5+i*deviation;
    double c0 = 0.33;
    double c1 = 0.33;
    class Thermo4PFM::CALPHADConcSolverTernary solver;
    solver.setup(c0, c1, hphi, RTinv, L_AB_L, L_AC_L, L_BC_L,
        L_AB_S, L_AC_S, L_BC_S, L_ABC_L, L_ABC_S, fA, fB, fC);
    nits[i] = solver.ComputeConcentration(&xdev[4*i], 1.e-8, 50);
}
}

    auto t2 = Clock::now();
    long int usec=
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout<<"Device time/us/solve: "<<(double)usec/(double)N<<std::endl;

    double tol = 1.e-8;
    for(int i=0;i<N;i++)
    {
        if( std::abs(xdev[4*i]-xhost[4*i])>tol ||
          std::abs(xdev[4*i+1]-xhost[4*i+1])>tol ||
          std::abs(xdev[4*i+2]-xhost[4*i+2])>tol ||
          std::abs(xdev[4*i+3]-xhost[4*i+3])>tol ||
          N<20)
        {
        std::cout<<"Device: x="<<xdev[4*i]<<","<<xdev[4*i+1]<<","
                               <<xdev[4*i+2]<<","<<xdev[4*i+3]<<std::endl;
        std::cout<<"Difference: "<<xdev[4*i]-xhost[4*i]<<", "
                                 <<xdev[4*i+1]-xhost[4*i+1]<<", "
                                 <<xdev[4*i+2]-xhost[4*i+2]<<", "
                                 <<xdev[4*i+3]-xhost[4*i+3]
                                 <<std::endl;
        std::cout<<"nits["<<i<<"]="<<nits[i]<<std::endl;
        }
    }
}

}
