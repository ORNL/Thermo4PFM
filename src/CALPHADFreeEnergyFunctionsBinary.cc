// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADEqConcSolverBinary.h"
#include "CALPHADFunctions.h"
#include "PhysicalConstants.h"
//#include "FuncFort.h"
#include "functions.h"

#include <boost/property_tree/json_parser.hpp>

#include <iomanip>
#include <string>

namespace pt = boost::property_tree;

#ifdef HAVE_TLOT
void readLcoefficients(pt::ptree& db, double (&LmixPhase)[4][3])
{
#else
void readLcoefficients(pt::ptree& db, double (&LmixPhase)[4][2])
{
#endif
    {
        int i = 0;
        for (pt::ptree::value_type& v : db.get_child("L0"))
        {
            LmixPhase[0][i] = v.second.get_value<double>();
            i++;
        }
#ifdef HAVE_TLOGT
        if (i < 3) LmixPhase[0][2] = 0.;
#endif
    }

    {
        int i = 0;
        for (pt::ptree::value_type& v : db.get_child("L1"))
        {
            LmixPhase[1][i] = v.second.get_value<double>();
            i++;
        }
#ifdef HAVE_TLOGT
        if (i < 3) LmixPhase[1][2] = 0.;
#endif
    }

    // L2
    {
        auto child = db.get_child_optional("L2");
        if (child)
        {
            int i = 0;
            for (pt::ptree::value_type& v : db.get_child("L2"))
            {
                LmixPhase[2][i] = v.second.get_value<double>();
                i++;
            }
#ifdef HAVE_TLOGT
            if (i < 3) LmixPhase[2][2] = 0.;
#endif
        }
        else
        {
            LmixPhase[2][0] = 0.0;
            LmixPhase[2][1] = 0.0;
#ifdef HAVE_TLOGT
            LmixPhase[2][2] = 0.0;
#endif
        }
    }

    // L3
    {
        auto child = db.get_child_optional("L3");
        if (child)
        {
            int i = 0;
            for (pt::ptree::value_type& v : db.get_child("L3"))
            {
                LmixPhase[3][i] = v.second.get_value<double>();
                i++;
            }
#ifdef HAVE_TLOGT
            if (i < 3) LmixPhase[3][2] = 0.;
#endif
        }
        else
        {
            LmixPhase[3][0] = 0.0;
            LmixPhase[3][1] = 0.0;
#ifdef HAVE_TLOGT
            LmixPhase[3][2] = 0.0;
#endif
        }
    }
}

CALPHADFreeEnergyFunctionsBinary::CALPHADFreeEnergyFunctionsBinary(
    pt::ptree& calphad_db, boost::optional<pt::ptree&> newton_db,
    const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type,
    const bool with_third_phase)
    : d_energy_interp_func_type(energy_interp_func_type),
      d_conc_interp_func_type(conc_interp_func_type),
      d_with_third_phase(with_third_phase)
{
    d_fenergy_diag_filename = "energy.vtk";

    int N = 2;
    if (d_with_third_phase)
    {
        N = 3;
    }

    d_fA = new double[N];
    d_fB = new double[N];
    d_L0 = new double[N];
    d_L1 = new double[N];
    d_L2 = new double[N];
    d_L3 = new double[N];

    d_ceq_l = -1;
    d_ceq_a = -1;
    d_ceq_b = -1;

    readParameters(calphad_db);

    setupSolver(newton_db);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::setupSolver(
    boost::optional<pt::ptree&> newton_db)
{
    std::cout << "CALPHADFreeEnergyFunctionsBinary::setupSolver()..."
              << std::endl;
    d_solver = new CALPHADConcentrationSolverBinary(d_with_third_phase);

    if (newton_db) readNewtonparameters(newton_db.get());
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::readNewtonparameters(
    pt::ptree& newton_db)
{
    double tol         = newton_db.get<double>("tol", 1.e-8);
    double alpha       = newton_db.get<double>("alpha", 1.);
    int maxits         = newton_db.get<int>("max_its", 20);
    const bool verbose = newton_db.get<bool>("verbose", false);

    d_solver->SetTolerance(tol);
    d_solver->SetMaxIterations(maxits);
    d_solver->SetDamping(alpha);
    d_solver->SetVerbose(verbose);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::readParameters(pt::ptree& calphad_db)
{
    pt::ptree& species0_db = calphad_db.get_child("SpeciesA");
    std::string dbnameL("PhaseL");
    d_g_species_phaseL[0].initialize("L0", species0_db.get_child(dbnameL));
    std::string dbnameA("PhaseA");
    d_g_species_phaseA[0].initialize("A0", species0_db.get_child(dbnameA));
    std::string dbnameB("PhaseB");
    if (d_with_third_phase)
    {
        d_g_species_phaseB[0].initialize("B0", species0_db.get_child(dbnameB));
    }

    pt::ptree& speciesB_db = calphad_db.get_child("SpeciesB");
    d_g_species_phaseL[1].initialize("L1", speciesB_db.get_child(dbnameL));
    d_g_species_phaseA[1].initialize("A1", speciesB_db.get_child(dbnameA));
    if (d_with_third_phase)
    {
        d_g_species_phaseB[1].initialize("B1", speciesB_db.get_child(dbnameB));
    }

    // read Lmix coefficients
    std::string dbnamemixL("LmixPhaseL");
    std::string dbnamemixA("LmixPhaseA");
    std::string dbnamemixB("LmixPhaseB");
    pt::ptree Lmix0_db = calphad_db.get_child(dbnamemixL);
    readLcoefficients(Lmix0_db, d_LmixPhaseL);

    pt::ptree Lmix1_db = calphad_db.get_child(dbnamemixA);
    readLcoefficients(Lmix1_db, d_LmixPhaseA);

    if (d_with_third_phase)
    {
        pt::ptree Lmix2_db = calphad_db.get_child(dbnamemixB);
        readLcoefficients(Lmix2_db, d_LmixPhaseB);
    }

    // print database just read
    std::clog << "CALPHAD database..." << std::endl;
    pt::write_json(std::clog, calphad_db);
}

//-----------------------------------------------------------------------

double CALPHADFreeEnergyFunctionsBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
    const double l0 = lmixPhase(0, pi, temperature);
    const double l1 = lmixPhase(1, pi, temperature);
    const double l2 = lmixPhase(2, pi, temperature);
    const double l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &d_g_species_phaseL[0];
            break;
        case PhaseIndex::phaseA:
            g_species = &d_g_species_phaseA[0];
            break;
        case PhaseIndex::phaseB:
            g_species = &d_g_species_phaseB[0];
            break;
        default:
            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
                         "computeFreeEnergy(), undefined phase"
                      << "!!!" << std::endl;
            abort();
            return 0.;
    }

    double fe = conc[0] * g_species[0].fenergy(temperature)
                + (1. - conc[0]) * g_species[1].fenergy(temperature)
                + CALPHADcomputeFMixBinary(l0, l1, l2, l3, conc[0])
                + CALPHADcomputeFIdealMixBinary(
                      gas_constant_R_JpKpmol * temperature, conc[0]);

    // subtract -mu*c to get grand potential
    if (gp)
    {
        double deriv;
        computeDerivFreeEnergy(temperature, conc, pi, &deriv);
        fe -= deriv * conc[0];
    }

    return fe;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
    const double l0 = lmixPhase(0, pi, temperature);
    const double l1 = lmixPhase(1, pi, temperature);
    const double l2 = lmixPhase(2, pi, temperature);
    const double l3 = lmixPhase(3, pi, temperature);

    CALPHADSpeciesPhaseGibbsEnergy* g_species;

    switch (pi)
    {
        case PhaseIndex::phaseL:
            g_species = &d_g_species_phaseL[0];
            break;
        case PhaseIndex::phaseA:
            g_species = &d_g_species_phaseA[0];
            break;
        case PhaseIndex::phaseB:
            g_species = &d_g_species_phaseB[0];
            break;
        default:
            std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
                         "computeFreeEnergy(), undefined phase!!!"
                      << std::endl;
            abort();
            return;
    }

    double mu = (g_species[0].fenergy(temperature)
                    - g_species[1].fenergy(temperature))
                + CALPHADcomputeFMix_derivBinary(l0, l1, l2, l3, conc[0])
                + CALPHADcomputeFIdealMix_derivBinary(
                      gas_constant_R_JpKpmol * temperature, conc[0]);

    deriv[0] = mu;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computeSecondDerivativeFreeEnergy(
    const double temp, const double* const conc, const PhaseIndex pi,
    std::vector<double>& d2fdc2)
{
    assert(conc[0] >= 0.);
    assert(conc[0] <= 1.);

    const double l0_l = lmixPhase(0, pi, temp);
    const double l1_l = lmixPhase(1, pi, temp);
    const double l2_l = lmixPhase(2, pi, temp);
    const double l3_l = lmixPhase(3, pi, temp);
    const double rt   = gas_constant_R_JpKpmol * temp;

    d2fdc2[0]
        = (CALPHADcomputeFMix_deriv2Binary(l0_l, l1_l, l2_l, l3_l, conc[0])
            + CALPHADcomputeFIdealMix_deriv2Binary(rt, conc[0]));
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::setupValuesForTwoPhasesSolver(
    const double temperature, double* L0, double* L1, double* L2, double* L3,
    double* fA, double* fB, const PhaseIndex pi0, const PhaseIndex pi1)
{
    PhaseIndex pis[2] = { pi0, pi1 };

    // loop over two phases
    for (short i = 0; i < 2; i++)
    {
        L0[i] = lmixPhase(0, pis[i], temperature);
        L1[i] = lmixPhase(1, pis[i], temperature);
        L2[i] = lmixPhase(2, pis[i], temperature);
        L3[i] = lmixPhase(3, pis[i], temperature);

        switch (pis[i])
        {

            case PhaseIndex::phaseL:
                fA[i] = d_g_species_phaseL[0].fenergy(temperature);
                fB[i] = d_g_species_phaseL[1].fenergy(temperature);
                break;

            case PhaseIndex::phaseA:
                fA[i] = d_g_species_phaseA[0].fenergy(temperature);
                fB[i] = d_g_species_phaseA[1].fenergy(temperature);
                break;

            case PhaseIndex::phaseB:
                fA[i] = d_g_species_phaseB[0].fenergy(temperature);
                fB[i] = d_g_species_phaseB[1].fenergy(temperature);
                break;

            default:
                std::cerr << "CALPHADFreeEnergyFunctionsBinary::"
                             "setupValuesForTwoPhasesSolver: Undefined phase"
                          << std::endl;
        }
    }
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::setup(const double temperature)
{
    d_fA[0] = d_g_species_phaseL[0].fenergy(temperature);
    d_fA[1] = d_g_species_phaseA[0].fenergy(temperature);

    d_fB[0] = d_g_species_phaseL[1].fenergy(temperature);
    d_fB[1] = d_g_species_phaseA[1].fenergy(temperature);

    d_L0[0] = lmixPhase(0, PhaseIndex::phaseL, temperature);
    d_L1[0] = lmixPhase(1, PhaseIndex::phaseL, temperature);
    d_L2[0] = lmixPhase(2, PhaseIndex::phaseL, temperature);
    d_L3[0] = lmixPhase(3, PhaseIndex::phaseL, temperature);

    d_L0[1] = lmixPhase(0, PhaseIndex::phaseA, temperature);
    d_L1[1] = lmixPhase(1, PhaseIndex::phaseA, temperature);
    d_L2[1] = lmixPhase(2, PhaseIndex::phaseA, temperature);
    d_L3[1] = lmixPhase(3, PhaseIndex::phaseA, temperature);

    if (d_with_third_phase)
    {
        d_fA[2] = d_g_species_phaseB[0].fenergy(temperature);
        d_fB[2] = d_g_species_phaseB[1].fenergy(temperature);

        d_L0[2] = lmixPhase(0, PhaseIndex::phaseB, temperature);
        d_L1[2] = lmixPhase(1, PhaseIndex::phaseB, temperature);
        d_L2[2] = lmixPhase(2, PhaseIndex::phaseB, temperature);
        d_L3[2] = lmixPhase(3, PhaseIndex::phaseB, temperature);
    }
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyFunctionsBinary::computeCeqT(const double temperature,
    const PhaseIndex pi0, const PhaseIndex pi1, double* ceq, const int maxits,
    const bool verbose)
{
    if (verbose)
        std::cout << "CALPHADFreeEnergyFunctionsBinary::computeCeqT()"
                  << std::endl;
    assert(temperature > 0.);

    setupValuesForTwoPhasesSolver(
        temperature, d_L0, d_L1, d_L2, d_L3, d_fA, d_fB, pi0, pi1);
    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);
    CALPHADEqConcentrationSolverBinary eq_solver;
    eq_solver.SetMaxIterations(maxits);

    int ret = eq_solver.ComputeConcentration(
        ceq, RTinv, d_L0, d_L1, d_L2, d_L3, d_fA, d_fB);

    if (ret >= 0)
    {
        if (verbose)
        {
            std::cout << "CALPHAD, c_eq phase0=" << ceq[0] << std::endl;
            std::cout << "CALPHAD, c_eq phase1=" << ceq[1] << std::endl;
        }

        if (pi1 == PhaseIndex::phaseB)
        {
            d_ceq_b = ceq[1];
            d_ceq_l = 0.5 * (ceq[0] + d_ceq_l);
        }
        else
        {
            d_ceq_l = ceq[0];
            d_ceq_a = ceq[1];
        }
    }
    else
    {
        std::cout << "CALPHADFreeEnergyFunctionsBinary, WARNING: ceq "
                     "computation did not converge"
                  << std::endl;
    }

    return (ret >= 0);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies(
    const double temperature, const double hphi, const double heta,
    const double conc, double& fl, double& fa, double& fb)
{
    // std::cout<<"CALPHADFreeEnergyFunctionsBinary::computePhasesFreeEnergies()"<<endl;

    const int N = d_with_third_phase ? 3 : 2;

    double* c = new double[N];
    for (int ii = 0; ii < N; ii++)
    {
        c[ii] = conc;
    }
    // std::cout<<"d_ceq_l="<<d_ceq_l<<endl;
    // std::cout<<"d_ceq_a="<<d_ceq_a<<endl;
    // std::cout<<"d_ceq_b="<<d_ceq_b<<endl;
    if (d_ceq_l >= 0.) c[0] = d_ceq_l;
    if (d_ceq_a >= 0.) c[1] = d_ceq_a;
    if (d_ceq_b >= 0.) c[2] = d_ceq_b;

    setup(temperature);

    double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);
    int ret      = d_solver->ComputeConcentration(
        c, conc, hphi, heta, RTinv, d_L0, d_L1, d_L2, d_L3, d_fA, d_fB);

    if (ret < 0)
    {
        if (d_with_third_phase)
        {
            std::cerr << "d_ceq_l=" << d_ceq_l << std::endl;
            std::cerr << "d_ceq_a=" << d_ceq_a << std::endl;
            std::cerr << "d_ceq_b=" << d_ceq_b << std::endl;
            std::cerr << "ERROR in "
                         "CALPHADFreeEnergyFunctionsBinary::"
                         "computePhasesFreeEnergies()"
                         " ---"
                      << "conc=" << conc << ", hphi=" << hphi
                      << ", heta=" << heta << std::endl;
        }
        else
        {
            std::cerr << "ERROR in "
                         "CALPHADFreeEnergyFunctionsBinary::"
                         "computePhasesFreeEnergies()"
                         " ---"
                      << "conc=" << conc << ", hphi=" << hphi << std::endl;
        }
        abort();
    }

    assert(c[0] >= 0.);
    fl = computeFreeEnergy(temperature, &c[0], PhaseIndex::phaseL, false);

    assert(c[1] >= 0.);
    fa = computeFreeEnergy(temperature, &c[1], PhaseIndex::phaseA, false);

    fb = 0.;
    if (d_with_third_phase)
    {
        assert(c[2] >= 0.);
        assert(c[2] <= 1.);
        fb = computeFreeEnergy(temperature, &c[2], PhaseIndex::phaseB, false);
    }
    delete[] c;
}

//-----------------------------------------------------------------------

int CALPHADFreeEnergyFunctionsBinary::computePhaseConcentrations(
    const double temperature, const double* const conc, const double phi,
    const double eta, double* x)
{
    assert(x[0] >= 0.);
    assert(x[1] >= 0.);
    assert(x[0] <= 1.);
    assert(x[1] <= 1.);

    const double conc0 = conc[0];

    const double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);

    d_fA[0] = getFenergyPhaseL(0, temperature);
    d_fA[1] = getFenergyPhaseA(0, temperature);

    d_fB[0] = getFenergyPhaseL(1, temperature);
    d_fB[1] = getFenergyPhaseA(1, temperature);

    d_L0[0] = lmixPhase(0, PhaseIndex::phaseL, temperature);
    d_L1[0] = lmixPhase(1, PhaseIndex::phaseL, temperature);
    d_L2[0] = lmixPhase(2, PhaseIndex::phaseL, temperature);
    d_L3[0] = lmixPhase(3, PhaseIndex::phaseL, temperature);

    d_L0[1] = lmixPhase(0, PhaseIndex::phaseA, temperature);
    d_L1[1] = lmixPhase(1, PhaseIndex::phaseA, temperature);
    d_L2[1] = lmixPhase(2, PhaseIndex::phaseA, temperature);
    d_L3[1] = lmixPhase(3, PhaseIndex::phaseA, temperature);

    const char interp_func_type = concInterpChar(d_conc_interp_func_type);
    const double hphi           = interp_func(phi, interp_func_type);

    double heta = 0.0;
    // std::cout<<"d_ceq_a="<<d_ceq_a<<endl;
    // x[0] = ( d_ceq_l>=0. ) ? d_ceq_l : 0.5;
    // x[1] = ( d_ceq_a>=0. ) ? d_ceq_a : 0.5;

    if (d_with_third_phase)
    {
        assert(x[2] >= 0.);
        assert(x[2] <= 1.);
        // x[2] = ( d_ceq_b>=0. ) ? d_ceq_b : 0.5;

        d_fA[2] = getFenergyPhaseB(0, temperature);
        d_fB[2] = getFenergyPhaseB(1, temperature);

        d_L0[2] = lmixPhase(0, PhaseIndex::phaseB, temperature);
        d_L1[2] = lmixPhase(1, PhaseIndex::phaseB, temperature);
        d_L2[2] = lmixPhase(2, PhaseIndex::phaseB, temperature);
        d_L3[2] = lmixPhase(3, PhaseIndex::phaseB, temperature);

        heta = interp_func(eta, interp_func_type);
    }

    // conc could be outside of [0.,1.] in a trial step
    double c0 = conc[0] >= 0. ? conc[0] : 0.;
    c0        = c0 <= 1. ? c0 : 1.;
    int ret   = d_solver->ComputeConcentration(
        x, c0, hphi, heta, RTinv, d_L0, d_L1, d_L2, d_L3, d_fA, d_fB);
    if (ret == -1)
    {
        std::cerr << "ERROR, "
                     "CALPHADFreeEnergyFunctionsBinary::"
                     "computePhaseConcentrations() "
                     "failed for conc="
                  << conc0 << ", hphi=" << hphi << ", heta=" << heta
                  << std::endl;
        abort();
    }

    return ret;
}

//-----------------------------------------------------------------------

void CALPHADFreeEnergyFunctionsBinary::energyVsPhiAndC(const double temperature,
    const double* const ceq, const bool found_ceq, const double phi_well_scale,
    const std::string& phi_well_type, const int npts_phi, const int npts_c)
{
    std::cout << "CALPHADFreeEnergyFunctionsBinary::energyVsPhiAndC()..."
              << std::endl;

    double slopec = 0.;
    double fc0    = 0.;
    double fc1    = 0.;
    if (found_ceq)
    {
        // compute slope of f between equilibrium concentrations
        // to add slopec*conc to energy later on

        fc0    = computeFreeEnergy(temperature, &ceq[0], PhaseIndex::phaseL);
        fc1    = computeFreeEnergy(temperature, &ceq[1], PhaseIndex::phaseA);
        slopec = -(fc1 - fc0) / (ceq[1] - ceq[0]);
    }
    std::cout << std::setprecision(8) << "fc0: " << fc0 << "..."
              << ", fc1: " << fc1 << "..." << std::endl;
    std::cout << "CALPHADFreeEnergyFunctionsBinary: Use slope: " << slopec
              << "..." << std::endl;

    // reset cmin, cmax, deltac
    double cmin   = std::min(ceq[0], ceq[1]);
    double cmax   = std::max(ceq[0], ceq[1]);
    double dc     = cmax - cmin;
    cmin          = std::max(0.25 * cmin, cmin - 0.25 * dc);
    cmax          = std::min(1. - 0.25 * (1. - cmax), cmax + 0.25 * dc);
    cmax          = std::max(cmax, cmin + dc);
    double deltac = (cmax - cmin) / (npts_c - 1);

    std::ofstream tfile(d_fenergy_diag_filename.data(), std::ios::out);

    printEnergyVsPhiHeader(
        temperature, npts_phi, npts_c, cmin, cmax, slopec, tfile);

    for (int i = 0; i < npts_c; i++)
    {
        double conc = cmin + deltac * i;
        printEnergyVsPhi(&conc, temperature, phi_well_scale, phi_well_type,
            npts_phi, slopec, tfile);
    }
}

// Print out free energy as a function of phase
// for given composition and temperature
// File format: ASCII VTK, readble with Visit
void CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhiHeader(
    const double temperature, const int nphi, const int nc, const double cmin,
    const double cmax, const double slopec, std::ostream& os) const
{
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << "Free energy + " << slopec << "*c [J/mol] at T=" << temperature
       << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET STRUCTURED_POINTS" << std::endl;

    os << "DIMENSIONS   " << nphi << " " << nc << " 1" << std::endl;
    double asp_ratio_c = (nc > 1) ? (cmax - cmin) / (nc - 1) : 1.;
    os << "ASPECT_RATIO " << 1. / (nphi - 1) << " " << asp_ratio_c << " 1."
       << std::endl;
    os << "ORIGIN        0. " << cmin << " 0." << std::endl;
    os << "POINT_DATA   " << nphi * nc << std::endl;
    os << "SCALARS energy float 1" << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhi(
    const double* const conc, const double temperature,
    const double phi_well_scale, const std::string& phi_well_type,
    const int npts, const double slopec, std::ostream& os)
{
    // std::cout << "CALPHADFreeEnergyFunctionsBinary::printEnergyVsPhi()..." <<
    // std::endl;
    const double dphi = 1.0 / (double)(npts - 1);
    const double eta  = 0.0;

    // os << "# phi     f(phi)     for c=" << conc
    //           << " eta=" << eta
    //           << " and T=" << temperature << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double phi = i * dphi;

        double e       = fchem(phi, eta, conc, temperature);
        const double w = phi_well_scale * well_func(phi);

        os << e + w + slopec * conc[0] << std::endl;
    }
    // os << std::endl;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsEta(
    const double* const conc, const double temperature,
    const double eta_well_scale, const std::string& eta_well_type,
    const int npts, const double slopec, std::ostream& os)
{
    // std::cout << "CALPHADFreeEnergyFunctionsBinary::printEnergyVsEta()..." <<
    // std::endl;
    const double deta = 1.0 / (double)(npts - 1);
    const double phi  = 1.0;

    for (int i = 0; i < npts; i++)
    {
        const double eta = i * deta;

        double e = fchem(phi, eta, conc, temperature);

        double w = eta_well_scale * well_func(eta);

        os << e + w + slopec * conc[0] << std::endl;
    }
    // os << std::endl;
}

//=======================================================================
// compute free energy in [J/mol]
double CALPHADFreeEnergyFunctionsBinary::fchem(const double phi,
    const double eta, const double* const conc, const double temperature)
{
    const char interp_func_type = concInterpChar(d_conc_interp_func_type);
    const double hcphi          = interp_func(phi, interp_func_type);
    double heta                 = 0.0;
    if (d_with_third_phase)
    {
        heta = interp_func(eta, interp_func_type);
    }

    const double tol = 1.e-8;
    double fl        = 0.;
    double fa        = 0.;
    double fb        = 0.;
    if ((phi > tol) & (phi < (1. - tol)))
    {
        computePhasesFreeEnergies(
            temperature, hcphi, heta, conc[0], fl, fa, fb);
    }
    else
    {
        if (phi <= tol)
        {
            fl = computeFreeEnergy(temperature, conc, PhaseIndex::phaseL);
        }
        else
        {
            fa = computeFreeEnergy(temperature, conc, PhaseIndex::phaseA);
            if (d_with_third_phase)
                fb = computeFreeEnergy(temperature, conc, PhaseIndex::phaseB);
        }
    }

    const char interpf = energyInterpChar(d_energy_interp_func_type);
    const double hfphi = interp_func(phi, interpf);
    double e = (1.0 - hfphi) * fl + hfphi * ((1.0 - heta) * fa + heta * fb);

    return e;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsBinary::printEnergyVsComposition(
    const double temperature, const int npts)
{
    std::ofstream os("FvsC.dat", std::ios::out);

    const double dc = 1.0 / (double)(npts - 1);

    os << "#phi=0" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        double e = fchem(0., 0., &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }
    os << std::endl << std::endl;

    os << "#phi=1" << std::endl;
    for (int i = 0; i < npts; i++)
    {
        const double conc = i * dc;

        double e = fchem(1., 0., &conc, temperature);
        os << conc << "\t" << e << std::endl;
    }

    if (d_with_third_phase)
    {
        os << std::endl;

        os << "#eta=1" << std::endl;
        for (int i = 0; i < npts; i++)
        {
            const double conc = i * dc;

            double e = fchem(1., 1., &conc, temperature);
            os << conc << "\t" << e << std::endl;
        }
    }
}
