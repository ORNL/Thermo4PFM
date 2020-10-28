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
#ifndef included_CALPHADEqPhaseConcSolverTernary
#define included_CALPHADEqPhaseConcSolverTernary

#include "DampedNewtonSolver.h"

class CALPHADEqPhaseConcentrationSolverTernary : public DampedNewtonSolver
{
public:
    CALPHADEqPhaseConcentrationSolverTernary(const double c0, const double c1);

    virtual ~CALPHADEqPhaseConcentrationSolverTernary(){};

    int ComputeConcentration(double* const conc, const double RTinv,
        const double* const L_AB_L, const double* const L_AC_L,
        const double* const L_BC_L, const double* const L_AB_S,
        const double* const L_AC_S, const double* const L_BC_S,
        const double* const L_ABC_L, const double* const L_ABC_S,
        const double* const fA, const double* const fB, const double* const fC);

    int ComputeConcentration(double* const conc);

    void setup(const double RTinv, const double* const L_AB_L,
        const double* const L_AC_L, const double* const L_BC_L,
        const double* const L_AB_S, const double* const L_AC_S,
        const double* const L_BC_S, const double* const L_ABC_L,
        const double* const L_ABC_S, const double* const fA,
        const double* const fB, const double* const fC);

    virtual void RHS(const double* const x, double* const fvec);

    virtual void Jacobian(const double* const x, double** const fjac);

private:
    double RTinv_;
    double RT_;
    double hphi_;

    // energies of 3 species, in two phase each
    double fA_[2];
    double fB_[2];
    double fC_[2];

    // L coefficients for 2 possible phases (L and S)
    double L_AB_L_[4];
    double L_AC_L_[4];
    double L_BC_L_[4];
    double L_ABC_L_[3];

    double L_AB_S_[4];
    double L_AC_S_[4];
    double L_BC_S_[4];
    double L_ABC_S_[3];

    // nominal concentration
    double conc_[2];
};

#endif
