set(SOURCES ${CMAKE_SOURCE_DIR}/src/InterpolationType.cc
            ${CMAKE_SOURCE_DIR}/src/math_utilities.cc
            ${CMAKE_SOURCE_DIR}/src/NewtonSolver.cc
            ${CMAKE_SOURCE_DIR}/src/DampedNewtonSolver.cc
            ${CMAKE_SOURCE_DIR}/src/xlogx.cc
            ${CMAKE_SOURCE_DIR}/src/functions.cc
            ${CMAKE_SOURCE_DIR}/src/well_functions.cc
            ${CMAKE_SOURCE_DIR}/src/InterpolationType.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADFunctions.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADSpeciesPhaseGibbsEnergy.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADSpeciesPhaseGibbsEnergyExpansion.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADFreeEnergyFunctionsBinary.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADFreeEnergyFunctionsTernary.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADConcSolverBinary.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADConcSolverTernary.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADEqConcSolverBinary.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADEqConcSolverTernary.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADConcSolverTernary.cc
            ${CMAKE_SOURCE_DIR}/src/CALPHADEqPhaseConcSolverTernary.cc
            ${CMAKE_SOURCE_DIR}/src/KKSdiluteBinaryConcentrationSolver.cc
            ${CMAKE_SOURCE_DIR}/src/KKSFreeEnergyFunctionDiluteBinary.cc
            ${CMAKE_SOURCE_DIR}/src/xlogx.cc)

add_library(thermo4pfm ${SOURCES})

set(INSTALL_DIR ${CMAKE_SOURCE_DIR})

set(PUBLIC_HEADERS
    PhysicalConstants.h
    KKSFreeEnergyFunctionDiluteBinary.h
    InterpolationType.h
    math_utilities.h
    FreeEnergyFunctions.h
    Phases.h
    NewtonSolver.h
    DampedNewtonSolver.h
    xlogx.h
    functions.h
    well_functions.h
    CALPHADFunctions.h
    CALPHADSpeciesPhaseGibbsEnergy.h
    CALPHADSpeciesPhaseGibbsEnergyExpansion.h
    CALPHADFreeEnergyFunctions.h
    CALPHADFreeEnergyFunctionsBinary.h
    CALPHADFreeEnergyFunctionsTernary.h
    CALPHADConcSolverBinary.h
    CALPHADConcSolverTernary.h
    CALPHADEqConcSolverBinary.h
    CALPHADEqConcSolverTernary.h
    CALPHADEqPhaseConcSolverTernary.h
    KKSdiluteBinaryConcentrationSolver.h
   )
install(TARGETS thermo4pfm DESTINATION ${INSTALL_DIR}/lib)
install(FILES ${PUBLIC_HEADERS}
        DESTINATION ${INSTALL_DIR}/include)
