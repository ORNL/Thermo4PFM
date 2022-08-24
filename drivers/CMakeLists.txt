find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/src )

add_executable(loopxlogx
               ${CMAKE_SOURCE_DIR}/drivers/loopxlogx.cc)

add_executable(loopCALPHADSpeciesPhaseGibbsEnergy
               ${CMAKE_SOURCE_DIR}/drivers/loopCALPHADSpeciesPhaseGibbsEnergy.cc)

add_executable(loopDeterminant
               ${CMAKE_SOURCE_DIR}/drivers/loopDeterminant.cc
               ${CMAKE_SOURCE_DIR}/src/Determinant.cc )

add_executable(loopCALPHADConcSolverBinary
               ${CMAKE_SOURCE_DIR}/drivers/loopCALPHADConcSolverBinary.cc)

add_executable(loopCALPHADConcSolverTernary
               ${CMAKE_SOURCE_DIR}/drivers/loopCALPHADConcSolverTernary.cc )

add_executable(loopCALPHADbinaryKKS
               ${CMAKE_SOURCE_DIR}/drivers/loopCALPHADbinaryKKS.cc)

add_executable(plotEnergyVsComposition
               ${CMAKE_SOURCE_DIR}/drivers/plotEnergyVsComposition.cc )
add_executable(plotEnergyVsCompositionThreePhase
               ${CMAKE_SOURCE_DIR}/drivers/plotEnergyVsCompositionThreePhase.cc )
add_executable(computeKKSTernary
               ${CMAKE_SOURCE_DIR}/drivers/computeKKSTernary.cc)
add_executable(computeTieLineTernary
               ${CMAKE_SOURCE_DIR}/drivers/computeTieLineTernary.cc)
add_executable(computeKKS3Ph2Sl
               ${CMAKE_SOURCE_DIR}/drivers/computeKKS3Ph2Sl.cc)

target_link_libraries(plotEnergyVsComposition thermo4pfm)
target_link_libraries(plotEnergyVsCompositionThreePhase thermo4pfm)
target_link_libraries(computeKKSTernary thermo4pfm)
target_link_libraries(computeTieLineTernary thermo4pfm)
target_link_libraries(computeKKS3Ph2Sl thermo4pfm)
target_link_libraries(loopxlogx thermo4pfm)
target_link_libraries(loopCALPHADSpeciesPhaseGibbsEnergy thermo4pfm)
target_link_libraries(loopCALPHADConcSolverBinary thermo4pfm)
target_link_libraries(loopCALPHADConcSolverTernary thermo4pfm)
target_link_libraries(loopCALPHADbinaryKKS thermo4pfm)

target_link_libraries(loopxlogx m)
if (${OPENMP_CXX_FOUND})
  target_link_libraries(plotEnergyVsComposition OpenMP::OpenMP_CXX)
  target_link_libraries(loopxlogx OpenMP::OpenMP_CXX)
  target_link_libraries(loopDeterminant OpenMP::OpenMP_CXX)
  target_link_libraries(loopCALPHADSpeciesPhaseGibbsEnergy OpenMP::OpenMP_CXX)
  target_link_libraries(loopCALPHADConcSolverBinary OpenMP::OpenMP_CXX)
  target_link_libraries(loopCALPHADConcSolverTernary OpenMP::OpenMP_CXX)
  target_link_libraries(loopCALPHADbinaryKKS OpenMP::OpenMP_CXX)
endif()
