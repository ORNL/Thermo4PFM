find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/src )

set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

add_executable(testPropertyTree
      ${TEST_DIR}/property_tree/testPropertyTree.cc)
add_executable(testAlCuDatabases
             ${TEST_DIR}/testAlCuDatabases.cc )
add_executable(testInterpolationFunctions
               ${TEST_DIR}/testInterpolationFunctions.cc )
add_executable(testCALPHADFunctions
               ${TEST_DIR}/testCALPHADFunctions.cc )
add_executable(testCALPHADSpeciesPhaseGibbsEnergyExpansion
               ${TEST_DIR}/testCALPHADSpeciesPhaseGibbsEnergyExpansion.cc)
add_executable(testCALPHADGibbsEnergy
               ${TEST_DIR}/testCALPHADGibbsEnergy.cc )
add_executable(testCALPHADbinaryKKS
               ${TEST_DIR}/testCALPHADbinaryKKS.cc )
add_executable(testCALPHADbinaryKKSThreePhase
               ${TEST_DIR}/testCALPHADbinaryKKSThreePhase.cc )
add_executable(testCALPHADbinaryKKS3Ph2Sl
                ${TEST_DIR}/testCALPHADbinaryKKS3Ph2Sl.cc )
add_executable(testCALPHADbinaryEquilibrium
               ${TEST_DIR}/testCALPHADbinaryEquilibrium.cc )
add_executable(testCALPHADConcSolverBinaryThreePhase
              ${TEST_DIR}/testCALPHADConcSolverBinaryThreePhase.cc )
add_executable(testCALPHADConcSolverBinary3Ph2Sl
                ${TEST_DIR}/testCALPHADConcSolverBinary3Ph2Sl.cc )
add_executable(testKKSdiluteBinaryKKS
               ${TEST_DIR}/testKKSdiluteBinaryKKS.cc )
add_executable(testKKSdiluteBinaryEquilibrium
               ${TEST_DIR}/testKKSdiluteBinaryEquilibrium.cc )
add_executable(testCALPHADJacobianBinary
               ${TEST_DIR}/testCALPHADJacobianBinary.cc )
add_executable(testCALPHADJacobianTernary
               ${TEST_DIR}/testCALPHADJacobianTernary.cc )
add_executable(testCALPHADternaryKKS
               ${TEST_DIR}/testCALPHADternaryKKS.cc )
add_executable(testCALPHADConcSolverTernary
               ${TEST_DIR}/testCALPHADConcSolverTernary.cc )
add_executable(testCALPHADternaryTieLine
               ${TEST_DIR}/testCALPHADternaryTieLine.cc)

if(${WITH_OPENMP_OFFLOAD})
  add_executable(testOpenMPoffload
                 ${TEST_DIR}/testOpenMPoffload.cc )
endif()
add_executable(testLoopCALPHADSpeciesPhaseGibbsEnergy
               ${TEST_DIR}/testLoopCALPHADSpeciesPhaseGibbsEnergy.cc)
add_executable(testLoopCALPHADbinaryEquilibrium
               ${TEST_DIR}/testLoopCALPHADbinaryEquilibrium.cc)
add_executable(testLoopCALPHADbinaryKKS
               ${TEST_DIR}/testLoopCALPHADbinaryKKS.cc)
add_executable(testLoopCALPHADternaryKKS
               ${TEST_DIR}/testLoopCALPHADternaryKKS.cc)
add_executable(testLoopCALPHADternaryTieLine
               ${TEST_DIR}/testLoopCALPHADternaryTieLine.cc)
add_executable(testLoopKKSdiluteBinaryKKS
               ${TEST_DIR}/testLoopKKSdiluteBinaryKKS.cc)

target_link_libraries(testAlCuDatabases thermo4pfm)
target_link_libraries(testInterpolationFunctions thermo4pfm)
target_link_libraries(testCALPHADFunctions thermo4pfm)
target_link_libraries(testCALPHADSpeciesPhaseGibbsEnergyExpansion thermo4pfm)
target_link_libraries(testCALPHADGibbsEnergy thermo4pfm)
target_link_libraries(testCALPHADbinaryKKS thermo4pfm)
target_link_libraries(testCALPHADbinaryKKSThreePhase thermo4pfm)
target_link_libraries(testCALPHADbinaryKKS3Ph2Sl thermo4pfm)
target_link_libraries(testCALPHADbinaryEquilibrium thermo4pfm)
target_link_libraries(testCALPHADConcSolverBinaryThreePhase thermo4pfm)
target_link_libraries(testCALPHADConcSolverBinary3Ph2Sl thermo4pfm)
target_link_libraries(testCALPHADternaryKKS thermo4pfm)
target_link_libraries(testCALPHADConcSolverTernary thermo4pfm)
target_link_libraries(testCALPHADternaryTieLine thermo4pfm)
target_link_libraries(testKKSdiluteBinaryKKS thermo4pfm)
target_link_libraries(testKKSdiluteBinaryEquilibrium thermo4pfm)

target_link_libraries(testLoopCALPHADSpeciesPhaseGibbsEnergy thermo4pfm)
target_link_libraries(testLoopCALPHADbinaryKKS thermo4pfm)
target_link_libraries(testLoopCALPHADbinaryEquilibrium thermo4pfm)
target_link_libraries(testLoopCALPHADternaryKKS thermo4pfm)
target_link_libraries(testLoopCALPHADternaryTieLine thermo4pfm)
target_link_libraries(testCALPHADJacobianBinary thermo4pfm)
target_link_libraries(testCALPHADJacobianTernary thermo4pfm)
target_link_libraries(testLoopKKSdiluteBinaryKKS thermo4pfm)

if (${OPENMP_CXX_FOUND})
  target_link_libraries(testCALPHADSpeciesPhaseGibbsEnergyExpansion OpenMP::OpenMP_CXX)
  if(${WITH_OPENMP_OFFLOAD})
    target_link_libraries(testOpenMPoffload OpenMP::OpenMP_CXX)
  endif()
  target_link_libraries(testLoopCALPHADSpeciesPhaseGibbsEnergy OpenMP::OpenMP_CXX)
  target_link_libraries(testLoopCALPHADbinaryEquilibrium OpenMP::OpenMP_CXX)
  target_link_libraries(testLoopCALPHADbinaryKKS OpenMP::OpenMP_CXX)
  target_link_libraries(testLoopCALPHADternaryTieLine OpenMP::OpenMP_CXX)
  target_link_libraries(testLoopCALPHADternaryKKS OpenMP::OpenMP_CXX)
  target_link_libraries(testLoopKKSdiluteBinaryKKS OpenMP::OpenMP_CXX)
endif()

add_test(NAME testPropertyTree
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testPropertyTree
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/property_tree)

add_test(NAME testAlCuDatabases
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testAlCuDatabases
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testInterpolationFunctions
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testInterpolationFunctions)
add_test(NAME testCALPHADFunctions
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADFunctions)
add_test(NAME testCALPHADSpeciesPhaseGibbsEnergyExpansion
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADSpeciesPhaseGibbsEnergyExpansion
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADGibbsEnergy
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADGibbsEnergy
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADbinaryKKS
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADbinaryKKS
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADbinaryKKSThreePhase
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADbinaryKKSThreePhase
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADbinaryKKS3Ph2Sl
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADbinaryKKS3Ph2Sl
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADbinaryEquilibrium
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADbinaryEquilibrium
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADConcSolverBinaryThreePhase
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADConcSolverBinaryThreePhase
          WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADConcSolverBinary3Ph2Sl
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADConcSolverBinary3Ph2Sl
        WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testKKSdiluteBinaryKKS
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testKKSdiluteBinaryKKS
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testKKSdiluteBinaryEquilibrium
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testKKSdiluteBinaryEquilibrium
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADJacobianBinary
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADJacobianBinary
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADJacobianTernary
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADJacobianTernary
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADternaryKKS
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADternaryKKS
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADConcSolverTernary
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADConcSolverTernary
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testCALPHADternaryTieLine
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADternaryTieLine
         WORKING_DIRECTORY ${TEST_DIR})
if(${WITH_OPENMP_OFFLOAD})
  add_test(NAME testOpenMPoffload
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
           ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
           ${CMAKE_CURRENT_BINARY_DIR}/testOpenMPoffload)
endif()
add_test(NAME testLoopCALPHADSpeciesPhaseGibbsEnergy
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testLoopCALPHADSpeciesPhaseGibbsEnergy
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testLoopCALPHADbinaryEquilibrium
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testLoopCALPHADbinaryEquilibrium
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testLoopCALPHADbinaryKKS
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testLoopCALPHADbinaryKKS
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testLoopCALPHADternaryKKS
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testLoopCALPHADternaryKKS
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testLoopCALPHADternaryTieLine
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testLoopCALPHADternaryTieLine
         WORKING_DIRECTORY ${TEST_DIR})
add_test(NAME testLoopKKSdiluteBinaryKKS
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROCS_FLAG}
         ${MPIEXEC_NUMPROCS} ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/testLoopKKSdiluteBinaryKKS
         WORKING_DIRECTORY ${TEST_DIR})
