find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/src )

set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

add_executable(testPropertyTree
      ${TEST_DIR}/property_tree/testPropertyTree.cc)

add_executable(testCALPHADFunctions
               ${TEST_DIR}/testCALPHADFunctions.cc
               ${CMAKE_SOURCE_DIR}/src/xlogx.cc
               ${CMAKE_SOURCE_DIR}/src/CALPHADFunctions.cc )
add_executable(testCALPHADGibbsEnergy
               ${TEST_DIR}/CALPHADGibbsEnergy/testCALPHADGibbsEnergy.cc
               ${CMAKE_SOURCE_DIR}/src/CALPHADSpeciesPhaseGibbsEnergy.cc
               ${CMAKE_SOURCE_DIR}/src/CALPHADSpeciesPhaseGibbsEnergyExpansion.cc
               )

add_test(NAME testPropertyTree
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testPropertyTree
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/property_tree)

add_test(NAME testCALPHADFunctions
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADFunctions)
add_test(NAME testCALPHADGibbsEnergy
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADGibbsEnergy
         WORKING_DIRECTORY ${TEST_DIR}/CALPHADGibbsEnergy)
