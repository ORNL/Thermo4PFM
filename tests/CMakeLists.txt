find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/src )

add_executable(testPropertyTree
      ${CMAKE_SOURCE_DIR}/tests/property_tree/testPropertyTree.cc)

add_executable(testCALPHADFunctions
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADFunctions.cc
               ${CMAKE_SOURCE_DIR}/src/xlogx.cc
               ${CMAKE_SOURCE_DIR}/src/CALPHADFunctions.cc )

add_test(NAME testPropertyTree
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testPropertyTree
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/property_tree)

add_test(NAME testCALPHADFunctions
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADFunctions)
