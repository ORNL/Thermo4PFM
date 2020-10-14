find_package(PythonInterp 3.2 REQUIRED)

add_executable(testPropertyTree
      ${CMAKE_SOURCE_DIR}/tests/property_tree/testPropertyTree.cc)

add_test(NAME testPropertyTree
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testPropertyTree
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/property_tree)
