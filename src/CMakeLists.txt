set(SOURCES ${CMAKE_SOURCE_DIR}/src/CALPHADFunctions.cc
            ${CMAKE_SOURCE_DIR}/src/xlogx.cc)

add_library(thermo4pfm ${SOURCES})

install(TARGETS thermo4pfm DESTINATION ${CMAKE_SOURCE_DIR}/lib)

