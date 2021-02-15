include_directories( ${CMAKE_SOURCE_DIR}/src )

add_executable(plotEnergyVsComposition
               ${CMAKE_SOURCE_DIR}/drivers/plotEnergyVsComposition.cc )

target_link_libraries(plotEnergyVsComposition thermo4pfm)
