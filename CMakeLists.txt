cmake_minimum_required(VERSION 3.1)
project (CALPHAD4PFM CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the location of additional CMAKE modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(OpenMP)

# clang-format (optional)
set(WITH_CLANG_FORMAT FALSE CACHE BOOL "Indent code with clang-format")
if(${WITH_CLANG_FORMAT})
  find_package(CLANG_FORMAT)
  if(${CLANG_FORMAT_FOUND})
    message(STATUS "Indent with clang-format")
    file(GLOB_RECURSE FORMAT_SOURCES src/*.cc src/*.h tests/*.cc tests/*.h)
    add_custom_target(format
      COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file ${FORMAT_SOURCES}
      DEPENDS ${FORMAT_SOURCES})
  endif(${CLANG_FORMAT_FOUND})
endif()

set(WITH_TLOGT4L FALSE CACHE BOOL "Include T*log(T) in L coeffs")
if(${WITH_TLOGT4L})
  add_definitions(-DMAX_POL_T_INDEX=3)
  add_definitions(-DLMIX_WTLOGT)
else()
  add_definitions(-DMAX_POL_T_INDEX=2)
endif(${WITH_TLOGT4L})

enable_testing()

add_subdirectory(src)

add_subdirectory(tests)

add_subdirectory(drivers)
