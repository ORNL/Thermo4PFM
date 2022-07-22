cmake_minimum_required(VERSION 3.9)
project (CALPHAD4PFM CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the location of additional CMAKE modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

SET(CMAKE_INSTALL_PREFIX "$(CMAKE_SOURCE_DIR}/install" CACHE PATH
    "Directory to install Thermo4PFM")

find_package(Boost)
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
endif()

#OpenMP
find_package(OpenMP REQUIRED)
set(WITH_OPENMP_OFFLOAD FALSE CACHE BOOL "Build with OpenMP offload")
if(${WITH_OPENMP_OFFLOAD})
  if(OpenMP_CXX_VERSION LESS "4.5")
    message(SEND_ERROR "OpenMP 4.5 or later required")
  endif()
  message(STATUS "Build with OpenMP offload")
  add_definitions(-DHAVE_OPENMP_OFFLOAD)
endif()

# clang-format (optional)
set(WITH_CLANG_FORMAT FALSE CACHE BOOL "Indent code with clang-format")
if(${WITH_CLANG_FORMAT})
  find_package(CLANG_FORMAT)
  if(${CLANG_FORMAT_FOUND})
    message(STATUS "Indent with clang-format")
    file(GLOB_RECURSE FORMAT_SOURCES src/*.cc src/*.h tests/*.cc tests/*.h drivers/*.cc)
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

set(SHARED_LIB FALSE CACHE BOOL "Whether or not build static lib")

enable_testing()

add_subdirectory(src)

add_subdirectory(tests)

add_subdirectory(drivers)
