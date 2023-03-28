cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(GenCosmic)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS EG)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---VMC
# set by VMC standalone, in ROOT case we don't use it
set(VMC_LIBRARIES)
if(NOT ROOT_vmc_FOUND)
  find_package(VMC REQUIRED)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.)

set(HEADERS GeneratorCosmics.h)

ROOT_GENERATE_DICTIONARY(G__GeneratorCosmics ${HEADERS} LINKDEF GeneratorCosmicsLinkDef.h)

#---Create a shared library with geneated dictionary
add_library(GeneratorCosmics SHARED GeneratorCosmics.cxx G__GeneratorCosmics.cxx)
target_link_libraries(GeneratorCosmics ${ROOT_LIBRARIES} ${VMC_LIBRARIES})


set_target_properties(GeneratorCosmics
  PROPERTIES
  PUBLIC_HEADER "${HEADERS}" )


install(TARGETS GeneratorCosmics
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libGeneratorCosmics_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libGeneratorCosmics.rootmap
      DESTINATION lib)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(GeneratorCosmics PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)
