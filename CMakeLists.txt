cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(AEGIS)

add_subdirectory(MICROCERN)

add_subdirectory(TEPEMGEN)

add_subdirectory(GeneratorCosmics)

add_subdirectory(GeneratorParam)

if (DEFINED ENV{HIJING_ROOT})
   add_subdirectory(THijing)
endif (DEFINED ENV{HIJING_ROOT})

if (DEFINED ENV{DPMJET_ROOT})
   add_subdirectory(TDPMjet)
endif (DEFINED ENV{DPMJET_ROOT})

#Install(DIRECTORY Geometry gconfig DESTINATION share/Detectors/)
