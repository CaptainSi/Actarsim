#--------------------CMake file for ActarSim in Geant 4.9.6.p02--------------
# Author: D. Perez-Loureiro 28/08/2013
# Modified and updated for newer cmake, root and geant4 versions 2015-2017
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ACTARSim)
set(CMAKE_MACOSX_RPATH 1)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Cmake Path for installation
#
#CMAKE_INSTALL_PREFIX
#----------------------------------------------------------------------------
# Load some basic macros which are needed later on
include(FindROOT.cmake)

#---------------------------------------------------------------------------
# Create the directory for the ROOT files
#
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/root_files")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/root_files/sim_files")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/root_files/dig_files")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/root_files/out_files")

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
# Setup ROOT include directories

include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
#Setup ROOT linking directories
#
link_directories( ${ROOT_LIBRARY_DIR})

#----------------------------------------------------------------------------
#Setup of compiler flags
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
   EXEC_PROGRAM("sw_vers -productVersion | cut -d . -f 1-2" OUTPUT_VARIABLE MAC_OS_VERSION)
   MESSAGE("--- Found a Mac OS X System ${MAC_OS_VERSION}")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

      STRING(COMPARE EQUAL "10.5" "${MAC_OS_VERSION}" MAC_OS_10_5)
      IF(MAC_OS_10_5)
        SET(CMAKE_CXX_FLAGS "-m64")
        SET(CMAKE_Fortran_FLAGS "-m64")
      ENDIF(MAC_OS_10_5)

      SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
      SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3 -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O2")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3 -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")



      endif (CMAKE_COMPILER_IS_GNUCXX)

endif (CMAKE_SYSTEM_NAME MATCHES Darwin)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Set everything needed for the root dictonary and create the
# dictionary
set(ActarSim_HEADERS ${PROJECT_SOURCE_DIR}/include/ActarSimData.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimTrack.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimSimpleTrack.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimSilHit.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimSilRingHit.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimSciHit.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimSciRingHit.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimPlaHit.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimPrimaryInfo.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimCrossSectionVariable.hh
		     ${PROJECT_SOURCE_DIR}/include/ActarSimBeamInfo.hh
		     )
set(ActarSim_LINKDEF ${PROJECT_SOURCE_DIR}/src/ActarSimLinkDef.hh )
set(ActarSim_DICTIONARY ActarSimDict.cc)
ROOT_GENERATE_DICTIONARY("${ActarSim_HEADERS}" "${ActarSim_LINKDEF}" "${ActarSim_DICTIONARY}" "${ROOT_INCLUDE_DIR}")


#----------------------------------------------------------------------------
# Set the library version in the main CMakeLists.txt
SET(ActarSim_MAJOR_VERSION 0)
SET(ActarSim_MINOR_VERSION 0)
SET(ActarSim_PATCH_VERSION 0)
SET(ActarSim_VERSION "${ActarSim_MAJOR_VERSION}.${ActarSim_MINOR_VERSION}.${ActarSim_PATCH_VERSION}")
SET(ActarSim_LIBRARY_PROPERTIES ${ActarSim_LIBRARY_PROPERTIES}
    VERSION "${ActarSim_VERSION}"
    SOVERSION "${ActarSim_MAJOR_VERSION}"
    SUFFIX ".sl"
)





add_library(actar SHARED ActarSim.cc ${ActarSim_DICTIONARY} ${sources} ${headers})
target_link_libraries(actar ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
set_target_properties(actar PROPERTIES ${ActarSim_LIBRARY_PROPERTIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(actarsim ActarSim.cc ${ActarSim_DICTIONARY} ${sources} ${headers})
target_link_libraries(actarsim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_dependencies(actarsim actar)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build ActarSim. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(ACTARSIM_SCRIPTS
gui.mac
vis_actarTPCDEMO.mac
vis_actarTPC.mac
vis_maiko.mac
#vis_kine.mac
vis_alpha.mac
vis_simpleBox.mac
vis_simpleTube.mac
kine_batch.mac
12C_batch.mac
digitizationMacro.C
digit.h
readerPads.C
runDigitization.C
runReaderPads.C
MTrack.h
Reducer.C
Analysis_Reduce.C
Analysis_Angle.C
Analysis_Profile.C
 )

foreach(_script ${ACTARSIM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add DOXYGEN DOCUMENTATION
# The documentation is created in the binary (build) directory subdir /html
# To read the documentation open the index.html file using a browser
#
#
#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#add_custom_target(doc ALL
#${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#COMMENT "Generating API documentation with Doxygen" VERBATIM
#)
#endif(DOXYGEN_FOUND)

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(ActarSim DEPENDS ActarSim)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS actarsim DESTINATION bin)
#install(TARGETS actarsim DESTINATION lib)

install(FILES ${headers}
#        DESTINATION "./include"
        DESTINATION "${PROJECT_BINARY_DIR}/include"
        )
