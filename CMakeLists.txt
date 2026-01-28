#----------------------------------------------------------------------------
# Setup the project

cmake_minimum_required(VERSION 3.16...3.21)
project(PositronSource)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/processes_OC_FS/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/processes_OC_FS/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/processes_OC_FS/include/*.hh)
file(GLOB_RECURSE macros RELATIVE ${PROJECT_SOURCE_DIR} macros/*.mac)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries

add_executable(positronSource positronSource.cc ${sources} ${headers})
target_link_libraries(positronSource ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all macro/scripts to the build directory.
# This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

foreach(_file ${macros})
	configure_file(
    	${PROJECT_SOURCE_DIR}/${_file}
    	${PROJECT_BINARY_DIR}/${_file}
    	COPYONLY
	)
endforeach()

set(my_SCRIPTS
  execute.sh
  G4_scan_conventional.sh
  G4_scan.sh
  G4_scan_crystalline.sh
  G4_scan_GT.sh
)
foreach(_script ${my_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone

#add_custom_target(PositronSource DEPENDS positronsource)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX

install(TARGETS positronSource DESTINATION bin)

#----------------------------------------------------------------------------
# Create directiories in the binary path

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/input) #directory where to put input file to be read

#----------------------------------------------------------------------------
# Copy directiories in the binary path

file(COPY data DESTINATION ${PROJECT_BINARY_DIR}) #directory where to put all potentials 
#file(COPY input DESTINATION ${PROJECT_BINARY_DIR}) #directory where to put input file to be read

#----------------------------------------------------------------------------

