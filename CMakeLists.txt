#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  main configuration file of cRacklet Library
#
# @section LICENSE
#
# cRacklet - A spectral boundary integral method for interface fracture simulation
# Copyright (©) 2012 - 2013 Fabian Barras
#               2014 - ongoing EPFL (Ecole Polytechnique Fédérale de Lausanne)
#               Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
# 
# cRacklet is the result of a collaboration between the Computational Solid Mechanics 
# Laboratory (LSMS) of Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland 
# and the Department of Aerospace Engineering of the University of Illinois at 
# Urbana-Champaign, United States of America.
# 
# cRacklet is free software: you can redistribute it and/or modify it under the terms 
# of the GNU General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version.
# 
# cRacklet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.  
# If not, see <http://www.gnu.org/licenses/>.

#===============================================================================
cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_C_COMPILER "/usr/bin/gcc-12")
#set(CMAKE_CXX_COMPILER "/usr/bin/g++-12")
project(cRacklet)
enable_language(CXX Fortran)

cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  cmake_policy(SET CMP0074 NEW)
endif()

set(cRacklet_CMAKE_DIR "${cRacklet_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${cRacklet_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${cRacklet_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "-std=c++0x" CACHE STRING "" FORCE)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
mark_as_advanced(BUILD_SHARED_LIBS)

if(NOT CRACKLET_TARGETS_EXPORT)
  set(CRACKLET_TARGETS_EXPORT cRackletTargets)
endif()

include(GNUInstallDirs)

include(cRackletExtraCompilationProfiles)

#===============================================================================
# Release file target
#===============================================================================
set(GET_RELEASE_INFO ${PROJECT_SOURCE_DIR}/get_release_info.sh)
set(RELEASE_FILE ${PROJECT_BINARY_DIR}/release.cc)
set_property(SOURCE ${RELEASE_FILE} PROPERTY GENERATED 1)

find_program(uuencode_exec "uuencode")
mark_as_advanced(uuencode_exec)

#===============================================================================
# External Library
#===============================================================================
option(CRACKLET_MULTI_THREADED "Use multi-threaded version of cRacklet using OpenMP" OFF)

find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

include_directories(${GSL_INCLUDE_PATH})
include_directories(${FFTW_INCLUDE_PATH})

if(CRACKLET_MULTI_THREADED)
  add_definitions(-DCRACKLET_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

#===============================================================================
# Kernels pre-computation, now done in cRackletSimulationMacro
#===============================================================================
SET(KERNELS_DIR ${PROJECT_SOURCE_DIR}/pre-computed-kernels)

add_subdirectory(${KERNELS_DIR})

include(cRackletSimulationMacros)
include(cRackletTools)

#===============================================================================
# Python interface
#===============================================================================
option(CRACKLET_PYTHON_INTERFACE "Activate python binding" OFF)

if(CRACKLET_PYTHON_INTERFACE)
  find_package(PythonInterp 3)
  add_external_package(pybind11 VERSION 2.2.4 CONFIG)
  add_subdirectory(python)
endif()

#===============================================================================
# Library
#===============================================================================

add_subdirectory(src)

#===============================================================================
# Documentation
#===============================================================================

option(CRACKLET_DOCUMENTATION "Build documentation" OFF)

if(CRACKLET_DOCUMENTATION)
  add_subdirectory(doc)
endif()

#===============================================================================
# Release
#===============================================================================

add_custom_command(
  OUTPUT ${RELEASE_FILE}
  COMMAND ${GET_RELEASE_INFO} ${PROJECT_BINARY_DIR} > ${RELEASE_FILE}
  DEPENDS ${GET_RELEASE_INFO}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMENT "Generating the Release files"
   )

add_custom_target(RELEASE_INFO DEPENDS ${RELEASE_FILE})

#===============================================================================
# Examples and tests
#===============================================================================

option(CRACKLET_TESTS "Activate tests" OFF)
option(CRACKLET_EXAMPLES "Activate examples" ON)

if(CRACKLET_TESTS)
  enable_testing()
  include(CTest)
  include(cRackletTests)
  mark_as_advanced(BUILD_TESTING)
  add_subdirectory(tests)
endif()

if(CRACKLET_EXAMPLES)
  add_subdirectory(examples)
endif()

#===============================================================================
# Install
#===============================================================================

include(cRackletInstall)
