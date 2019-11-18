#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  main configuration file of cRacklet Library
#
# @section LICENSE
#
# cRacklet - A spectral boundary integral method for interface fracture simulation
# Copyright (©) 2012 - 2013 Fabian Barras
#               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
cmake_minimum_required(VERSION 2.8)
project(cRacklet)
enable_language(CXX Fortran)

set(cRacklet_CMAKE_DIR "${cRacklet_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${cRacklet_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${cRacklet_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "-std=c++0x" CACHE STRING "" FORCE)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
mark_as_advanced(BUILD_SHARED_LIBS)

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
option(CRACKLET_USE_LIBSURFER "Use LibSurfer tools to generate fractal surfaces" OFF)
option(CRACKLET_MULTI_THREADED "Use multi-threaded version of cRacklet using OpenMP" OFF)

set(LIBSURFER_FOLDER /home/fabarras/Softwares/LibSurfer)
set(LIBSURFER_BINARIES /home/fabarras/Softwares/LibSurfer/build)

find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

include_directories(${GSL_INCLUDE_PATH})
include_directories(${FFTW_INCLUDE_PATH})

if(CRACKLET_USE_LIBSURFER)
  add_definitions(-DCRACKLET_USE_LIBSURFER)
  find_package(LIBSURFER REQUIRED)
  include_directories(${LIBSURFER_INCLUDE_DIRS})
endif()

if(CRACKLET_MULTI_THREADED)
  add_definitions(-DCRACKLET_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

#===============================================================================
# Kernels pre-computation
#===============================================================================
set(KERNELS_DIR ${PROJECT_SOURCE_DIR}/pre-computed-kernels)

add_subdirectory(${KERNELS_DIR})

set(KERNELS_EXT
  h11
  h22
  k12
 )

function(copy_kernel_files poisson_ratio is_bottom)
  foreach(_ext ${KERNELS_EXT})
    file(COPY ${KERNELS_DIR}/nu_.${poisson_ratio}_${_ext}.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    if(is_bottom)
      file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/nu_.${poisson_ratio}_${_ext}.dat ${CMAKE_CURRENT_BINARY_DIR}/nu_.${poisson_ratio}_${_ext}b.dat)
    endif()
  endforeach()
endfunction()

#===============================================================================
# Python interface
#===============================================================================
option(CRACKLET_PYTHON_INTERFACE "Activate python binding" OFF)

if(CRACKLET_PYTHON_INTERFACE)
  add_subdirectory(python)
endif()

#===============================================================================
# Library
#===============================================================================

add_subdirectory(src)

add_custom_command(
  OUTPUT ${RELEASE_FILE}
  COMMAND ${GET_RELEASE_INFO} ${PROJECT_BINARY_DIR} > ${RELEASE_FILE}
  DEPENDS ${GET_RELEASE_INFO}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMENT "Generating the Release files"
   )

add_custom_target(RELEASE_INFO DEPENDS ${RELEASE_FILE})

option(CRACKLET_TESTS "Activate tests" OFF)
option(CRACKLET_SIMULATIONS "Activate simulations" ON)

if(CRACKLET_TESTS)
  enable_testing()
  include(CTest)
  include(cRackletTests)
  mark_as_advanced(BUILD_TESTING)
  add_subdirectory(tests)
endif()

if(CRACKLET_SIMULATIONS)
  add_subdirectory(sim)
endif()
