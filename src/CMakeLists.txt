#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  configuration file of cRacklet sources
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

list(APPEND CRACKLET_LIBRARY_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}")
set(CRACKLET_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR} ${CRACKLET_LIBRARY_INCLUDE_DIRS}
  CACHE INTERNAL "Internal include directories to link with cRacklet as a subproject")

set(cRacklet_SRC
  data_register.cc
  spectral_model.cc
  crack_profile.cc
  interfacer.cc
  convolution_manager.cc
  spectral_convolution_manager.cc
  cohesive_law.cc
  cohesive_law_viscoelastic.cc
  rate_and_state_law.cc
  coulomb_law.cc
  regularized_coulomb_law.cc
  data_dumper.cc
  simulation_driver.cc
  )

set(CRACKLET_LIBRARY_PUBLIC_HDRS
  cohesive_law.hh
  cohesive_law_inline_impl.hh
  cohesive_law_viscoelastic.hh
  cohesive_law_viscoelastic_formulations.hh
  cohesive_law_viscoelastic_inline_impl.hh
  contact_law.hh
  convolution_manager.hh
  convolution_manager_inline_impl.hh
  coulomb_law.hh
  cRacklet_common.hh
  crack_profile.hh
  data_dumper.hh
  data_dumper_inline_impl.hh
  data_register.hh
  data_register_inline_impl.hh
  interface_law.hh
  interfacer.hh
  interfacer_inline_impl.hh
  rate_and_state_law.hh
  rate_and_state_formulations.hh
  regularized_coulomb_law.hh
  ring_buffer.hh
  ring_buffer_inline_impl.hh
  simulation_driver.hh
  simulation_driver_inline_impl.hh
  spectral_convolution_manager.hh
  spectral_model.hh
  spectral_model_inline_impl.hh
  )

set_source_files_properties(${RELEASE_FILE} PROPERTIES GENERATED 1)

add_library(cRacklet ${cRacklet_SRC} ${RELEASE_FILE})
add_dependencies(cRacklet RELEASE_INFO)

target_link_libraries(cRacklet PUBLIC ${FFTW_LIBRARIES} ${GSL_LIBRARIES})

target_include_directories(cRacklet PUBLIC ${FFTW_INCLUDE_PATH})

target_include_directories(cRacklet PUBLIC ${GSL_INCLUDE_PATH})

if(CRACKLET_USE_LIBSURFER)
  target_link_libraries(cRacklet PUBLIC ${LIBSURFER_LIBRARIES})
endif()

if(CRACKLET_MULTI_THREADED)
  add_definitions(${OpenMP_CXX_FLAGS})
  target_link_libraries(cRacklet PUBLIC ${FFTW_OPENMP_LIBRARY} ${OpenMP_CXX_FLAGS})
endif()

#target_include_directories(cRacklet
#  PUBLIC ${CRACKLET_INCLUDE_DIRS}
#  )

target_include_directories(cRacklet
  PRIVATE   $<BUILD_INTERFACE:${CRACKLET_INCLUDE_DIRS}>
  INTERFACE $<INSTALL_INTERFACE:include/cracklet>
  )

# small trick for build includes in public
set_property(TARGET cRacklet APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CRACKLET_INCLUDE_DIRS}>)


set_target_properties(cRacklet
  PROPERTIES
  INTERFACE_CXX_STANDARD 14
  )

if(CRACKLET_LIBRARY_PUBLIC_HDRS)
  set_property(TARGET cRacklet PROPERTY PUBLIC_HEADER ${CRACKLET_LIBRARY_PUBLIC_HDRS})
endif()

if(CRACKLET_LIBRARY_PRIVATE_HDRS)
  set_property(TARGET cRacklet PROPERTY PRIVATE_HEADER ${CRACKLET_LIBRARY_PRIVATE_HDRS})
endif()


install(TARGETS cRacklet
  EXPORT ${CRACKLET_TARGETS_EXPORT}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cracklet/
  )

if("${CRACKLET_TARGETS_EXPORT}" STREQUAL "cRackletTargets")

  install(EXPORT cRackletTargets
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME}
    COMPONENT dev)

  #Export for build tree
  export(EXPORT cRackletTargets
  FILE "${CMAKE_BINARY_DIR}/cRackletTargets.cmake")
  export(PACKAGE cRacklet)
endif()
