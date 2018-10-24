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

set(cRacklet_SRC
  data_register.cc
  spectral_model.cc
  crack_profile.cc
  interfacer.cc
  convolution_manager.cc
  spectral_convolution_manager.cc
  cohesive_law.cc
  rate_and_state_law.cc
  coulomb_law.cc
  regularized_coulomb_law.cc
  data_dumper.cc
  simulation_driver.cc
  )

set_source_files_properties(${RELEASE_FILE} PROPERTIES GENERATED 1)

add_library(cRacklet ${cRacklet_SRC} ${RELEASE_FILE})
add_dependencies(cRacklet RELEASE_INFO)

target_link_libraries(cRacklet ${FFTW_LIBRARIES} ${GSL_LIBRARIES})

if(CRACKLET_USE_LIBSURFER)
  target_link_libraries(cRacklet ${LIBSURFER_LIBRARIES})
endif()

if(CRACKLET_MULTI_THREADED)
  add_definitions(${OpenMP_CXX_FLAGS})
  target_link_libraries(cRacklet ${FFTW_OPENMP_LIBRARY} ${OpenMP_CXX_FLAGS})
endif()
