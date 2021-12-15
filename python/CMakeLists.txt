#===============================================================================
# @file   CMakeLists.txt
#
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Mon Feb 03 15:10:20 2020
#
# @brief  main configuration file of cracklet python library
#
# @section LICENSE
#
# cracklet - A spectral boundary integral method for interface fracture simulation
# Copyright (©) 2012 - 2013 Fabian Barras
#               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# 
# cracklet is the result of a collaboration between the Computational Solid Mechanics 
# Laboratory (LSMS) of Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland 
# and the Department of Aerospace Engineering of the University of Illinois at 
# Urbana-Champaign, United States of America.
# 
# cracklet is free software: you can redistribute it and/or modify it under the terms 
# of the GNU General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version.
# 
# cracklet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.  
# If not, see <http://www.gnu.org/licenses/>.

#===============================================================================

set(PYCRACKLET_SRCS
  py_cracklet.cc
  py_interface_law.cc
  py_data_register.cc
  py_spectral_model.cc
  py_simulation_driver.cc
  py_interfacer.cc
  py_contact_law.cc
  py_coulomb_law.cc
  py_regularized_coulomb_law.cc
  py_cohesive_law.cc
  py_cohesive_law_viscoelastic.cc
  py_rate_and_state_law.cc
  py_data_dumper.cc
  py_crack_profile.cc
  )

pybind11_add_module(py11_cracklet ${PYCRACKLET_SRCS})
target_include_directories(py11_cracklet
  SYSTEM BEFORE
  PRIVATE ${PYBIND11_INCLUDE_DIR}
  PRIVATE ${pybind11_INCLUDE_DIR}
  PRIVATE ${PYTHON_INCLUDE_DIRS})

target_link_libraries(py11_cracklet PUBLIC cRacklet)
set_target_properties(py11_cracklet PROPERTIES
  DEBUG_POSTIFX ""
  LIBRARY_OUTPUT_DIRECTORY cracklet)

file(COPY cracklet DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(_python_install_dir
  ${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

install(TARGETS py11_cracklet
  LIBRARY DESTINATION ${_python_install_dir})

install(DIRECTORY cRacklet
  DESTINATION ${_python_install_dir}
  FILES_MATCHING PATTERN "*.py")
