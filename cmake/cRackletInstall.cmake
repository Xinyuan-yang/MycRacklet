#===============================================================================
# @file   cRackletInstall.txt
#
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Wed June 24 15:18:20 2020
#
# @brief Create the files that allows users to link with cRacklet in an other
# cmake project
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

configure_file(cmake/cRackletBuildTreeSettings.cmake.in
  "${PROJECT_BINARY_DIR}/cRackletBuildTreeSettings.cmake" @ONLY)

# Create the cRackletConfig.cmake and cRackletConfigVersion files
get_filename_component(CONF_REL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)
configure_file(cmake/cRackletConfig.cmake.in "${PROJECT_BINARY_DIR}/cRackletConfig.cmake" @ONLY)

# Copy environement file for python interface

configure_file(cmake/cRacklet_environement.sh.in
  ${PROJECT_BINARY_DIR}/cRacklet_environement.sh  @ONLY)

configure_package_config_file(cmake/cRackletConfig.cmake.in
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME}
  )

install(FILES
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindFFTW.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindSphinx.cmake
  ${PROJECT_SOURCE_DIR}/cmake/Modules/FindGSL.cmake
  ${PROJECT_BINARY_DIR}/cRackletConfig.cmake
  ${PROJECT_SOURCE_DIR}/cmake/cRackletSimulationMacros.cmake
  ${PROJECT_SOURCE_DIR}/cmake/cRackletTests.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME}
  COMPONENT dev)
