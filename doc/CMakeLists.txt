#===============================================================================
# @file   CMakeLists.txt
#
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date creation: Fri May 22 2020
# @date last modification: Fri Apr 2 2021
#
# @brief  configuration file for documentation
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
#
# @section DESCRIPTION
#
#===============================================================================

find_package(Sphinx REQUIRED)

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/dev-docs/source/conf.py"
  "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/dev-docs/source/Doxyfile"
  "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/dev-docs/source/input_def.in"
  "${CMAKE_CURRENT_BINARY_DIR}/input_def"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/dev-docs/source/xml_output_def.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xml_output_def"
  @ONLY)

add_custom_target(dev_doc
  ${SPHINX_EXECUTABLE}
  -j 4
  -q -b html
  -c "${CMAKE_CURRENT_BINARY_DIR}"
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/dev-docs/source"
  "${SPHINX_HTML_DIR}"
  COMMENT "Building HTML documentation with Sphinx")
