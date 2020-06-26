#===============================================================================
# @file   CMakeLists.txt
#
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date creation: Fri May 22 2020
# @date last modification: Fri May 22 2020
#
# @brief  configuration file for doxygen
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

set(DOXYGEN_INPUT_DOX ${CMAKE_CURRENT_BINARY_DIR}/cracklet.dox)
set(DOXYGEN_HTML_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(DOXYGEN_OUTPUT ${DOXYGEN_HTML_DIR}/index.html)

# ---------------------------------------------------------------------------- #
# Doxygen                                                                      #
# ---------------------------------------------------------------------------- #
find_package(Doxygen REQUIRED)

set(DOXYGEN_WARNINGS NO)
set(DOXYGEN_QUIET YES)
if(CMAKE_VERBOSE_MAKEFILE)
  set(DOXYGEN_WARNINGS YES)
  set(DOXYGEN_QUIET NO)
endif(CMAKE_VERBOSE_MAKEFILE)

if (CMAKE_VERSION VERSION_GREATER 3.9.5)

  #set(DOXYGEN_WARNINGS YES)
  #set(DOXYGEN_QUIET NO)
  set(DOXYGEN_STRIP_FROM_PATH ${PROJECT_SOURCE_DIR})
  set(DOXYGEN_STRIP_FROM_INC_PATH ${PROJECT_SOURCE_DIR})
  set(DOXYGEN_ALIASES
    "rst=\\verbatim embed:rst"
    "endrst=\\endverbatim"
    )
  set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
  set(DOXYGEN_WARN_IF_DOC_ERROR YES)
  set(DOXYGEN_WARN_AS_ERROR NO)
  set(DOXYGEN_RECURSIVE YES)
  #set(DOXYGEN_EXAMPLE_PATH "${PROJECT_SOURCE_DIR}/examples")
  set(DOXYGEN_EXAMPLE_RECURSIVE YES)
  set(DOXYGEN_SOURCE_BROWSER NO)
  set(DOXYGEN_CLANG_ASSISTED_PARSING NO)
  #set(DOXYGEN_CLANG_OPTIONS          )
  set(DOXYGEN_CLANG_DATABASE_PATH ${CMAKE_BINARY_DIR})
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_GENERATE_HTML YES)
  set(DOXYGEN_GENERATE_HTMLHELP NO)
  set(DOXYGEN_GENERATE_LATEX NO)
  set(DOXYGEN_GENERATE_XML NO)
  set(DOXYGEN_HTML_OUTPUT html)
  set(DOXYGEN_ENABLE_PREPROCESSING YES)
  set(DOXYGEN_MACRO_EXPANSION YES)
  set(DOXYGEN_COLLABORATION_GRAPH NO)
  set(DOXYGEN_UML_LOOK YES)
  set(DOXYGEN_TEMPLATE_RELATIONS YES)
  set(DOXYGEN_CALL_GRAPH YES)
  set(DOXYGEN_CALLER_GRAPH YES)
  set(DOXYGEN_DOT_GRAPH_MAX_NODES 500)
  set(DOXYGEN_SHOW_FILES NO)
  set(DOXYGEN_LOOKUP_CACHE_SIZE 9)
  
  set(_SRCS
    ${PROJECT_SOURCE_DIR}
    )
    
  doxygen_add_docs(doxygen-doc
    ${_SRCS}
    COMMENT "Building HTML documentation with Doxygen in ${DOXYGEN_HTML_DIR}"
    )

else()

  message(${PROJECT_BINARY_DIR})
  
  make_directory(${DOXYGEN_HTML_DIR})

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cracklet.dox.in
    ${DOXYGEN_INPUT_DOX}
    )

  add_custom_command(
    OUTPUT ${DOXYGEN_OUTPUT}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT_DOX}
    DEPENDS ${DOXYGEN_INPUT_DOX}
    COMMENT "Building HTML documentation with Doxygen in ${DOXYGEN_HTML_DIR}"
    )

  add_custom_target(doxygen-doc ALL
    DEPENDS ${DOXYGEN_OUTPUT}
    )

  add_custom_target(doxygen-doc-forced
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
    COMMENT "Building HTML documentation with Doxygen (forced) in ${DOXYGEN_XML_DIR}"
    )
endif()
