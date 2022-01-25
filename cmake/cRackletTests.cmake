#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  Macro for cRacklet tests execution
#
# @section LICENSE
#
# Copyright (©) 2012 - 2013 Fabian Barras
#               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#===============================================================================
set(cRacklet_DIFF_SCRIPT ${cRacklet_CMAKE_DIR}/diff.sh)
#===============================================================================
set(_test_flags
  PYTHON
  )

macro(register_cRacklet_test test_name description)
  cmake_parse_arguments(_register_test
    "${_test_flags}"
    "${_test_one_variables}"
    "${_test_multi_variables}"
    ${ARGN}
    )
  
  set(_arguments -n "${test_name}")
  
  if(_register_test_PYTHON)  
    # Source the python environement
    list(APPEND _arguments -E "${PROJECT_BINARY_DIR}/cRacklet_environement.sh")
    configure_file("${test_name}.py" "${test_name}.py" COPYONLY)
    list(APPEND _arguments -e "${test_name}.py")
  else()
    # Create executable for cpp tests
    add_executable(${test_name} ${test_name}.cc)
    target_link_libraries(${test_name} cRacklet)
  endif()
  
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified")
    list(APPEND _arguments -r "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified")
  endif()
  
  add_test(NAME ${test_name} COMMAND ${cRacklet_DIFF_SCRIPT} ${_arguments})

endmacro()
