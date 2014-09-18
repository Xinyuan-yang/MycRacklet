#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  Macro for cRacklet tests execution
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#===============================================================================
set(cRacklet_DIFF_SCRIPT ${cRacklet_CMAKE_DIR}/diff.sh)
#===============================================================================
macro(register_cRacklet_test test_name description)

  add_executable(${test_name} ${test_name}.cc)
  target_link_libraries(${test_name} cRacklet)
  
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
    add_test(${test_name} ${cRacklet_DIFF_SCRIPT} ${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
  endif()

endmacro()