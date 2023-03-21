#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  configuration file of cRacklet simulations folder
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

option(CRACKLET_EXAMPLES_BIMATERIAL "Activate bimaterial (Aluminium-Homalite) simulations" OFF)
option(CRACKLET_EXAMPLES_2D_HETEROG "Activate 2-D heterogenous interface simulations" ON)
option(CRACKLET_EXAMPLES_3D_ASPERTITY "Activate 3-D interface with asperities simulations" OFF)
option(CRACKLET_EXAMPLES_CRACK_FRONT_WAVES "Activate simulations of the crack front waves study" OFF) 
option(CRACKLET_EXAMPLES_FRICTION "Activate simulation of frictional contact problems" ON)
option(CRACKLET_EXAMPLES_VISCOELASTIC "Activate simulation of viscoelastic interfaces" OFF)
option(CRACKLET_EXAMPLES_ROUGH_SURFACES "Activate simulation of microcontacts failure" OFF)
option(CRACKLET_EXAMPLES_PROFILING "Executables to analyze code performance" ON)

if(CRACKLET_EXAMPLES_BIMATERIAL)
  add_subdirectory(001_UIUC_AluHoma)
endif()

if(CRACKLET_EXAMPLES_2D_HETEROG)
  add_subdirectory(002_2d_heterog)
endif()

if(CRACKLET_EXAMPLES_3D_ASPERTITY)
  add_subdirectory(003_3d_asperity)
endif()

if(CRACKLET_EXAMPLES_CRACK_FRONT_WAVES)
  add_subdirectory(004_crack_front_waves)
endif()

if(CRACKLET_EXAMPLES_FRICTION)
  add_subdirectory(005_friction)
endif()

if(CRACKLET_EXAMPLES_VISCOELASTIC)
  add_subdirectory(006_viscoelastic)
endif()

if(CRACKLET_EXAMPLES_ROUGH_SURFACES)
  add_subdirectory(007_rough_surfaces)
endif()

if(CRACKLET_EXAMPLES_PROFILING)
  add_subdirectory(profiling)
endif()

if(CRACKLET_PYTHON_INTERFACE)
  option(CRACKLET_EXAMPLES_PYTHON "Activate examples on how to use the python interface" OFF)

  if(CRACKLET_EXAMPLES_PYTHON)
    add_subdirectory(python)
  endif()
endif()
