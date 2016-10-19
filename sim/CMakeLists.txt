#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  configuration file of cRacklet simulations folder
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

option(UIUC_ALU_HOMA "Activate Aluminium-Homalite simulations" ON) 
option(2D_HETEROG "Activate 2-D heterogenous interface simulations" ON)
option(3D_ASPERTITY "Activate 3-D interface with asperities simulations" ON)

if(UIUC_ALU_HOMA)
  add_subdirectory(001_UIUC_AluHoma)
endif()

if(2D_HETEROG)
  add_subdirectory(002_2d_heterog)
endif()

if(3D_ASPERTITY)
  add_subdirectory(003_3d_asperity)
endif()
