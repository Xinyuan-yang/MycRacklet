#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
# @author Damien Spielmann <damien.spielmann@epfl.ch>
#
# @date   Mon Mar 2 9:12:20 2015
#
# @brief  configuration file of cRacklet tests
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

option(BIMA_FRACT "Bimaterial Alu-Homalite fracture test" ON)

if(BIMA_FRACT)
  add_subdirectory(bima_fract)
endif()

option(BIMA_FRACT_3D "Bimaterial Alu-Homalite fracture test in 3D" ON)

if(BIMA_FRACT_3D)
  add_subdirectory(bima_fract_3d)
endif()

option(ASPERITY_3D "Testing 3D fracture event with a circular asperity" ON)

if(ASPERITY_3D)
  add_subdirectory(asperity_3d)
endif()

option(CRACKLET_ENGINE "Testing core features of cRacklet library" ON)

if(CRACKLET_ENGINE)
  add_subdirectory(cRacklet_engine)
endif()

option(FRICTION "Testing frictional models of cRacklet library" ON)

if(FRICTION)
  add_subdirectory(friction)
endif()

option(CONSTANT_CRACK_SPEED "Testing constant crack speed algorithm" ON)

if(CONSTANT_CRACK_SPEED)
  add_subdirectory(cst_speed)
endif()
