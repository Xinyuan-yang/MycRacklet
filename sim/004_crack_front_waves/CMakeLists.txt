#===============================================================================

# @file   CMakeLists.txt
# @author Fabian Barras <fabian.barras@mn.uio.no>
# @date   Thu Jan 16 14:52:14 2020
#
# @brief  configuration file for crack front waves simulations
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

#===============================================================================

add_cracklet_simulation(cst_speed_fracture_2d_3d cst_speed_fracture_2d_3d.cc
  NU_TOP 35
  NU_BOT 35)

#Copy the given kernel files in the build directory
#For the bottom material with \nu=0.35
#For the top material with \nu=0.35
#The kernel file need to be available in the folder {cRacklet_SOURCE_DIR}/pre-computed-kernels
#Otherwise new kernel files can be generated with {cRacklet_SOURCE_DIR}/pre-computed-kernels/invert_serial.f

copy_input_file(input_file.cra) #Copy the given input file in the build directory

