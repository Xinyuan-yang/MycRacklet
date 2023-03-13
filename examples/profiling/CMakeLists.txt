#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date   Wed Aug 19 9:12:20 2015
#
# @brief  configuration file for profiling folder
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

include_directories (${cRacklet_SOURCE_DIR}/src)

add_executable(size_checking size_checking.cc)
target_link_libraries (size_checking cRacklet)

add_executable(printKernels printKernels.cc)
target_link_libraries (printKernels cRacklet)

copy_kernel_files(33 TRUE)
copy_kernel_files(35 TRUE)
copy_kernel_files(35 FALSE)
