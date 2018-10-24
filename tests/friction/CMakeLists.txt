#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date   Wed May 3 17:12:20 2017
#
# @brief  configuration file of cRacklet frictional models
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

register_cRacklet_test(test_rate_and_state_PMMA "Rate and state instability along PMMA")
register_cRacklet_test(test_rate_and_state_PMMA_regularized "Regularized Rate and state formulation")

copy_kernel_files(33 TRUE)
copy_kernel_files(33 FALSE)
