#===============================================================================
# @file   CMakeLists.txt
#
# @author Fabian Barras <fabian.barras@epfl.ch>
#
# @date   Mon Sept 1 9:12:20 2014
#
# @brief  configuration file for UIUC simulations folder
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

add_cracklet_simulation(mode_II_rate_and_state mode_II_rate_and_state.cc
  NU_TOP 33
  NU_BOT 33)

add_cracklet_simulation(mode_III_rate_and_state mode_III_rate_and_state.cc)

add_cracklet_simulation(travelling_pulse travelling_pulse.cc)

add_cracklet_simulation(mode_II_slip_weakening mode_II_slip_weakening.cc)

add_cracklet_simulation(mode_III_slip_weakening mode_III_slip_weakening.cc)
