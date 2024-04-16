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
add_cracklet_simulation(mode_II_rate_and_state mode_II_rate_and_state.cc
  NU_TOP 33
  NU_BOT 33)

add_cracklet_simulation(mode_III_rate_and_state mode_III_rate_and_state.cc)

add_cracklet_simulation(travelling_pulse travelling_pulse.cc)
copy_input_file(input_pulse_PMMA.dat)

add_cracklet_simulation(mode_II_slip_weakening mode_II_slip_weakening.cc)

add_cracklet_simulation(mode_II_slip_weakening_new mode_II_slip_weakening_new.cc)

add_cracklet_simulation(mode_II_slip_weakening_new_coulomb mode_II_slip_weakening_new_coulomb.cc)

add_cracklet_simulation(mode_III_slip_weakening mode_III_slip_weakening.cc)

add_cracklet_simulation(mode_III_slip_weakening_old mode_III_slip_weakening_old.cc)

add_cracklet_simulation(mode_III_slip_weakening_new mode_III_slip_weakening_new.cc)

add_cracklet_simulation(mode_III_slip_weakening_new_with_load mode_III_slip_weakening_new_with_load.cc)

add_cracklet_simulation(mode_III_slip_weakening_new_coulomb mode_III_slip_weakening_new_coulomb.cc)

add_cracklet_simulation(mode_III_slip_weakening_load mode_III_slip_weakening_load.cc)

add_cracklet_simulation(mode_III_slip_weakening_load_equivalent mode_III_slip_weakening_load_equivalent.cc)

add_cracklet_simulation(mode_III_slip_weakening_v_const mode_III_slip_weakening_v_const.cc)

add_cracklet_simulation(mode_III_slip_weakening_restart mode_III_slip_weakening_restart.cc)

add_cracklet_simulation(mode_III_slip_weakening_non_uniform_load mode_III_slip_weakening_non_uniform_load.cc)

add_cracklet_simulation(mode_III_slip_weakening_non_uniform_load_increasing mode_III_slip_weakening_non_uniform_load_increasing.cc)

add_cracklet_simulation(mode_III_slip_weakening_ref_dual mode_III_slip_weakening_ref_dual.cc)
add_cracklet_simulation(mode_III_slip_weakening_ref_dual_equivalent mode_III_slip_weakening_ref_dual_equivalent.cc)

add_cracklet_simulation(mode_III_slip_weakening_new_with_load_equivalent mode_III_slip_weakening_new_with_load_equivalent.cc)

add_cracklet_simulation(mode_III_slip_weakening_new_input mode_III_slip_weakening_new_input.cc)

add_cracklet_simulation(mode_III_slip_weakening_new_with_load_increasing mode_III_slip_weakening_new_with_load_increasing.cc)

add_cracklet_simulation(mode_III_slip_weakening_new_with_load_increasing_2 mode_III_slip_weakening_new_with_load_increasing_2.cc)

add_cracklet_simulation(mode_III_slip_weakening_new_with_load_equivalent_increasing mode_III_slip_weakening_new_with_load_equivalent_increasing.cc)

add_cracklet_simulation(mode_III_slip_weakening_nucleation mode_III_slip_weakening_nucleation.cc)

add_cracklet_simulation(mode_III_slip_weakening_equivalent_GL mode_III_slip_weakening_equivalent_GL.cc)

add_cracklet_simulation(mode_III_slip_weakening_GL mode_III_slip_weakening_GL.cc)

add_cracklet_simulation(mode_III_slip_weakening_exp mode_III_slip_weakening_exp.cc)

add_cracklet_simulation(mode_III_slip_weakening_linear mode_III_slip_weakening_linear.cc)

add_cracklet_simulation(mode_III_slip_weakening_pointload mode_III_slip_weakening_pointload.cc)

add_cracklet_simulation(mode_III_slip_weakening_dual_cohesive mode_III_slip_weakening_dual_cohesive.cc)

add_cracklet_simulation(mode_III_slip_weakening_fluidpression mode_III_slip_weakening_fluidpression.cc)

copy_input_file(input_pulse_PMMA.dat)
