"""**
 * @file   modeI.py
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Mon Oct 26 14:34:07 2020
 *
 * @brief  Python script for a pure mode I solicitation of an interface between similar materials 
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * 
 * cRacklet is the result of a collaboration between the Computational Solid Mechanics 
 * Laboratory (LSMS) of Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland 
 * and the Department of Aerospace Engineering of the University of Illinois at 
 * Urbana-Champaign, United States of America.
 * 
 * cRacklet is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License as published by the Free Software Foundation, 
 * either version 3 of the License, or (at your option) any later version.
 * 
 * cRacklet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program.  
 * If not, see <http://www.gnu.org/licenses/>.
 */
"""

import numpy as np

from py_cRacklet import SpectralModel
from py_cRacklet import InterfacerLinearCoupledCohesive
from py_cRacklet import CohesiveLaw
from py_cRacklet import DataDumper
from py_cRacklet import DataRegister
from py_cRacklet import DataFields
from py_cRacklet import SimulationDriver

def main():

    # Geometry and material parameters
    
    nb_time_steps = 15000
    nb_elements = 2048
    dom_size = 0.3
    dx = dom_size / nb_elements
    initial_crack_size = 50 * dx
    propagation_domain = 0.9 * dom_size

    nu = 0.35
    E = 5.3e9
    cs = 1263

    # Cut of the loaded material kernels
    
    tcut = 100
    
    # Loading case

    load = 3e6
    psi = 0
    phi = 0

    # Cohesive parameters

    crit_n_open = 0.02e-3
    crit_s_open = 0.02e-3

    max_s_str = 5e6
    max_n_str = 5e6

    # Griffith crack length
    
    G_length = crit_n_open*max_n_str/(load*load*np.pi)*E/(1-nu*nu)

    print("Griffith crack length is {:.5f}".format(G_length))
    print("Domain size is {:.2f}".format(dom_size))
    
    # Creating the model
    
    model = SpectralModel(nb_elements,nb_time_steps,dom_size,nu,E,cs,tcut,"Mode I debonding at Homalite-Homalite interface")

    # Simulation driver

    sim_driver = SimulationDriver(model);
    
    interfacer = InterfacerLinearCoupledCohesive(model)

    DataRegister.registerParameterReal("critical_normal_opening",crit_n_open)
    DataRegister.registerParameterReal("critical_shear_opening",crit_s_open)
    DataRegister.registerParameterReal("max_normal_strength",max_n_str)
    DataRegister.registerParameterReal("max_shear_strength",max_s_str)
    
    interfacer.createUniformInterface()
    interfacer.createThroughCrack(0,initial_crack_size)
    interfacer.createThroughWall(propagation_domain,dom_size)
    
    cohesive_law = model.getInterfaceLaw()

    sim_driver.initConstantLoading(load,psi,phi)
    
    # Dumper configuration
    
    dumper = DataDumper(model)
    
    st_diag_id = 'ST_Diagram_id.cra'
    st_diag_nor_trac = 'ST_Diagram_normal_tractions.cra'
    st_diag_nor_velo = 'ST_Diagram_normal_velocitiy_jumps.cra'
    top_u = 'top_displ_snapshot.cra'
    
    dumper.initVectorDumper(st_diag_nor_trac,DataFields._interface_tractions,1)
    dumper.initDumper(st_diag_nor_velo,DataFields._normal_velocity_jumps)
    dumper.initDumper(st_diag_id,DataFields._id_crack)
    dumper.initDumper(top_u,DataFields._top_displacements)

    print_freq = 0.02 * nb_time_steps

    # Launch the crack up to dynamic propagation

    sim_driver.launchCrack(0.,G_length,0.1)

    end_propagation = 0.9*propagation_domain

    x_tip = 0
    t = 0
    
    while (x_tip < end_propagation):

        sim_driver.solveStep()
        x_tip = DataRegister.getCrackTipPosition(0,nb_elements) * dx
        
        if ((t%print_freq)==0):
            print("Process at {:.3f} %".format(t/nb_time_steps*100))
            print("Crack position at {:.3f} %".format(x_tip/dom_size*100))

            dumper.dumpAll()

        t += 1


def plotEvolution():

    print("To be implemented")
        
if __name__ == '__main__':

    main()
    plotEvolution()
