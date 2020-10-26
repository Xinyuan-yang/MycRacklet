"""**
 * @file   bimaterial.py
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Wed Jul 1 11:40:07 2020
 *
 * @brief  Python script for the study of dynamic debonding at a planar interface btw Aluminium (mtl) and Homalite (poly) 
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
from py_cRacklet import RegularizedCoulombLaw
from py_cRacklet import InterfacerLinearCoupledCohesive
from py_cRacklet import CohesiveLaw
from py_cRacklet import DataDumper
from py_cRacklet import DataFields

def main():

    # Geometry and material parameters
    
    nb_time_steps = 7000
    nb_elements = 4096
    dom_size = 1
    crack_size = 0.05
    nu_mtl = 0.33
    nu_poly = 0.35
    E_mtl = 71e9
    E_poly = 5.3e9
    cs_mtl = 3100
    cs_poly = 1263

    # Cut of the loaded material kernels
    
    tcut_mtl = 100
    tcut_poly = 100
    
    # Loading case

    load = 3e6
    psi = 75
    phi = 0

    # Cohesive parameters

    crit_n_open = 0.02e-3
    crit_s_open = 0.02e-3

    max_s_str = 5e6
    max_n_str = 5e6

    # Friction parameters

    regularized_time_scale = 0.1
    coef_frict = 0.25

    contact_law = RegularizedCoulombLaw(coef_frict,regularized_time_scale,nb_elements)

    # Creating the model
    
    model = SpectralModel(nb_elements,nb_time_steps,dom_size,nu_mtl,nu_poly,E_mtl,E_poly,cs_mtl,cs_poly,tcut_mtl,tcut_poly,"Mixed-mode debonding at Aluminium Homalite interface")

    # Stable time step coefficient
    
    beta = 0.4

    model.initModel(beta)
    model.setLoadingCase(load,psi,phi)

    interfacer = InterfacerLinearCoupledCohesive(model)
    interfacer.createThroughCenteredCrack(crack_size,crit_n_open,crit_s_open,max_n_str,max_s_str)

    cohesive_law = model.getInterfaceLaw()
    cohesive_law.preventSurfaceOverlapping(contact_law)

    model.updateLoads()
    model.initInterfaceFields()

    # Dumper configuration
    
    dumper = DataDumper(model)
    
    st_diag_id = 'ST_Diagram_id.cra'
    st_diag_nor_trac = 'ST_Diagram_normal_tractions.cra'
    st_diag_shear_velo = 'ST_Diagram_shear_velocitiy_jumps.cra'
    top_u = 'top_displ_snapshot.cra'
    bot_u = 'bot_displ_snapshot.cra'
    tractions = 'trac_snapshots.cra'
    point_his = 'Points_history.cra'
    points_int = (np.arange(0,1,0.1)*nb_elements).astype(int)

    energy = 'Energy.cra'
    
    dumper.initVectorDumper(st_diag_nor_trac,DataFields._interface_tractions,1)
    dumper.initDumper(st_diag_shear_velo,DataFields._shear_velocity_jumps)
    dumper.initDumper(st_diag_id,DataFields._id_crack)
    dumper.initDumper(top_u,DataFields._top_displacements)
    dumper.initDumper(bot_u,DataFields._bottom_displacements)  
    dumper.initDumper(tractions,DataFields._interface_tractions)  
    dumper.initPointsDumper(point_his,points_int)  
    dumper.initIntegratorsDumper(energy)  

    print_freq = 0.1 * nb_time_steps

    for t in np.arange(0,nb_time_steps):
        model.updateDisplacements()
        model.fftOnDisplacements()
        model.computeStress()
        model.computeInterfaceFields()
        model.increaseTimeStep()

        dumper.dump(st_diag_id)
        dumper.dump(st_diag_nor_trac)
        dumper.dump(st_diag_shear_velo)
        dumper.dump(point_his)
        dumper.dump(energy)        

        if ((t%print_freq)==0):
            print("Process at " + str(t/nb_time_steps*100) + " %")

            dumper.dump(top_u)
            dumper.dump(bot_u)
            dumper.dump(tractions)
            
if __name__ == '__main__':

    main()
