#!/usr/bin/env python3

# @file   test_pybind_3d.py
#
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Fri Jan 29 17:27:10 2021
#
# @brief  pytest for 3D interfaces
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

import sys
import pytest
import numpy as np
import cracklet as cra
from fixtures import model_3D

def test_creation_model_3D(model_3D):

    sim_driver = cra.SimulationDriver(model_3D);

    beta = 0.2 # By default
    time = 0 / 512 / 1263 * beta
    
    assert([512,32] == model_3D.getNbElements())
    assert(time == model_3D.getTime())
    assert(0 == model_3D.getCurrentTimeStep())

    dx = 1 / 512
    
    assert([dx,dx] == model_3D.getElementSize())
    assert(3 == model_3D.getDim())

if __name__ == "__main__":
    pytest.main(sys.argv)
