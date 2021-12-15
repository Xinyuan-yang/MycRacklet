#!/usr/bin/env python3

# @file   fixtures.py
#
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Wed Dec 15 15:02:10 2021
#
# @brief  fixtures for pytest
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

import pytest
import numpy as np
import cracklet as cra

# Create a fixture model class for a 2D interface
@pytest.fixture
def model_2D():
    nb_time_steps = 100
    nb_x = 512
    dom_size = 1
    nu = 0.33
    E = 5.3e9
    cs = 1263
    tcut = 100
    model = cra.SpectralModel(nb_x,nb_time_steps,dom_size,nu,E,cs,tcut,"Fixture 2D Model")
    return model

# Create a fixture model class for a 3D interface
@pytest.fixture
def model_3D():
    nb_time_steps = 100
    nb_x = 512
    nb_z = 32
    dom_size_x = 1
    dom_size_z = 0.0625
    nu = 0.33
    E = 5.3e9
    cs = 1263
    tcut = 100
    model = cra.SpectralModel([nb_x,nb_z],nb_time_steps,[dom_size_x,dom_size_z],nu,E,cs,tcut,"Fixture 2D Model")
    return model
