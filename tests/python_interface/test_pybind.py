#!/usr/bin/env python

# @file   test_pybind.py
#
# @author Thibault Roch <thibault.roch@epfl.ch>
#
# @date   Fri Jan 29 17:27:10 2021
#
# @brief  configuration file of cRacklet python interface test
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

# Create a fixture model class
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

def test_creation_model_2D(model_2D):
    
    sim_driver = cra.SimulationDriver(model_2D);

    beta = 0.2 # By default
    time = 1 / 512 / 1263 * beta
    
    assert(512 == model_2D.getNbElements()[0])
    assert(time == model_2D.getTime())
    assert(1 == model_2D.getCurrentTimeStep())

    dx = 1 / 512
    
    assert(dx == model_2D.getElementSize()[0])
    assert(3 == model_2D.getDim())

def test_uniform_loading(model_2D):

    load = 100
    
    sim_driver = cra.SimulationDriver(model_2D);
    interfacer = cra.InterfacerLinearCoupledCohesive(model_2D)
    cra.DataRegister.registerParameterReal("critical_normal_opening",0.02e-3)
    cra.DataRegister.registerParameterReal("critical_shear_opening",0.02e-3)
    cra.DataRegister.registerParameterReal("max_normal_strength",5e6)
    cra.DataRegister.registerParameterReal("max_shear_strength",5e6)
    
    interfacer.createUniformInterface()
    sim_driver.initConstantLoading(load,0,0)

    assert([0,load,0] == model_2D.getUniformLoading())

# 3D
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

def test_creation_model_3D(model_3D):

    sim_driver = cra.SimulationDriver(model_3D);

    beta = 0.2 # By default
    time = 1 / 512 / 1263 * beta
    
    assert([512,32] == model_3D.getNbElements())
    assert(time == model_3D.getTime())
    assert(1 == model_3D.getCurrentTimeStep())

    dx = 1 / 512
    
    assert([dx,dx] == model_3D.getElementSize())
    assert(3 == model_3D.getDim())

if __name__ == "__main__":
    pytest.main(sys.argv)
