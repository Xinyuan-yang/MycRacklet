#!/usr/bin/env python
import sys
import pytest
import numpy as np
import py_cRacklet as cra

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
