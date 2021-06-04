Python interface
================

Compile the library
-------------------

In order to use the python interface build for cRacklet, you need pybind11 and python3.

During configuration, activate CRACKLET_PYTHON_INTERFACE. Running::

  make

in **build** folder will compile the python library called py_cRacklet in **build/python/** folder. Please add this path to your python path. Otherwise you can do ::

  make install
  
You can activate the python examples with the option **CRACKLET_EXAMPLES_PYTHON**.

Use the python interface
------------------------

Once the library is compiled and its path added to your python path, you can import cRacklet with::

  import cracklet

One can run a complete simulation using the python interface. The required steps to run a simulation are identical to the **c++** counterpart. First a :py:class:`SpectralModel <py_cRacklet.SpectralModel>` needs to be created::

  model = py_cRacklet.SpectralModel(nb_elements,nb_time_steps,dom_size,nu,E,cs,tcut,simulation_summary,output_dir)

A :py:class:`SimulationDriver <py_cRacklet.SimulationDriver>` can then be associated to the model::

  sim_driver = cracklet.SimulationDriver(model)
  
The parameters can be registered using :py:func:`registerParameterReal <py_cRacklet.DataRegister::registerParameterReal>`, and the :py:class:`Interfacer <py_cRacklet.InterfacerLinearCoupledCohesive>` can then be used to create the interface (here, an example with an homogeneous interface following a linear coupled cohesive law, with a broken area)::

  cracklet.DataRegister.registerParameterReal("critical_normal_opening",crit_n_open)
  cracklet.DataRegister.registerParameterReal("critical_shear_opening",crit_s_open)
  cracklet.DataRegister.registerParameterReal("max_normal_strength",max_n_str)
  cracklet.DataRegister.registerParameterReal("max_shear_strength",max_s_str)
  interfacer = py_cRacklet.InterfacerLinearCoupledCohesive(model)    
  interfacer.createUniformInterface()
  interfacer.createThroughCrack(crack_start,crack_end)

The loading can then be initiated with::

   sim_driver.initConstantLoading(load,psi,phi)
  
The :py:class:`DataDumper <py_cRacklet.DataDumper>` is used to dump the interfacial fields that are available in :py:enum:`DataFields <py_cRacklet.DataFields>` ::
   
  dumper = DataDumper(model)
  dumper.initDumper(output_file,py_cRacklet.DataFields._top_displacements)

One can then solve a step with ::
  
  sim_driver.solveStep()

and dump the results with::

  dumper.dumpAll()
