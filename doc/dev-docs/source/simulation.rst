Simulation
==========

Find cRacklet and add a simulation
----------------------------------

You can create your own **cMake** project and link it to **cRacklet** by adding::

  find_package(cRacklet REQUIRED)

to your **CMakeLists.txt**. After doing so you can create an executable with::

  add_cracklet_simulation(executable source
                          NU_TOP nu_t
                          NU_BOT nu_b)

with NU_TOP and NU_BOT optional arguments that, if provided, will copy the corresponding kernels to your simulation folder. nu_t and nu_b are the value of the Poisson modulus associated with respectively the top and the bottom material.

Only specific kernels are already generated with the distributed version of cRacklet. You can generate new kernels using the Fortran routine **invert_serial.f** found in the **pre-computed-kernels/** folder in the root of **cRacklet** repository. You will need to change the inputs used by **invert_serial.f** in **invert_serial.in**.

Simulation structure
--------------------

Create the model
^^^^^^^^^^^^^^^^

In order to run a **cRacklet** simulation, you will need first to create a :cpp:class:`SpectralModel <SpectralModel>` object. This class is the core of cRacklet librairy and contains the methods processing the different steps required to solve the elastodynamic response of the two semi-infinite half space as desribed in [2].

The model for a 2D simulation with similar material on top and bottom is the following::

  #include "spectral_model.hh"

  int main(){
  
  SpectralModel model(nb_elements,nb_time_steps,dom_size,nu,E,cs,tcut,simulation_summary,output_dir);

  // your code ....

  }
  
with ``nb_elements`` the number of discretization points, ``nb_time_steps`` the number of time steps, ``dom_size`` the real length of the interface, ``nu``, ``E``, ``cs`` the characteristics of the elastic bulk (respectively the Poisson ratio, the Young Modulus and the shear wave speed velocity), ``t_cut`` the cut in the kernels, ``simulation_summary`` a string containing a description of your simmulation and ``output_dir`` the folder in which the output will be written.

It is then required to init the model and setting the loading case::

  model.initModel();
  model.setLoadingCase(load,psi,phi);

with ``psi`` and ``phi`` the orientation of the load and ``load`` its absolute value.

Create the interface
^^^^^^^^^^^^^^^^^^^^

The second step is to create and define the initial interface properties, with the help of the :cpp:class:`Interfacer <Interfacer>` object which is templated by one of the possible :cpp:enum:`FractureLawType <FractureLawType>`. The corresponding headers have to be included (here an example for the linear coupled cohesive law)::
  
  #include "interfacer.hh"
  #include "cohesive_law.hh"

  int main(){

  // Model creation...

  Interfacer<_linear_coupled_cohesive> interfacer(model);
  
  }
  
Each interface law has build-in methods method allowing for the creation of specific interface layout, see the API documentation for more details. The simplest case is to create a spatially homogeneous interface. It requires some parameters to be registered in the model first: this can either be done manually with :cpp:func:`registerParameter <DataRegister::registerParameter>` of the :cpp:class:`DataRegister <DataRegister>` class or using the built-in reader :cpp:func:`readInputFile <DataRegister::readInputFile>` that parse an input file and register every parameters. As an example, to create an uniform interface with a cohesive law behavior, one can write::

  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("max_normal_strength",max_n_str);
  DataRegister::registerParameter("max_shear_strength",max_s_str);
  interfacer.createUniformInterface();

with ``critical_normal_opening``, ``critical_shear_opening``, ``max_normal_strength`` and ``max_shear_strength`` being parameters describing the interface response.

Additionaly, one can add a pre-inserted crack with the following call, setting the resistance to 0 between ``start_crack`` and ``end_crack``  ::

  interfacer.createThroughCrack(start_crack,end_crack);

Then, the loads can be properly applied and the interface fields (stresses, velocities, displacements) can be initiated::

  model.updateLoads();
  model.initInterfaceFields();

Setting up the dumper
^^^^^^^^^^^^^^^^^^^^^

A :cpp:class:`DataDumper <DataDumper>` object needs to be created. One can then select the fields to dump with the various method associted to the DataDumper class. Here we give an example where the dumper is configure to output the displacement field of the top elastic body in the ``output`` file. The fields that can be dumped are found in :cpp:enum:`DataFields <DataFields>` ::

  DataDumper dumper(model);
  dumper.initDumper("output", _top_displacements);

Solving steps
^^^^^^^^^^^^^

Solving the equilibrium of stresses at the interface requires to execute several functions from the object::

  model.updateDisplacements();
  model.fftOnDisplacements();
  model.computeStress();
  model.computeInterfaceFields();
  model.increaseTimeStep();

The output can be generated by calling::

  dumper.dumpAll()

Simulation Driver
^^^^^^^^^^^^^^^^^
  
Alternatively, one can ressort to :cpp:class:`SimulationDriver <SimulationDriver>` object to handle all the necessary operations required to solve a step. The constructor of the :cpp:class:`SimulationDriver <SimulationDriver>` object takes a :cpp:class:`SpectralModel <SpectralModel>` as an argument. The construction of the simulation driver include the model initialisation. The loading can be handle directly from the :cpp:class:`SimulationDriver <SimulationDriver>`. The whole solve step routine can then be reduce to a single call, as shown in the example below for a homogenous constant loading ::

  SimulationDriver sim_driver(model);  
  sim_driver.initConstantLoading(load,psi,phi);  
  sim_driver.solveStep();
