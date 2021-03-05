Simulation
----------

Find cRacklet and add a simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can create your own **cMake** project and link it to **cRacklet** by adding::

  find_package(cRacklet REQUIRED)

to your **CMakeLists.txt**. After doing so you can create an executable with::

  add_cracklet_simulation(executable source
                          NU_TOP nu_t
                          NU_BOT nu_b)

with NU_TOP and NU_BOT optional arguments that, if provided, will copy the corresponding kernels to your simulation folder. nu_t and nu_b are the value of the Poisson modulus associated with respectively the top and the bottom material.

Only specific kernels are already generated with the distributed version of cRacklet. You can generate new kernels using the Fortran routine **invert_serial.f** found in the **pre-computed-kernels/** folder in the root of **cRacklet** repository. You will need to change the inputs used by **invert_serial.f** in **invert_serial.in**.

Simulation structure
^^^^^^^^^^^^^^^^^^^^

In order to run a **cRacklet**, you will need first to create a SpectralModel object. This class is the core of cRacklet librairy and contains the methods processing the different steps required to solve the elastodynamic response of the two semi-infinite half space as desribed in [2].


Python interface
^^^^^^^^^^^^^^^^

In order to use the python interface build for cRacklet, you need pybind11 and python3.

During configuration, activate CRACKLET_PYTHON_INTERFACE. The make command will build a python library in the **build/python/** folder. Please add this path to your python path.

You can activate the python example with the option **CRACKLET_EXAMPLES_PYTHON**.

