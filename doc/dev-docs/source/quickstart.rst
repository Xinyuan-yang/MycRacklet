Quickstart
==========

Here is a quick introduction to get started with cRacklet.

Installation from source
------------------------

The following dependencies are required for cRacklet:

- a **C++ compiler**
- **cMake**
- **FFTW3**
- **GSL** library
  
Optional dependencies are:

- **pybind11** (for python binding, automatically installed)
- **python3** (for python binding)
- **OpenMP** (for multi-threaded parallel computing)
- **pytest** (for tests)
- **Doxygen** and **Sphinx** (for documentation)
- a **Fortran** compiler to generate new kernels
  
First clone the git repository::

  git clone https://gitlab.com/cracklet/cracklet
  
You can then compile cRacklet using cMake::

  mkdir build
  cd build
  ccmake ..
  make
  
Building the docs
-----------------

To build the documentation locally, activate the documentation option with **cMake** (CRACKLET_DOCUMENTATION). Then from the build directory::
  
  make dev_doc

Running the test
----------------

You need to activate the test options with cMake (CRACKLET_TESTS). You can then run the test with::

  make test
