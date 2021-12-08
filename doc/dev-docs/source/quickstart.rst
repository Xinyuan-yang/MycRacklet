Quickstart
==========

Here is a quick introduction to get started with cRacklet.

Installation from source
------------------------

The following dependencies are required for cRacklet:

 - a **C++ compiler**
 - **CMake**
 - **FFTW3**
 - **GSL** library

If you have Debian/Ubuntu based system you can install them with::

  apt install cmake gcc libfftw3-dev libgsl-dev

If you have a MacOS::

  brew install cmake gcc fftw gsl

f you have a windows system you can use windows sub system Linux and check apt instruction.
  
Optional dependencies are:

 - **pybind11** (for python binding, automatically installed)
 - **python3** (for python binding)
 - a compiler that supports **OpenMP** (for multi-threaded parallel computing)
 - **pytest** (for tests)
 - **Doxygen** and **Sphinx** (for documentation)
 - a **Fortran** compiler to generate new kernels
  
First clone the git repository::

  git clone https://gitlab.com/cracklet/cracklet
  
You can then compile cRacklet using cMake::

  cd cracklet
  mkdir build
  cd build
  ccmake .. or cmake ..
  [ Set the options that you need ]
  make
  
Building the docs
-----------------

To build the documentation locally, activate the documentation option with **CMake** (CRACKLET_DOCUMENTATION). Then from the build directory::
  
  make dev_doc

Running the test
----------------

You need to activate the test options with **CMake** (CRACKLET_TESTS). You can then run the test with::

  make test
