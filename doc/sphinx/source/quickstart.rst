Quickstart
----------

Here is a quick introduction to get started with cRacklet.

Installation from source
^^^^^^^^^^^^^^^^^^^^^^^^

The following dependencies are required for cRacklet:

- a **C++ compiler**
- **cMake**
- **FFTW3**
- **GSL** library
- **pybind11** (automatically installed)
  
Optional dependencies are:

- **OpenMP** (for multi-threaded parallel computing)
- **pytest** (for tests)
- **Doxygen** and **Sphinx** (for documentation)

First clone the git repository::

  git clone htpps://c4science.ch/source/cRacklet.git

You can then compile cRacklet using cMake::

  mkdir build
  cd build
  ccmake ..
  make
  
Building the docs
^^^^^^^^^^^^^^^^^

TBA

Running the test
^^^^^^^^^^^^^^^^^

You need to activate the test options with cMake (CRACKLET_TESTS)
