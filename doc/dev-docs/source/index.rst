.. cRacklet documentation master file, created by
   sphinx-quickstart on Wed Feb  3 18:37:50 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cRacklet - A boundary integral library for interfacial rupture simulations
==========================================================================

*cRacklet* cRacklet is a C++ boundary integral library based on a spectral formulation of the dynamic wave equations in two semi-infinite linearly elastic solids. Its implementation is specially tailored for the modeling of dynamic crack/rupture propagation along a planar interface bonding the two half-spaces. The main benefit of this spectral method is a numerical discretization limited to the mid-plane, thereby providing a very fine description of the dynamic rupture processes, unattainable with finite-element or finite-difference schemes. For more details about the method and its applications, refer to "related publications" section.

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   ./quickstart
   ./simulation
   ./examples
   ./api_reference

Seeking help - Bug reports
--------------------------

You can ask your questions or report any bugs you find here https://gitlab.com/tiburoch/cracklet/-/issues

Code contribution
-----------------

Any contribution to cRacklet, bugfixes or new features, is welcome with a gitlab pull request.
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
