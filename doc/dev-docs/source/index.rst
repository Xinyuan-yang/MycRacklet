.. cRacklet documentation master file, created by
   sphinx-quickstart on Wed Feb  3 18:37:50 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cRacklet - A boundary integral library for interfacial rupture simulations
##########################################################################

.. image:: https://gitlab.com/cracklet/cracklet/badges/master/pipeline.svg
   :target: https://gitlab.com/cracklet/cracklet/-/commits/master

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gl/cracklet%2Ftutorials/master?filepath=supershear%2Fsupershear.ipynb
      
*cRacklet* cRacklet is a C++ boundary integral library based on a spectral formulation of the dynamic wave equations in two semi-infinite linearly elastic solids `Geubelle and Rice (1995) <http://www.sciencedirect.com/science/article/pii/002250969500043I>`_  `Breitenfeld and Geubelle (1998) <https://link.springer.com/article/10.1023/A:1007535703095>`_ . Its implementation is specially tailored for the modeling of dynamic crack/rupture propagation along a planar interface bonding the two half-spaces. The main benefit of this spectral method is a numerical discretization limited to the mid-plane, thereby providing a very fine description of the dynamic rupture processes, unattainable with finite-element or finite-difference schemes. For more details about the method and its applications, refer to the publications section.

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   ./quickstart
   ./simulation
   ./interface_law
   ./python_interface
   ./examples
   ./tutorials
   ./api_reference
   ./publications
   ./authors

Seeking help - Bug reports
--------------------------

You can ask your questions or report any bugs you find here https://gitlab.com/cracklet/cracklet/-/issues

Code contribution
-----------------

Any contribution to cRacklet, bugfixes or new features, is welcome with a merge request on GitLab.
   
Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
