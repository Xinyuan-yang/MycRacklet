Examples
========

The directory ``examples/`` in cRacklet' root repository contains example scripts
dedicated to various aspects of cRacklet:

Bi-material interface
---------------------

Dynamic debonding at a planar interface btw Aluminium (mtl) and Homalite (poly). This example is a 1D interface loaded in mixed mode (Mode I and Mode II). 

Super-shear transition and constant driving speed
-------------------------------------------------

This folder contains two script, dealing with a 1D interface between two similar bodies. The interface follows a linear coupled cohesive law. The properties of the interface are spatially heterogeneous and consist of successive weak and strong stripes.

andrews_transition.cc
^^^^^^^^^^^^^^^^^^^^^

Study super-shear transition following Andrews(1976) formalism (accelerating mode II shear crack that transition to super-shear velocity via the nucleation of a daughter crack ahead of the primary crack)

cst_speed_fracture.cc
^^^^^^^^^^^^^^^^^^^^^

Dynamic fracture at a constant crack speed in presence of heterogeneities. The loading is updated dynamically such that the rupture propagates at a given velocity (sub-Rayleigh)

3D Interface in the presence of a strong asperity
-------------------------------------------------

3d interface in the presence of one asperity.

Rate and state friction laws
----------------------------

Rupture nucleation on an interface following a generic rate and state friction law. The interface is initially sliding at a steady-state velocity and nucleation occurs via a slight perturbation of the velocity field.
