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

Seeking help
------------

You can ask your questions on `c4science <https://c4science.ch>`_ using this
`form
(ADD LINK)
If you do not have an account, you can create one `here
<https://c4science.ch/auth/start/?next=%2F>`_.      

Contribution
------------

Code
....

To contribute code to cRacklet, you can use `Arcanist
<https://secure.phabricator.com/book/phabricator/article/arcanist/>`_ to send
code differentials. In a nutshell, the process to contribute is:

1. Create a branch for the modifications you wish to submit
2. Work on your branch (commits + run tests)
3. ``arc diff`` to send your code for review
4. Commit any requested changes
5. ``arc diff`` to send your modifications

For reviewers:

1. Checkout a code differential using ``arc patch D???``
2. Accept the code differential on `c4science <https://c4science.ch>`_.
3. ``arc land`` to merge the differential
4. Profit with ``arc anoid``

Bug reports
...........

You can also contribute to Tamaas by reporting any bugs you find `here
(ADD LINK)
if you have an account on `c4science <https://c4science.ch>`_.
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
