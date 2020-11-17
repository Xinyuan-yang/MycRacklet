title: 'cRacklet: a spectral boundary integral method library for interfacial rupture simulation'
tags:
  - interfaces
  - boundary integral
  - dynamic rupture
  - and something...
authors:
  - name: Thibault Roch
    orcid: 0000-0002-2495-8841
    affiliation: 1
  - name: Fabian Barras
    orcid: 0000-0003-1109-0200
    affiliation: "1,2"
  - name: Philippe H Geubelle
    orcid: 0000-0002-4670-5474
    affiliation: 3
  - name: Jean-François Molinari
    orcid: 0000-0002-1728-1844
    affiliation: 1
affiliations:
 - name: Civil Engineering Institute, École Polytechnique Fédérale de Lausanne, Switzerland
   index: 1
 - name: The Njord Centre Department of Physics, Department of Geosciences, University of Oslo, Norway
   index: 2
 - name: Department of Aerospace Engineering of the University of Illinois at Urbana-Champaign, United States of America
   index: 3
date: 17 November 2020
bibliography: paper.bib

# Summary

Boundary integral methods.



The numerical method allows for a detailed description of the evolution of the interacial fields especially in the failure zone captured with the aid of a failure model.

# Statement of need

`cRacklet` is a C++ library with a Python interface developed as a collaboration between the Computational Solid Mechanics Laboratory at EPFL and the the Department of Aerospace Engineering of the niversity of Illinois at Urbana-Champaign that implements a spectral formulation of the elastodynamics boundary integral relations between the displacements and the corresponding traction stress acting at a planar interface between two solids. [@geubelle_spectral_1995] , [@breitenfeld_numerical_1998]. The interfacial stress in this formulation is partly given by a convolution of the history of the displacement fields and these convolutions are computed efficiently using pre-computed kernel and parallel computation using FFTW3/OPENMPI. 

cRacklet is aimed at researchers studying rupture dynamics ...

# Features

cRacklet allows for planar rupture interface simulations loaded in any directions. cRacklet handle the simulation of interfaces bonded between dissimilar elastic solids. Any stress or material heterogeneity along the fracture plane can be resolved using cRacklet. Several interfacial behavior are included in the library, such as:
    - slip-weakening [ref] with the possibility to coupled with standard Coulomb
    - several formulation of the rate and state dependant friction laws [Ruina & Dieterich + regularized?]. 

We are not aware of any public software package including the implementation of the so-called spectral formulation of the elastodynamics boundary integral equations.

# Publications

The following publications have been made possible with cRacklet:

# Acknowledgements

# References