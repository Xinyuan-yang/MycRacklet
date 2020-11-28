---
title: 'cRacklet: a spectral boundary integral method library for interfacial rupture simulation'
tags:
  - boundary integral
  - dynamic rupture
  - elastodynamics
  - friction
  - c++
  - python
authors:
  - name: Thibault Roch
    orcid: 0000-0002-2495-8841
    affiliation: 1
  - name: Fabian Barras
    orcid: 0000-0003-1109-0200
    affiliation: "1, 2"
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
date: 28 November 2020
bibliography: paper.bib

---

# Summary

The study of dynamically propagating rupture along faults is of prime importance in fields ranging from engineering to geosciences. Numerical simulations of these phenomena are computationally costly and challenging. A fine discretization in time and space is required to accurately represent singularities and discontinuities near rupture edges but simulations must also involve larger length scales such as the propagation length. Also, the behavior of such interfaces can be highly non-linear thus increasing the problem complexity. Conventional numerical approaches for fracture problem, for instance, the use of cohesive elements in finite-element method [@ortiz_finite-deformation_1999], requires to discretize the whole body containing the fault and is consequently computationally expensive. The use of boundary integral methods, reducing the dimensionality of the problem, enable to focus the computational efforts on the fracture plane and allows for a detailed description of the interfacial field's evolution.

# Statement of need

`cRacklet` is a C++ library with a Python interface [@pybind11] developed as a collaboration between the Computational Solid Mechanics Laboratory at EPFL and the Department of Aerospace Engineering of the University of Illinois at Urbana-Champaign that implements a spectral formulation of the elastodynamics boundary integral relations between the displacements and the corresponding traction stress acting at a planar interface between two homogeneous elastic solids [@geubelle_spectral_1995], [@breitenfeld_numerical_1998]. The stresses acting on the interfaces are partly computed with a convolution of the history of the interfacial displacement in the Fourier domain, computed efficiently with FFTW3/OpenMP. The prescription of an interfacial behavior allows solving for the equilibrium at a given time step. Time integration is achieved using an explicit time-stepping scheme. cRacklet is aimed at researchers interested in interfacial dynamics, ranging from nucleation problem to dynamic propagation of rupture fronts.

# Features

cRacklet allows for planar rupture interface simulations loaded in any combination of normal traction, in-plane and out-of-plane shear solicitations. cRacklet handles the simulation of interfaces bonded between dissimilar elastic solids. Any stress or material heterogeneity along the fracture plane can be resolved using cRacklet. Several interfacial behaviors are included in the library, such as:

- Slip-weakening laws [@ida_cohesive_1972] [@palmer_growth_1973]. This behavior can be coupled with a classical Coulomb friction law or a regularized one [@prakash_frictional_1998] to handle friction emerging from the contact of the two surfaces.

- Rate and state dependant friction laws, including the original formulation by [@dieterich_modeling_1979] and [@ruina_slip_1983] and the N-shape formulation proposed in [@aldam_critical_2017]

Existing softwares dealing with interfacial rupture are based on a simplified formulation, known as the quasi-dynamic approximation [@qdyn]. We are not aware of any public software package including the implementation of the full elastodynamics equations using the so-called boundary integral spectral method.

# Performance

Comparison with Akantu?

# Example

Add a nice figure with interesting behavior... (Space-time diagram of something...)

# Publications

The following publications have been made possible with cRacklet:

- @barras_study_2014

- @barras_interplay_2017

- @brener_unstable_2018

- @barras_emergence_2019

- @barras_emergence_2020

- @fekak_crack_2020

# Acknowledgments

We acknowledge the financial support of the Swiss National Science Foundation (grants #162569 and ...)

# References
