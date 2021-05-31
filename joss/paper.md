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
date: 18 May 2021
bibliography: paper.bib

---

# Summary

The study of dynamically propagating rupture along faults is of prime importance in fields ranging from engineering to geosciences. Numerical simulations of these phenomena are computationally costly and challenging. A fine spatial discretization is needed to represent accurately the singular fields associated with the rupture edges. Besides, the problems of interest usually involve a larger length scale along which rupture will propagate (i.e. a tectonic fault). The physical phenomena in play also occur at different timescale, from the slow process or rupture nucleation to the fast travel of crack front close the elastic wave speeds. Large and finely discretized spatio-temporal domains are required, which are computationally costly. In addition, the behavior of such interfaces can be highly non-linear thus increasing the problem complexity. The use of boundary integral methods reduces the dimensionality of the problem. This enables to focus the computational efforts on the fracture plane and allows for a detailed description of the interfacial field's evolution.

# Statement of need

`cRacklet` is a C++ library with a Python interface [@pybind11] developed as a collaboration between the Computational Solid Mechanics Laboratory at EPFL and the Department of Aerospace Engineering of the University of Illinois at Urbana-Champaign.  `cRacklet` implements a spectral formulation of the elastodynamics boundary integral relations between the displacements and the corresponding traction stress acting at a planar interface between two homogeneous elastic solids [@geubelle_spectral_1995], [@breitenfeld_numerical_1998]. The formulation implemented is the *independent* one, which considers the top and bottom solids separately [@breitenfeld_numerical_1998]. The stresses acting on the interfaces are connected to the history of the displacement on the interface in the Fourier domain via a convolution. The convolutions are efficiently computed within a shared-memory parallel framework using FFTW3/OpenMP. The prescription of an interfacial behavior allows solving for the equilibrium of the interface. Time integration is achieved using an explicit time-stepping scheme. `cRacklet` is aimed at researchers interested in interfacial dynamics, ranging from nucleation problems to dynamic propagation of rupture fronts. While the spectral boundary integral formulation is not novel and is already used and cited in the scientific literature, we believe that `cRacklet` will be a useful addition to the community. `cRacklet` is efficient, accessible (C++ or Python), and suited to study a broad class of problems (fracture and friction).

# Features

`cRacklet` allows for planar rupture interface simulations (in 2D or 3D) loaded in any combination of normal traction, in-plane, and out-of-plane shear solicitations. `cRacklet` handles the simulation of interfaces bonded between dissimilar elastic solids. Any stress or material heterogeneity along the fracture plane can be resolved using `cRacklet`. Several interfacial behaviors are included in the library, such as:

- Slip-weakening law [@ida_cohesive_1972] [@palmer_growth_1973]: the cohesive strength is a linearly decreasing function of the opening gap. This law can be coupled with a friction law to handle surface interactions in the case of contact between the solids. Two implementations are available, the classical Coulomb friction law and a regularized one [@prakash_frictional_1998].

- Rate and state dependant friction laws: the frictional resistance is a function of the slip velocity and the history of the interface (the state variable). Several formulations are implemented, including the original ones by [@dieterich_modeling_1979] and [@ruina_slip_1983]. More novel formulations such as rate and state friction with velocity-strengthening behaviors (i.e. N-shaped) are also available, see [@barsinai_2014] for example.

# Performance

We illustrate in \autoref{fig:scalability} the scaling capability of `cRacklet` and compare it to Amdahl's law. The scaling study shows that approximately $80$ to $85$ of the program is parallelized (this includes the computation of the Fourier transform of the displacements, the convolution and the invert transform of the stresses back to the real domain.).
 
![Time required to solve $1e5$ time step with $2^{15}$ discretization points, as a function of the number of threads. Computations run using the computational facilities of EPFL, here on a node composed of 2 Intel Broadwell processors running at $2.6 GHz$ with 14 cores each. The dashed grey lines correspond to Amdahl's law for the theoretical speedup, respectively with $85%$ (upper bound) and $80%$ (lower bound) of the program parallelized.\label{fig:scalability}](scalability.png){ width=80% }

# Example

The successive panels of \autoref{fig:evolution} illustrate the nucleation and propagation of a frictional crack at the interface between two solids with a spatially heterogeneous strength. In \autoref{fig:evolution} (a), the initial configuration of the system is shown: the areas in white are sticking (i.e. no velocity) while colored area are sliding (blue is for low slip velocity and red for larger ones). The shear load is increased with time in the following panels. The slip velocity increases and previously sticking part of the interface starts sliding. In \autoref{fig:evolution} (d), the expansion of frictional cracks has broken almost the entire interface.

![Snapshot of the slip velocity at the interface between two elastic solids under shear loading. The initial strength is highly heterogeneous. Loading and time has increased between the snapshots, starting from (a) to (d). White area correspond to sticking condition (no velocity) while colored ones are sliding. Low velocities are in blue and large ones in red.](scalability.png){ width=80% }

# Publications

The following publications have been made possible with cRacklet:

- @barras_study_2014

- @barras_interplay_2017

- @barras_supershear_2018

- @brener_unstable_2018

- @barras_emergence_2019

- @barras_emergence_2020

- @fekak_crack_2020

- @rezakhani_finite_2020

- @brener_unconventional_2021

- @lebihain_instability_2021

- @roch_velocity_2021

# Acknowledgments

T.R, F. B. and J-F. M. acknowledge the financial support from the Rothschild Caesarea Foundation and from the Swiss National Science Foundation (grants #162569 "Contact mechanics of rough surfaces). P. H. G. acknowledge the financial support from...

# References
