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

The study of dynamically propagating rupture along interfaces is of prime importance in various fields and system sizes, including tribology ($nm$ to $\mu m$), engineering ($mm$ to $m$) and geophysics ($m$ to $km$). Numerical simulations of these phenomena are computationally costly and challenging, as they usually require the coupling of two different spatio-temporal scales. A fine spatial discretization is needed to represent accurately the singular fields associated with the rupture edges. Besides, the problems of interest usually involve a larger length scale along which rupture will propagate driven by long-range traveling elastic waves. The physical phenomena at play also occur at different timescales, from the slow process of rupture nucleation to the fast propagation of crack front close the elastic wave speeds. Large and finely discretized spatio-temporal domains are required, which are computationally costly. In addition, the behavior of such interfaces can be highly non-linear thus increasing the problem complexity. The use of boundary integral methods reduces the dimensionality of the problem. This enables to focus the computational efforts on the fracture plane and allows for a detailed description of the interfacial failure processes.

# Statement of need

`cRacklet` is a C++ library with a Python interface [@pybind11] initiated as a collaboration between the Computational Solid Mechanics Laboratory at EPFL and the Department of Aerospace Engineering of the University of Illinois at Urbana-Champaign.  `cRacklet` implements a spectral formulation of the elastodynamics boundary integral relations between the displacements and the corresponding traction stress acting at a planar interface between two homogeneous elastic solids [@geubelle_spectral_1995], [@breitenfeld_numerical_1998]. The formulation implemented is the *independent* one, which considers the top and bottom solids separately [@breitenfeld_numerical_1998]. The stresses acting at the interface are related to the history of interfacial displacements via a time convolution evaluated in the Fourier domain. The convolutions are efficiently computed within a shared-memory parallel framework using FFTW3/OpenMP. The prescription of an interfacial behavior allows for solving the continuity of tractions and displacements through the interface. Time integration is achieved using an explicit time-stepping scheme. `cRacklet` is aimed at researchers interested in interfacial dynamics, ranging from nucleation problems to dynamic propagation of rupture fronts. While the spectral boundary integral formulation is a well-established method that has been extensively referenced in the literature, we believe that `cRacklet` will be a useful addition to the community by gathering in the same framework various kinds of interfacial problems and constitutive laws. `cRacklet` is efficient, accessible (C++ or Python), and suited to study a broad class of problems (fracture and friction). We wish that `cRacklet` will become a link between model developers and users by providing both adaptability and usability.

# Features

1. `cRacklet` is versatile and can be used to study a broad class of problems focused on the behavior of an interface between two semi-infinite solids. The code is particularly suited to study planar dynamic fracture and friction. The interface can be either between two or three-dimensional solids. It can be loaded in any combination of normal traction, in-plane, and out-of-plane shear solicitations. `cRacklet` handles the simulation of interfaces bonded between dissimilar elastic solids. Any stress or material heterogeneity along the fracture plane can be resolved. Several interfacial behaviors are included in the library, such as:

   - Cohesive fracture law [@dugdale_yielding_1960] [@barenblatt_mathematical_1962]: the cohesive strength is a linearly decreasing function of the opening gap. This law can be coupled with a friction law to handle surface interactions in the case of post-failure contact between the solids. Two implementations are available, the classical Coulomb friction law and a regularized one [@prakash_frictional_1998].

   - Rate and state dependant friction laws: the frictional resistance is a function of the slip velocity and the history of the interface (the state variable). Several formulations are implemented, including the original ones by [@dieterich_modeling_1979] and [@ruina_slip_1983]. More novel formulations such as rate and state friction with velocity-strengthening behaviors (i.e. N-shaped) are also available, see [@barsinai_2014] for example.

2. `cRacklet` is accessible and adaptable. It provides access through both its C++ and Python API to several options to design the various kind of problems mentioned before. `cRacklet` is adaptable due to its object-oriented implementation: it is simple to implement additional behavior for the interface without having to deal with the technical core of the code that handles the computation of the stresses in the Fourier domain. `cRacklet` can also be loaded as an external library to easily interact with other existing computational software. `cRacklet` also has tutorials available on Binder [@binder] which allows for a quick and easy introduction to its functionality.

3. `cRacklet` is efficient: the Fourier transforms and the convolutions are computed within a shared-memory parallel framework using FFTW3/OpenMP. We illustrate in \autoref{fig:scalability} the scaling capability of `cRacklet` and compare it to Amdahl's law. The scaling study shows that approximately $85\%$ to $90\%$ of the program is parallelized: this includes the computation of the Fourier transform of the displacements, the convolution, and the invert transform of the stresses back to the real domain.

![Time required to solve $1\text{e}5$ time steps with $2^{12}$ discretization points, as a function of the number of threads. Computations run using the computational facilities of EPFL, here on a node composed of 2 Intel Broadwell processors running at $2.6\,\text{GHz}$ with 14 cores each. The dashed grey lines correspond to Amdahl's law for the theoretical speedup, respectively with $90\%$ (upper bound) and $85\%$ (lower bound) of the program parallelized. \label{fig:scalability}](scalability.png){ width=100% }

# Example

The onset of sliding between two rough surfaces in frictional contact is an illustrative example of a multiscale rupture problem. Macroscopic shearing is resisted by the microcontacts, i.e. by the sparse contacting junctions existing between the asperities of the two surfaces.

The successive panels of \autoref{fig:evolution} illustrate the nucleation and propagation of a frictional rupture at the interface between two solids, from the individual failure of the microcontacts in pannel (b) to the propagation of a macroscopic circular rupture in panel (d). The spatially heterogeneous strength used in this example is a representation of the heterogeneous map of contact between two rough surfaces. In \autoref{fig:evolution} (a), the initial configuration of the system is shown. The areas in white are sticking (i.e. no velocity) and correspond to asperities in contact. Colored areas are sliding (blue is for low slip velocity and red for larger ones). The shear load is increased with time in the following panels. The slip velocities increase and previously sticking parts of the interface start sliding (micro-contacts are broken). The inset of \autoref{fig:evolution} (b) is a zoomed view of the interface where rupture starts at the asperity scale. In \autoref{fig:evolution} (d), frictional cracks have expanded over almost the entire interface.

![Snapshot of the slip velocity at the interface between two elastic solids under shear loading. The initial strength is highly heterogeneous. Loading and time have increased between the snapshots, starting from (a) to (d). White areas correspond to sticking conditions (no velocity) while colored ones are sliding. Low velocities are in blue and large ones in red. This simulation involve $2^{24}$ points and was run on one node (with two 16-core Intel E5-2683v4 2.1 GHz and 512 GiB RAM) of the computing cluster \textit{Fram} from the Norwegian e-infrastructure for research and education. \label{fig:evolution}](evolution.png){ width=95% }

# Publications

The following publications have been made possible with `cRacklet`:

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

T.R, F. B., and J-F. M. acknowledge the financial support from the Swiss National Science Foundation (grants #162569 "Contact mechanics of rough surfaces) and from the Rothschild Caesarea Foundation. F.B. acknowledges  support  of  the  Swiss  National  Science  Foundation through the fellowship No. P2ELP2/188034. F.B. acknowledges the Norwegian e-infrastructure for research and education (UNINETT Sigma2) for computing resources through grant NN9814K.

# References
