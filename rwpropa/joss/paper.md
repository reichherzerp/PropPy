---
title: 'Correlated random walk propagation of cosmic rays in turbulence'
tags:
  - Python
  - astronomy
  - cosmic rays
  - transients
  - AGN
  - transport
authors:
  - name: P. Reichherzer^[first author] # note this makes a footnote saying 'first author'
    orcid: 0000-0003-4513-8241
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Ruhr-Universität Bochum, Universitätsstraße 150, 44801 Bochum, Germany
   index: 1
 - name: Université Paris-Saclay, F-91190 Gif-sur-Yvette, France
   index: 2
date: 20 September 2021
bibliography: paper.bib

---

# Summary

RWPropa is an open-source python software for propagating charged high-energy particles (cosmic rays) in a turbulent magnetic field. Its modular architecture comprises various sources, magnetic fields, propagators, observers, and analyzing modules that cover a wide range of applications.

Here, propagation is based on a correlated random walk (CRW), which makes each simulation step significantly faster than comparable codes that have to solve the equation of motion (EOM) each time. This novel approach is justified by the fact that a transport equation can be derived via the formulation of the CRW, which is used in analytical descriptions of particle transport [@Litvinenko:2015; @Tautz:2016]
$$
\frac{\partial f}{\partial t} + \sum_i \tau_i \frac{\partial^2 f}{\partial t^2}= \sum_i \kappa_i \frac{\partial^2 f}{\partial x_i^2}.
$$
$i$ indicates the three spatial directions, and $\tau_i$ denotes the time scale for particles to become diffusive. From the diffusion coefficients $\kappa_i$, the relevant parameters of the CRW can be determined. These considerations confirm the validity for the diffusive phase of particle transport that the CRW approach correctly describes. The initial propagation phase is not yet diffusive and is therefore only tested by comparison tests between CRW and EOM that show that statistical properties such as the running diffusion coefficient $\kappa_i(t) = <x_i^2>/2t$ and the escape times from regions such as are relevant and present in many astronomical environments are comparable in both approaches.

This makes RWPropa a high-performance software for the simulation of charged particles in turbulent magnetic fields, especially for compact objects and transient events with short time scales, such as gamma-ray bursts (GRBs), active galactic nuclei (AGN) flares, where the accurate description of the initial particle propagation is crucial. Fast simulations of transient events can help to analyze observations and provide information to evaluate the need for follow-up observations in the context of real-time multimessenger astrophysics [@AstroColibri:2021].

# Statement of need 

Understanding the transport of charged high-energy particles in turbulent magnetic fields is an essential component of resolving the long-standing question of their extragalactic origin. The transport properties of cosmic rays are relevant in many ways to understand their origin: 

* In their sources, the transport properties determine the residence time in the sources and thus the interaction processes leading to the production of secondary particles [@BeckerTjus:2020]. 
* Due to the enormous distance from sources to our galaxy, cosmic rays have to travel through the turbulent intergalactic medium [@2018arXiv181103062A]. 
* In our Galaxy, the Galactic magnetic field influences their trajectory and, finally, their arrival in the Earth's atmosphere [@Reichherzerb:2021].

To describe the transport of cosmic rays, analytical theories have been developed over the last century [@Jokipii:1966; @Zweibel:2013; @Schlickeiser:2015]. However, these theories are often limited by simplifying assumptions; hence propagation codes have been developed over the last decades to overcome these limitations [@Giacalone:1999; @Casse:2001; @Shukurov:2017; @Reichherzer:2020]. In the most basic propagation method, particles are moved stepwise, with the next step always determined based on the Lorentz force. For this, however, the magnetic field must be computed for each propagation step for all particle positions, a process that is typically time-consuming in numerical simulations [@Schlegel:2020]. A much more efficient method (diffusive approach) is based on the statistical properties of the particles and exploits their theoretical description via a transport equation [@CRPropa:2017]. In the limit of infinitely large times, a diffusive transport occurs for all charged particles in isotropic turbulence, which can be described by the diffusion tensor. 

However, by definition, this approach can only model the transport of charged particles over large time scales so that the particles have enough time to become diffusive, which is a major drawback, especially when modeling transport in compact sources, where diffusion does not necessarily occur [@Reichherzerp:2021].

RWPropa was established to tackle this issue and to meet the need for realistic and fast simulations of the compact sources of cosmic rays, such as GRB, AGN flares. With this publication, we present RWPropa, which applies the novel approach of CRW, where statistical aspects are used for speed-up while also providing a good description of the initial phase, where classic diffusion approaches fail. The properties of the CRW can be determined directly from the diffusion tensor and the gyration radius of the particle.

The comparison of the three different approaches illustrates the good agreement of the CRW with the solution of the equation of motion. 

# Comparison

In principle, RWPropa can be applied wherever other propagation codes for charged particles such as CRPropa [@CRPropa:2016; @CRPropa:2021], DRAGON [@Dragon:2017], GALPROP [@Galprop:1998] are already in use. However, the advantages of RWPropa are especially in the improved performance and the accurate description of statistical transport properties also for the initial propagation, which is not possible for pure diffusive propagation approaches. 

As an example, the propagation of charged particles in blobs of blazar jets is discussed in the following. Due to the high performance and the good statistical description, even at early times, the software is excellently suited for the calculation of escape times of charged particles from certain zones (blob in the example), which in turn are required in (semi)analytical calculations. Also, particle distributions and arrival times can be simulated efficiently. 

Since CRPropa is the only code that supports both EOM and diffusive with anisotropic diffusion coefficients, this software is used for comparison simulations with the CRW approach. As a simulation setting, we consider $10^3$ protons isotropically emitted from a point source with an energy of $3\cdot 10^{15}$ eV. We consider two different magnetic field configurations in the following:

1. only an isotropic tubular magnetic field with magnetic field strength $B = 1$ Gaus. A Kolmogorov spectrum is assumed for the turbulence since this is a good description in many astronomical environments.
2. in addition to the turbulent field from 1. a directed background magnetic field in the $z$-direction is added. This provides the anisotropic diffusion.

\autoref{fig:comparison} shows a comparison of the simulation results for the calculated running diffusion coefficients for the three different propagation methods. 


![Comparison between different propagation approaches for the computation of running diffusion coefficients. $10^3$ protons with $E=10^{15}$ eV simulated in magnetic field configurations described in the text. Detailed explanations and the code respectively the data for the reproducibility are available in tutorial \footnote{\url{https://gitlab.ruhr-uni-bochum.de/reichp2y/rwpropa/-/blob/master/rwpropa/tutorials/Tutorial_4---Tutorial_Comparison_CRPropa.ipynb}}. \label{fig:comparison}](figure_comparison.pdf)

The diffusion approach employs a static diffusion coefficient which, by definition, is only valid in the limit of large times, as \autoref{fig:comparison} illustrates. This incorrect description of the diffusive approach of the initial propagation leads to wrong escape times of the particles from the sources region. Consequently, the number of secondary particles is significantly underestimated. The CRW approach has, despite its simplification and the associated performance improvement in comparison to the EOM approach, a comparable accuracy in the description of the relevant physics.

# References




