---
title: 'PropPy -- Correlated random walk propagation of cosmic rays in magnetic turbulence'
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
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: J. Becker Tjus
    orcid: 0000-0002-1748-7367
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Ruhr-Universität Bochum, D-44801 Bochum, Germany
   index: 1
 - name: Ruhr Astroparticle and Plasma Physics Center, D-44780 Bochum, Germany
   index: 2
 - name: Université Paris-Saclay, F-91190 Gif-sur-Yvette, France
   index: 3
date: 14 January 2022
bibliography: paper.bib

---
# Summary 
PropPy is an open-source Python software package for propagating charged high-energy particles (cosmic rays, CRs) in a turbulent magnetic field. Its modular architecture comprises various modules for sources, magnetic fields, propagators, and observers covering a wide range of applications.

When compared to codes that solve the equation of motion (EOM) in each propagation step, our propagation is based on a correlated random walk (CRW) in Cartesian (for isotropic diffusion) or cylindrical (for anisotropic diffusion) coordinates, which makes each simulation step significantly faster. This novel approach is justified by the fact that a transport equation can be derived via the formulation of the CRW (see theory section below), which is used in analytical descriptions of particle transport [@Litvinenko2015; @Tautz2016]:
\begin{equation}\label{eq:telegraph}
\frac{\partial f}{\partial t} + \sum_i \tau_i \frac{\partial^2 f}{\partial t^2}= \sum_i \kappa_i \frac{\partial^2 f}{\partial x_i^2},
\end{equation}
where $i$ indicates the three spatial directions, $\tau_i$ denotes the time scale for particles to become diffusive, and $\kappa_i$ is the diffusion coefficient, from which the relevant parameters of the CRW can be determined. 

Besides the analytical verification of the CRW ansatz, comparison simulations between PropPy and an established cosmic-ray propagation software, CRPropa, are presented.
These tests show that both approaches are comparable in terms of the statistical properties such as the running diffusion coefficient and the escape times from regions such as are relevant and present in many astrophysical environments.

This makes PropPy a high-performance software package for the simulation of charged particles in turbulent magnetic fields. This is especially true for compact objects and transient events with short time scales, such as gamma-ray bursts (GRBs), active galactic nuclei (AGN) flares, where the accurate description of the initial particle propagation is crucial. Fast simulations of transient events can help analyze observations and provide information to evaluate the need for follow-up observations in the context of real-time multimessenger astrophysics [@AstroColibri2021].

# Statement of need 
Understanding the transport of charged high-energy particles in turbulent magnetic fields is essential for resolving the long-standing question of their extragalactic origin. The transport properties of cosmic rays are relevant in many ways: 

* In cosmic-ray sources, the transport properties determine their residence time in the sources and thus the interaction processes leading to the production of secondary particles [@BeckerTjus2020]. 
* Due to the enormous distance from sources to our galaxy, cosmic rays have to travel through the turbulent intergalactic medium [@2018arXiv181103062A; @Schlegel:2020]. 
* In our galaxy, the galactic magnetic field influences their trajectory and, finally, their arrival in the Earth's atmosphere [@Reichherzerb2021].

Analytical theories have been developed over the last century [@Jokipii_1966; @Zweibel2013; @Schlickeiser2015; @Shalchi:2021ApJ] to describe the transport of cosmic rays. However, these theories are often limited by strongly simplifying assumptions concerning the transport of charged particles in turbulent magnetic fields. To overcome these limitations, propagation codes with dedicated cosmic-ray-transport simulations have been developed over the last decades [@Giacalone1999; @Casse2001; @Shukurov2017; @Reichherzer2020; @Reichherzer2021b]. In EOM propagation methods, particles are moved stepwise, with the next step always determined based on the solution of the EOM with the external force as the Lorentz force only taking into account magnetic fields. Note the magnetic field must be computed for each propagation step for all particle positions, a process that is typically time-consuming in numerical simulations. This is especially relevant when the particles are highly diffusive, i.e.\,when the size of the propagation environment $L$ exceeds the gyro radius of the particle $r_g\ll L$. A much more efficient method, the diffusive approach, is based on the statistical properties of the particles and exploits their theoretical description via a transport equation [@CRPropa2017]. In the limit of infinitely large times, diffusive transport occurs for all charged particles in isotropic turbulence. In the transport equation, the diffusion tensor implicitly contains all statistic properties. A major drawback of this approach is that can only model the transport of charged particles over large time scales so that the particles have enough time to become diffusive [@BeckerTjus2022].

To tackle this issue and meet the need for realistic and fast simulations of the sources of cosmic rays, we present the PropPy software. Our software applies the approach of the CRW, where statistical aspects are used for speed-up while also providing a good description of the initial phase. Additionally, the properties of the CRW can be determined directly from the diffusion tensor and the gyration radius of the particle.

In principle, the CRW propagation method implemented in PropPy can be applied wherever other propagation codes for charged particles such as CRPropa [@CRPropa2016; @CRPropa2021], DRAGON [@Dragon2017], GALPROP [@Galprop1998] are already in use. However, the advantages of PropPy are especially in the high performance and the accurate description of statistical transport properties also for the initial transport regime, which is not possible for pure diffusive propagation approaches. 


# Comparison
Simulations are used for describing as accurately as possible the particle transport that has an impact on numerous observable multimessenger signatures. In the following comparison, we focus on the transport properties in these sources, which are described by the diffusion coefficient.

Since CRPropa is the only code that supports both EOM and diffusive propagation methods with anisotropic diffusion coefficients, this software (version: CRPropa 3.1.7) is used for comparison simulations with PropPy. 

We compare the performance of PropPy with the two different propagation methods implemented in CRPropa, which are:

\begin{enumerate}
  \item Solving the EOM, using either the Boris-Push (BP) or the Cash-Karp (CK) algorithm. 
  \item Solving Stochastic Differential Equations (SDE). For this method, no turbulence has to be generated, but only the diffusion coefficient has to be inputted, which already contains the information on how the particles move statistically in the turbulence.
\end{enumerate}

The comparison of the three different approaches diffusive, EOM, and CRW shows that CRW simulation results are in good agreement with EOM simulations, while being considerably faster.


# Acknowledgements

We acknowledge support from funding from the German Science Foundation DFG, within the Collaborative Research Center SFB1491 "Cosmic Interacting Matters - From Source to Signal".
Special thanks to L. Schlegel, F. Schüssler, J. Suc, and E.G. Zweibel for valuable discussions.

# References
