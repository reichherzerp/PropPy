---
title: 'Correlated random walk propagation of cosmic rays in magnetic turbulence'
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
  - name: J. Becker Tjus #^[first author] # note this makes a footnote saying 'first author'
    orcid: 0000-0003-4513-8241
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

PropPy is an open-source python software for propagating charged high-energy particles (cosmic rays) in a turbulent magnetic field. Its modular architecture comprises various sources, magnetic fields, propagators, observers, and analyzing modules covering a wide range of applications.

Here, propagation is based on a correlated random walk (CRW) in Cartesian (for isotropic diffusion) or cylindrical (for anisotropic diffusion) coordinates, which makes each simulation step significantly faster than comparable codes that have to solve the equation of motion (EOM) in each propagation time. This novel approach is justified by the fact that a transport equation can be derived via the formulation of the CRW (see theory section below), which is used in analytical descriptions of particle transport [@Litvinenko2015; @Tautz2016]
\begin{equation}\label{eq:telegraph}
\frac{\partial f}{\partial t} + \sum_i \tau_i \frac{\partial^2 f}{\partial t^2}= \sum_i \kappa_i \frac{\partial^2 f}{\partial x_i^2}.
\end{equation}
$i$ indicates the three spatial directions, and $\tau_i$ denotes the time scale for particles to become diffusive. From the diffusion coefficients $\kappa_i$, the relevant parameters of the CRW can be determined. These considerations confirm the validity for the diffusive phase of particle transport that the CRW approach correctly describes. The initial propagation phase is not yet diffusive and is therefore only tested by comparison tests between CRW and EOM that show that statistical properties such as the running diffusion coefficient $\kappa_i(t) = \langle x_i^2\rangle /2t$ and the escape times from regions such as are relevant and present in many astronomical environments are comparable in both approaches.

This makes PropPy a high-performance software for the simulation of charged particles in turbulent magnetic fields, especially for compact objects and transient events with short time scales, such as gamma-ray bursts (GRBs), active galactic nuclei (AGN) flares, where the accurate description of the initial particle propagation is crucial. Fast simulations of transient events can help to analyze observations and provide information to evaluate the need for follow-up observations in the context of real-time multimessenger astrophysics [@AstroColibri2021].


# Statement of need 

Understanding the transport of charged high-energy particles in turbulent magnetic fields is essential for resolving the long-standing question of their extragalactic origin. The transport properties of cosmic rays are relevant in many ways to understand their origin: 

* In their sources, the transport properties determine the residence time in the sources and thus the interaction processes leading to the production of secondary particles [@BeckerTjus2020]. 
* Due to the enormous distance from sources to our galaxy, cosmic rays have to travel through the turbulent intergalactic medium [@2018arXiv181103062A]. 
* In our Galaxy, the Galactic magnetic field influences their trajectory and, finally, their arrival in the Earth's atmosphere [@Reichherzerb2021].

Analytical theories have been developed over the last century [@Jokipii_1966; @Zweibel2013; @Schlickeiser2015] to describe the transport of cosmic rays. However, these theories are often limited by simplifying assumptions; hence propagation codes have been developed over the last decades to overcome these limitations [@Giacalone1999; @Casse2001; @Shukurov2017; @Reichherzer2020]. In the most basic propagation method, particles are moved stepwise, with the next step always determined based on the Lorentz force. For this, however, the magnetic field must be computed for each propagation step for all particle positions, a process that is typically time-consuming in numerical simulations [@Schlegel2020]. A much more efficient method (diffusive approach) is based on the statistical properties of the particles and exploits their theoretical description via a transport equation [@CRPropa2017]. In the limit of infinitely large times, a diffusive transport occurs for all charged particles in isotropic turbulence, which can be described by the diffusion tensor. 

However, by definition, this approach can only model the transport of charged particles over large time scales so that the particles have enough time to become diffusive, which is a major drawback, especially when modeling transport in compact sources, where diffusion does not necessarily occur [@Reichherzerp2021].

PropPy was established to tackle this issue and meet the need for realistic and fast simulations of the compact sources of cosmic rays, such as GRB and AGN flares. With this publication, we present PropPy, which applies the novel approach of CRW, where statistical aspects are used for speed-up while also providing a good description of the initial phase, where classic diffusion approaches fail. Furthermore, the properties of the CRW can be determined directly from the diffusion tensor and the gyration radius of the particle.

The comparison of the three different approaches illustrates the good agreement of the CRW with the solution of the equation of motion. 


# Theory

Let us generally assume particle transport in one dimension. The following derivation is discussed in other contexts in the literature, such as when describing animal trails (see e.g., @Codling2008 for a review).

During the CRW, the following two substeps are performed in each propagation step:
\begin{enumerate}
    \item Particles that point in direction $x$ will turn around with the probability $\xi\tau_\mathrm{s}$ and otherwise continue along that direction with the probability $1-\xi \tau_\mathrm{s}$.
    \item The particles move the distance $\chi$ with the speed $ \chi/\tau_\mathrm{s} \equiv v$ along the direction established in the 1. substep. 
\end{enumerate}
If we divide the particle distribution per position at time $t$ into one distribution in positive direction $\alpha(x,t)$ and one in negative direction $\beta(x,t)$, the following sub distributions result one propagation step later
\begin{equation}\label{eq:alpha_def}
\alpha(x, t + \tau_\mathrm{s}) = (1-\xi \tau_\mathrm{s})\alpha(x-\chi,t) + \xi \tau_\mathrm{s} \beta(x-\chi,t), 
\end{equation}
\begin{equation}\label{eq:beta_def}
\beta(x, t + \tau_\mathrm{s}) = \xi \tau_\mathrm{s}\alpha(x+\chi,t) + (1-\xi \tau_\mathrm{s}) \beta(x+\chi,t).
\end{equation}
As particles move either in positive or negative direction, it yields $f(x,t) = \alpha(x,t)+\beta(x,t)$. For simplicity we define $\alpha(x,t) \equiv \alpha$ and $\beta(x,t) \equiv \beta$. 
Expanding \autoref{eq:alpha_def} and \autoref{eq:beta_def} for small steps $\tau_s, \chi \xrightarrow{} 0$ yields
\begin{equation}\label{eq:alpha}
\frac{\partial \alpha}{\partial t} = - v \frac{\partial \alpha}{\partial x} + \xi(\beta-\alpha),
\end{equation}
\begin{equation}\label{eq:beta}
\frac{\partial \beta}{\partial t} = v \frac{\partial \beta}{\partial x} - \xi(\beta-\alpha).
\end{equation}
Adding component-wise \autoref{eq:alpha} and \autoref{eq:beta} yields
\begin{equation}\label{eq:alpha_plus_beta}
\frac{\partial (\alpha + \beta)}{\partial t} = v \frac{\partial (\beta-\alpha)}{\partial x},
\end{equation}
with the time derivative
\begin{equation}\label{eq:alpha_plus_beta_deriv}
\frac{\partial^2 (\alpha + \beta)}{\partial t^2} = v \frac{\partial^2 (\beta-\alpha)}{\partial t \, \partial x}.
\end{equation}
Substracting component-wise \autoref{eq:alpha} from \autoref{eq:beta} and derivating with respect to $x$ yields
\begin{equation}\label{eq:beta_minus_alpha}
\frac{\partial^2 (\beta - \alpha)}{\partial t\, \partial x} = v\frac{\partial^2 (\alpha+\beta)}{\partial x^2} - 2\xi \frac{\partial (\beta - \alpha)}{\partial x}.
\end{equation}
Inserting \autoref{eq:beta_minus_alpha} into \autoref{eq:alpha_plus_beta_deriv} results in
\begin{equation}\label{eq:alpha_plus_beta_3}
\frac{\partial^2 (\alpha + \beta)}{\partial t^2} = v^2 \frac{\partial^2 (\alpha+\beta)}{\partial x^2} - 2\xi v\frac{\partial (\beta - \alpha)}{\partial x}.
\end{equation}
Inserting \autoref{eq:alpha_plus_beta} into \autoref{eq:alpha_plus_beta_3}
\begin{equation}\label{eq:alpha_plus_beta_rewritten}
\frac{1}{2\xi}\frac{\partial^2 (\alpha + \beta)}{\partial t^2} = \frac{v^2}{2\xi} \frac{\partial^2 (\alpha+\beta)}{\partial x^2} - \frac{\partial (\alpha + \beta)}{\partial t},
\end{equation}
Finally, with $f = \alpha + \beta$, we have
\begin{equation}\label{eq:alpha_plus_beta_rewritten}
\frac{1}{2\xi}\frac{\partial^2 f}{\partial t^2} = \frac{v^2}{2\xi} \frac{\partial^2 f}{\partial x^2} - \frac{\partial f}{\partial t}.
\end{equation}
Generalizing this approach for three dimensions, assuming local homogeneity, and connecting diffusion coefficients with the CRW parameters leads to \autoref{eq:telegraph}. Therefore, the statistics of particles that follow CRW can be described with the telegraph equation and thus agrees with analytical theories of particle transport of cosmic rays [@Litvinenko2015; @Tautz2016].


# Comparison

In principle, PropPy can be applied wherever other propagation codes for charged particles such as CRPropa [@CRPropa2016; @CRPropa2021], DRAGON [@Dragon2017], GALPROP [@Galprop1998] are already in use. However, the advantages of PropPy are especially in the improved performance and the accurate description of statistical transport properties also for the initial propagation, which is not possible for pure diffusive propagation approaches. 

For example, the propagation of charged particles in blobs of blazar jets is discussed in the following. Due to the high performance and the good statistical description, even at early times, the software is excellently suited for calculating escape times of charged particles from certain zones (blob in the example), which in turn are required in (semi)analytical calculations. Also, particle distributions and arrival times can be simulated efficiently. 

Since CRPropa is the only code that supports both EOM and diffusive with anisotropic diffusion coefficients, this software is used for comparison simulations with the CRW approach. As a simulation setting, we consider $10^3$ protons isotropically emitted from a point source with an energy of $3\cdot 10^{15}$ eV. 


The synthetic turbulence is generated as a superposition of planar waves of different amplitudes, wave numbers, and direction. Here, there are two possible approaches:
1. The complete turbulence can be generated in advance of the simulation and stored on a large grid by using inverse discrete Fourier transform. During run-time, the local magnetic field is computed via a interpolation of the surrounding grid points that store the magnetic field information. Here, we use tri-linear interpolation as it is fast and sufficiently accurate [@Schlegel2020]. We store the turbulence on $1024^3$ grid points.
2. The summation of different amplitudes, wave numbers, and directions can also be performed during run-time at the exact position where it is needed. Numerous constraints of the grid method are avoided here, with the disadvantage that the simulations take longer. We use 1000 wave modes, which was determined to be sufficient in convergence tests.

We consider three different propagation methods implemented in CRPropa:

1. Solving EOM with the Boris-Push (BP) method [@CRPropa2021]. 
2. Solving EOM with the Cash-Karp (CK) method [@CRPropa2016].
3. Solving Stochastic Differential Equations (SDE) [@CRPropa2017]. For this method, no turbulence has to be generated, but only the diffusion coefficient has to be passed, which already contains the information how the particles move statistically in the turbulence.

\autoref{fig:comparison} shows a comparison of the simulation results for the calculated running diffusion coefficients for the different propagation methods and turbulence generation methods. 

The left and right panels differ only in the simulation length considered. In the left panel only trajectories up to $10^{14}$ m are considered, whereas in the right panel trajectories up to $10^{17}$ m are considered. Since the mean free path length indicates the transition between ballistic to diffusive propagation, the left panel shows ballistic particle propagation and the right panel diffusive propagation.

The top panel shows the running diffusion coefficients as a function of time. The middle panel shows the effective diffusion coefficient at $10^{14$} m on the left and the converged diffusion coefficient on the right, since the running diffusion coefficient remains constant beginning at the diffusive limit, as can be seen in the top panel. 
The lowest panel shows the required processor time of the simulation as a function of the step size. The same processor was used for all simulations for better comparability. 


![Comparison between different propagation approaches for the computation of running diffusion coefficients. $10^3$ protons with $E=10^{15}$ eV simulated in magnetic field configurations described in the text. Configuration 1 without the ordered background magnetic field is used for the left panel and configuration 2 for the right panel. Detailed explanations and the code, respectively, the data for the reproducibility are available in tutorial 4 of the public repository. \label{fig:comparison}](comparison_compact_source.pdf)

The comparisons yield the following conclusions:
- Both EOM-based propagation approaches BP and CK and the correlated random walk from PropPy can correctly model the initial ballistic propagation phase. The Diffusive approach (SDE) cannot by construction, since it assumes diffusive particle transport. 
- The diffusive approach and the correlated random walk approach can use relatively large step sizes to model the correct statistical behavior. The latter only needs to resolve the mean-free paths sufficiently well, which is given by choosing the step size approximately 10 times smaller than the mean-free path. The EOM-based methods must choose step sizes small enough to resolve both the gyration motion and the scales of tubulence sufficiently well. This can be seen as the diffusion coefficients in the middle panel in the right figure converge to a constant value only when the step size is smaller than the gyration radii and smaller than the correlation lengths.
- The combination of smaller simulation times at the same step sizes with the lower step size requirements, by which larger step sizes lead to comparable results, translates into huge speedups for the diffusive method and the correlated random walk method implemented here in PropPy compared to the EOM-based methods. 


# Acknowledgements

I want to thank the audience in my [conference contribution](https://indico.cern.ch/event/1037017/contributions/4514419/) on the software and users, who helped with valuable feedback. 
A special thanks to L. Schlegel, F. Schüssler, and E.G. Zweibel for valuable discussions.


# References




