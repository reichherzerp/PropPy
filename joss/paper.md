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

# Introduction




# Statement of need 

Understanding the transport of charged high-energy particles in turbulent magnetic fields is essential for resolving the long-standing question of their extragalactic origin. The transport properties of cosmic rays are relevant in many ways: 

* In cosmic-ray sources, the transport properties determine their residence time in the sources and thus the interaction processes leading to the production of secondary particles [@BeckerTjus2020]. 
* Due to the enormous distance from sources to our Galaxy, cosmic rays have to travel through the turbulent intergalactic medium [@2018arXiv181103062A]. 
* In our Galaxy, the Galactic magnetic field influences their trajectory and, finally, their arrival in the Earth's atmosphere [@Reichherzerb2021].

Analytical theories have been developed over the last century [@Jokipii_1966; @Zweibel2013; @Schlickeiser2015] to describe the transport of cosmic rays. However, these theories are often limited by simplifying assumptions. To overcome these limitations, propagation codes have been developed over the last decades to overcome these limitations with dedicated cosmic-ray-transport simulations [@Giacalone1999; @Casse2001; @Shukurov2017; @Reichherzer2020; @Reichherzer2021b]. In EOM propagation methods, particles are moved stepwise, with the next step always determined based on the Lorentz force. Note the magnetic field must be computed for each propagation step for all particle positions, a process that is typically time-consuming in numerical simulations [@Schlegel2020]. A much more efficient method, the diffusive approach, is based on the statistical properties of the particles and exploits their theoretical description via a transport equation [@CRPropa2017]. In the limit of infinitely large times, a diffusive transport occurs for all charged particles in isotropic turbulence, which can be described by the diffusion tensor. A major drawback of this approach is that can only model the transport of charged particles over large time scales so that the particles have enough time to become diffusive. This is especially relevant for modeling transport in compact sources, where diffusion does not necessarily occur [@Reichherzerp2021].

To tackle this issue and meet the need for realistic and fast simulations of the compact sources of cosmic rays, we present in this paper PropPy. Our software applies the approach of CRW, where statistical aspects are used for speed-up while also providing a good description of the initial phase. Additionally, the properties of the CRW can be determined directly from the diffusion tensor and the gyration radius of the particle.

The comparison of the three different approaches diffusive, EOM, and CRW shows that CRW simulation results are in good agreement with EOM simulations, while being considerably faster.

# Theory

First we assume particle transport in one dimension, where they can move either in positive or negative direction along the $x$-axis. The following derivation has been discussed in various contexts in the literature, such as when describing animal trails (see e.g., @Codling2008 for a review), but can also be applied for cosmic-ray propagation (see e.g., @Seta2019).

During the CRW, the following two substeps are performed in each propagation step:
\begin{enumerate}
    \item Particles that point in positive direction will turn around with the probability $\xi\tau_\mathrm{s}$ and otherwise continue along that direction with the probability $1-\xi \tau_\mathrm{s}$. The same applies for particles that point in negative direction.
    \item The particles move the distance $\chi$ with the speed $ \chi/\tau_\mathrm{s} \equiv v$ along the direction established in the first substep. 
\end{enumerate}
If we divide the particle distribution per position at time $t$ into one distribution in positive direction $\alpha(x,t)$ and one in negative direction $\beta(x,t)$, the following sub distributions result one propagation step later
\begin{equation}\label{eq:alpha_def}
\alpha(x, t + \tau_\mathrm{s}) = (1-\xi \tau_\mathrm{s})\alpha(x-\chi,t) + \xi \tau_\mathrm{s} \beta(x-\chi,t), 
\end{equation}
\begin{equation}\label{eq:beta_def}
\beta(x, t + \tau_\mathrm{s}) = \xi \tau_\mathrm{s}\alpha(x+\chi,t) + (1-\xi \tau_\mathrm{s}) \beta(x+\chi,t).
\end{equation}
As particles move either in positive or negative direction, the total particle distribution yields $f(x,t) = \alpha(x,t)+\beta(x,t)$. For simplicity we define $\alpha(x,t) \equiv \alpha$ and $\beta(x,t) \equiv \beta$. 
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
Inserting \autoref{eq:alpha_plus_beta} into \autoref{eq:alpha_plus_beta_3} yields
\begin{equation}\label{eq:alpha_plus_beta_rewritten}
\frac{1}{2\xi}\frac{\partial^2 (\alpha + \beta)}{\partial t^2} = \frac{v^2}{2\xi} \frac{\partial^2 (\alpha+\beta)}{\partial x^2} - \frac{\partial (\alpha + \beta)}{\partial t}.
\end{equation}
Finally, with $f = \alpha + \beta$, we have
\begin{equation}\label{eq:alpha_plus_beta_rewritten}
\frac{1}{2\xi}\frac{\partial^2 f}{\partial t^2} = \frac{v^2}{2\xi} \frac{\partial^2 f}{\partial x^2} - \frac{\partial f}{\partial t}.
\end{equation}
In fact, when we generalize this approach for three dimensions, assuming local homogeneity, and connecting diffusion coefficients with the CRW parameters, this leads to the telegraph equation
\begin{equation}\label{eq:telegraph}
\frac{\partial f}{\partial t} + \sum_i \tau_i \frac{\partial^2 f}{\partial t^2}= \sum_i \kappa_i \frac{\partial^2 f}{\partial x_i^2}.
\end{equation}
Therefore, the statistics of particles that follow CRW can be described with the telegraph equation and thus agrees with analytical theories of particle transport of cosmic rays [@Litvinenko2015; @Tautz2016].


# Comparison
In principle, the CRW propagation method implemented in PropPy can be applied wherever other propagation codes for charged particles such as CRPropa [@CRPropa2016; @CRPropa2021], DRAGON [@Dragon2017], GALPROP [@Galprop1998] are already in use. However, the advantages of PropPy are especially in the high performance and the accurate description of statistical transport properties also for the initial transport regime, which is not possible for pure diffusive propagation approaches. 



Charged particles (cosmic rays) are accelerated to high energies in astrophysical sources until the gyration radius exceeds the system size according to the Hillas criterion, and the cosmic rays can no longer be confined by the accelerator. Since strong magnetic fields with a significant amount of turbulence typically prevail in these sources, the description of particle propagation in the sources is nontrivial and complicate the analytical description of transport. As an example, we consider the transport of charged particles in AGN jets, an environment for which the above codes are not optimized but whose underlying transport mechanisms can still be applied. 



Simulations are used for describing as accurately as possible the particle transport that has an impact on numerous observable multimessenger signatures. In the following comparison, we focus on the transport properties in these sources, which are described by the diffusion coefficient.



Since CRPropa is the only code that supports both EOM and diffusive propagation methods with anisotropic diffusion coefficients, this software (version: CRPropa 3.1.7) is used for comparison simulations with PropPy. 

While there are numerous possible sources covering a large parameter space of physical properties relevant to particle transport, for this comparison between PropPy propagation and CRPropa modules, we use typical parameters used in the literature for AGN plasmoids (see, e.g. [@Hoerbe2020, @BeckerTjus2022] and references therein):

- particle energies: $E=100\,$PeV
- isotropic 3d Kolmogorov turbulence
- magnetic field strength: $B_\mathrm{rms} = 1\,$Gaus
- correlation length turbulence: $l_\mathrm{c} \sim 10^{11}\,$m


With these parameters, we can derive the expected diffusion coefficient from theory [@Subedi2017]. These parameters result in gyroradii of the charged cosmic rays
\begin{equation}
r_\mathrm{g} = \frac{\sqrt{2}E}{q\,c\,B} = \frac{141\,\mathrm{PeV}}{q\,c \cdot 1\mathrm{G}} \approx 4.72\cdot10^{12}\,\mathrm{m}.
\end{equation}

Note that the factor $\sqrt{2}$ is introduced because of the isotropic directions of the magnetic field vectors in the turbulence. Particles are in the quasi-ballistic transport regime ($r_\mathrm{g} \gg l_\mathrm{c}$), where they experience only minor deflections. The expected diffusion coefficient $\kappa$ is 

\begin{equation}
\kappa_\mathrm{theory} = \frac{r_\mathrm{g}^2 \cdot c}{2l_\mathrm{c}} = \frac{(4.72\cdot10^{12}\,\mathrm{m})^2 \cdot c}{2\cdot 10^{11}\,m} \approx 3.34\cdot10^{22}\,\frac{\mathrm{m^2}}{\mathrm{s}}.
\end{equation}

This theoretical diffusion coefficient serves as an input for the CRPropa SDE and the PropPy simulation, and as a reference for the numerical simulations.

This diffusion coefficient results in expected mean-free paths of
\begin{equation}
\lambda_\mathrm{theory} = \frac{3 \kappa_\mathrm{theory}}{c} \approx 3.34\cdot10^{14}\,\mathrm{m}.
\end{equation}

Particles become diffusive at trajectory lengths of about $\lambda$, which is why the simulations are stopped after trajectory lengths of $10^{17}$ m to have some buffer and a clear plateau in the running diffusion coefficients.

As a simulation setting, $10^3$ protons with $E=100\,$PeV are emitted isotropically from a point source. The simulations and the presented results can be reproduced via the simulation and analysis scripts provided in the comparison folder of PropPy.

The summation of planar waves with different wave numbers, amplitudes, and directions generates the synthetic turbulence. Here, there are two possible approaches:
- The complete turbulence can be generated in advance of the simulation and stored on a large grid by using an inverse discrete Fourier transform. During run-time, the local magnetic field is computed via interpolation of the surrounding grid points that store the magnetic field information. Here, the tri-linear interpolation is used as it is fast and sufficiently accurate [@Schlegel2020]. The turbulence is stored on $1024^3$ grid points.
- The summation of different amplitudes, wavenumbers, and directions can also be performed during run-time at the exact position where it is needed. Numerous constraints of the first method, the grid method, are avoided in this plane-wave (PW) approach, with the disadvantage that the simulations take longer. 1000 wave modes are used, which was determined to be sufficient in convergence tests.


Three different propagation methods implemented in CRPropa are considered:
- Solving EOM with the Boris-Push (BP) method [@CRPropa2021]. 
- Solving EOM with the Cash-Karp (CK) method [@CRPropa2016].
- Solving Stochastic Differential Equations (SDE) [@CRPropa2017]. For this method, no turbulence has to be generated, but only the diffusion coefficient has to be inputted, which already contains the information on how the particles move statistically in the turbulence.

\autoref{fig:comparison} shows a comparison of the simulation results for the calculated running diffusion coefficients for the different methods of propagation and turbulence generation. 


The left and right panels differ only in the simulation length. In the left panel, only trajectories up to $10^{14}$ m are considered, whereas in the right panel, trajectories up to $10^{17}$ m are displayed. Since the mean-free path length indicates the transition between ballistic to diffusive propagation, the left panel shows ballistic particle propagation and the right panel diffusive propagation.


The top panel shows the running diffusion coefficients as a function of time. The middle panel shows the effective diffusion coefficient at $10^{14}$ m on the left and the converged diffusion coefficient on the right, since the running diffusion coefficient remains constant beginning at the diffusive limit, as can be seen in the top panel. 

The lowest panel shows the required processor time of the simulation as a function of the step size. The same processor was used for all simulations for better comparability. 


![Comparison between different propagation approaches for the computation of running diffusion coefficients. $10^3$ protons with $E=10^{15}$ eV simulated in Kolmogorov turbulence generated with the grid and PW method described in the text. \textit{Left panel} shows the ballistic propagation regime at the beginning of particle trajectories ($\leq 10^{14}$ m). \textit{Right panel} shows also long trajectory lengths, for which the particle transport becomes diffusive ($\gg \lambda$). \textit{Upper panel:} running diffusion coefficients as functions of trajectory lengths for different propagation methods. \textit{Middle panel:} diffusion coefficients at $10^{14}$ m (left) and converged ones (right) as functions of different step sizes of the propagation methods. Log deviation is defined as $|\mathrm{log}(\kappa_\mathrm{sim}) - \mathrm{log}(\kappa_\mathrm{theory})|$. \textit{Lower panel:} Simulation time per processor as functions of step sizes. Color-coded with the deviation from theoretical predictions. Usage of 260 threads for step sizes smaller than $10^{13}$ m. Note that the theory prediction is $\kappa_\mathrm{theo} = 3.34\cdot10^{22}\,\mathrm{m^2}/\mathrm{s}$. Averaging over 20 seeds of the same simulation with BP and PW with $s_\mathrm{step} = 10^{10}\,$m gives $\kappa_\mathrm{sim} = (3.41 \pm 0.16)\cdot10^{22}\,\mathrm{m^2}/\mathrm{s}$.
 \label{fig:comparison}](comparison_compact_source.pdf)


The comparisons yield the following results:

- The EOM-based propagation approaches BP and CK, as well as the CRW method from PropPy can correctly model the initial ballistic transport phase. The diffusive approach (SDE) can not describe this initial ballistic phase by construction since it always assumes diffusive particle transport. 

- The diffusive approach and the CRW approach can use relatively large step sizes to model the correct statistical behavior. The latter only needs to resolve the mean-free paths sufficiently well, which is guaranteed if the step size at least ten times smaller than the mean-free path. In the case of the EOM-based methods step sizes must be small enough to resolve both the gyration motion and the scales of turbulence sufficiently well. This can be seen as the diffusion coefficients in the middle right panel converge to a constant value only when the step sizes are smaller than the gyration radii and smaller than the correlation length of the turbulence.

- Smaller simulation times for given step sizes in combination with the fewer step size requirements translate into significant increasion in speed for the diffusive method and the CRW method compared to the EOM-based methods. 


# Conclusion
PropPy is an open-source python software for propagating charged high-energy particles (cosmic rays, CRs) in a turbulent magnetic field. Its modular architecture comprises various modules for sources, magnetic fields, propagators, and observers covering a wide range of applications.

When compared to codes that solve the equation of motion (EOM) in each propagation step, our propagation is based on a correlated random walk (CRW) in Cartesian (for isotropic diffusion) or cylindrical (for anisotropic diffusion) coordinates, which makes each simulation step significantly faster. This novel approach is justified by the fact that a transport equation can be derived via the formulation of the CRW (see theory section below), which is used in analytical descriptions of particle transport [@Litvinenko2015; @Tautz2016]:
\begin{equation}\label{eq:telegraph}
\frac{\partial f}{\partial t} + \sum_i \tau_i \frac{\partial^2 f}{\partial t^2}= \sum_i \kappa_i \frac{\partial^2 f}{\partial x_i^2},
\end{equation}
where $i$ indicates the three spatial directions, $\tau_i$ denotes the time scale for particles to become diffusive, and $\kappa_i$ is the diffusion coefficient, from which the relevant parameters of the CRW can be determined. 

Besides the analytical verification of the CRW ansatz, comparison simulations between PropPy and an established cosmic-ray propagation software, CRPropa, are presented.
These tests show that both approaches are comparable in terms of the statistical properties such as the running diffusion coefficient and the escape times from regions such as are relevant and present in many astrophysical environments.

This makes PropPy a high-performance software for the simulation of charged particles in turbulent magnetic fields. This is especially true for compact objects and transient events with short time scales, such as gamma-ray bursts (GRBs), active galactic nuclei (AGN) flares, where the accurate description of the initial particle propagation is crucial. Fast simulations of transient events can help analyze observations and provide information to evaluate the need for follow-up observations in the context of real-time multimessenger astrophysics [@AstroColibri2021].



# Acknowledgements

PR wants to thank the audience in his [conference contribution](https://indico.cern.ch/event/1037017/contributions/4514419/) on the software and users, who helped with valuable feedback. We thank for funding from the German Science Foundation DFG, within the Collaborative Research Center SFB1491 "Cosmic Interacting Matters - From Source to Signal".
Special thanks to L. Schlegel, F. Schüssler, J. Suc, and E.G. Zweibel for valuable discussions.

# References
