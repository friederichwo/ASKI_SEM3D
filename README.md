# ASKI_SEM3D_GEMINI

[![DOI](https://zenodo.org/badge/842171779.svg)](https://zenodo.org/doi/10.5281/zenodo.13321935)

This software performs teleseismic full waveform inversion on a spherical earth based on waveform sensitivity kernels using SPECFEM3D Cartesian as forward solver and GEMINI for remote wavefield injection. It is derived from the original ASKI package but specialized to SPECFEM as forward method, simplified, refurbished and siginificantly extended to support hybrid teleseismic full waveform inversion with remote wavefield injection into the regional SPECFEM domain.

The code uses teleseismic wavefields computed in a hybrid way with GEMINI and SPECFEM for a set of earthquakes as well as wavefields propgated away from the receiver positions to compute waveform sensitivity kernels for all available combinations of earthquakes and receivers. For efficiency reasons, the wavefields are  highly compressed by converting them into Fourier coefficients at a small set of frequencies. Kernels are computed as simple products in the frequency domain and pre-integrated over one SPECFEM element (``computeKernelsDmspaceParallel``). It also calculates Fourier coefficients of observed and synthetic seismograms at the receivers, the latter obtained from a SPECFEM forward run (``tranformMeasuredDataAsciiSeis`` and ``transformSpecfemSynthetics``). Data, synthetics and kernels are fed into a big linear system of equations extended by model regularization equations that is solved in a least-squares sense using a conjugate-gradient algorithm for velocity perturbations in the subsurface (``solveCglsKernelSystem``). The velocity model is updated and fed into SPECFEM for another iteration. 

Besides these core programs, there are many useful Python tools in the folder ``askipy`` that help in mastering the challenging workflow of full waveform inversion and also allow displaying velocity perturbations, kernels, and wavefield. The complex interplay between ASKI, SPECFEM and GEMINI is described in an extensive, 50 pages LaTeX documentation (see the ``doc`` folder) that provides step-by-step instructions how to use the code.

------------
Contributors
------------
Wolfgang Friederich, Florian Schumacher (now Schorn), Marcel Paffrath, Thomas Möller, Kasper D. Fischer, Samir Lamara

----------
References
----------
Schumacher, F., Friederich, W. and S. Lamara. A flexible, extendable, modular and computationally efficient approach to scattering-integral-based seismic full waveform inversion, Geophysical Journal International, 204, 2, 1100-1119, 2016, https://doi.org/10.1093/gji/ggv505.

Schumacher, F. and W. Friederich. ASKI: a modular toolbox for scattering-integral-based seismic full waveform inversion and sensitivity analysis utilizing external forward codes, Software X, 2016,  http://dx.doi.org/10.1016/j.softx.2016.10.005.

Friederich, W. and J. Dalkolmo. Complete synthetic seismograms for a spherically symmetric earth by a numerical computation of Green's function in the frequency domain,	Geophys. J. Int., 122, 537--550, 1995, https://doi.org/10.1111/j.1365-246X.1995.tb07012.x.

Komatitsch, D., Erlebacher, G., Göddeke, D. and D. Mich\'ea. High-order finite-element seismic wave propagation modeling with MPI on a large GPU cluster, J. Comput. Phys., 229, 28, 7692-7714, 2010, https://doi.org/10.1016/j.jcp.2010.06.024.

Komatitsch, D. and J. P. Vilotte. The spectral-element method: an efficient tool to simulate the seismic response of 2D and 3D geological structures, Bull. seism. Soc. Am., 88, 2, 368-392, 1998.

Komatitsch, D. and J. Tromp. Introduction to the spectral-element method for 3-{D} seismic wave propagation, Geophys. J. Int., 139, 3, 806-822, 1999, https://doi.org/10.1046/j.1365-246x.1999.00967.x.

Monteiller, V., Chevrot, S., Komatitsch, D. N. Fuji. A hybrid method to compute short-period synthetic seismograms of teleseismic body waves in a 3D regional model,
Geophysical Journal International, 192,1 230-247, 2013, https://doi.org/10.1093/gji/ggs006.
