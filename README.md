# SphericalFRI

This repository contains a collection of Matlab routines to accompany the paper
[_Sampling Sparse Signals on the Sphere: Algorithms and
Applications_](http://arxiv.org/abs/1502.07577) by Ivan Dokmanić and Yue M. Lu

__Dependecies__: This code requires some routines from the Spherical Harmonic
Transform package SSHT by Jason McEwen and Yves Wiaux. SSSH can be obtained at
[this link](http://www.jasonmcewen.org/codes/ssht/). Most of these dependencies will be removed in the future.

## Main functions

| Filename          | Description                                                    |
| ----------------- | -------------------                                            |
| SphereFRI.m       | Runs the FRI estimation algorithm                              |
| DiracSpectrum.m   | Computes the spectrum of a collection of Diracs                |
| RandomDiracs.m    | Creates a random collection of spherical Diracs                |
| SphereFRI_SN.m    | Estimation routine for the shot noise removal                  |
| SelectEms.m       | Selects data matrix elements for some m                        |
| LegendreCoeffs.m  | Coefficients of Legendre polynomials                           |
| PolypartCoeffs.m  | Coefficients  of the polynomial part of the spherical harmonic |
| CramerRaoBound.m  | Computation of the Cramer Rao lower bound                      |
| CmMatrix.m        | Matrix of polynomial coefficients corresponding to m           |

__Note__: The code to generate the spherical harmonics currently can't generate harmonics of high orders (to be replaced soon).


## Examples and paper figures

| Filename             | Description                                                       |
| -----------------    | -------------------                                               |
| Example1_Diracs.m    | Simple reconstruction of a collection of Diracs                   |
| Example2_Green.m     | Some plots of the Helmholtz Green's function ove r a rigid sphere |
| Example4_Green.m     | Figure 7 from the paper (source localization)                     |
| Example5_Green.m     | Figure 8 from the paper (ratios of Green's functions)             |
| Example6_Diffusion.m | Diffusion kernel and aliasing error                               |
| Example7_Diffusion.m | Diffusion reconstruction example                                  |
| Example8_ShotNoise.m | Shot noise removal example                                        |

## Authors

Ivan Dokmanić is with Laboratory for Audiovisual Communications
([LCAV](http://lcav.epfl.ch)) at [EPFL](http://www.epfl.ch), Switzerland. Yue M. Lu is with Signals, Information, and Networks Group ([SING][http://lu.seas.harvard.edu]) at Harvard University, Cambridge, MA, USA.

<img src="http://lcav.epfl.ch/files/content/sites/lcav/files/images/Home/LCAV_anim_200.gif">
