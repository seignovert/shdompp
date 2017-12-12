# [Plane-parallel SHDOM for Data Assimilation (_SHDOMPPDA_)](http://nit.colorado.edu/~evans/shdomppda/index.html)

## Overview of _SHDOMPP_ and _SHDOMPPDA_
This web page is the distribution point for two closely related radiative transfer models. _SHDOMPP_ is the plane-parallel version of SHDOM. _SHDOMPPDA_ is a version of _SHDOMPP_ for data assimilation. _SHDOMPPDA_ includes forward, tangent linear, and adjoint models of _SHDOMPP_ packaged in a Fortran 90 module with a easy to use interface.

_SHDOMPP_ computes unpolarized radiative transfer in a plane-parallel medium for either collimated solar and/or thermal emission sources of radiation. Currently Lambertian and RPV land surface BRDF surface reflection are supported. The atmosphere optical properties are specified with an input file listing the layer interface heights and temperatures, and layer optical depth, single scattering albedo, and Legendre series coefficients of the phase function. Two output files are generated: one with the upwelling and downwelling hemispheric fluxes and actinic fluxes at the layer boundaries and another with radiances at specified heights and directions. The _SHDOMPP_ calculation may be either monochromatic or spectrally integrated with a k-distribution.

_SHDOMPPDA_ inputs pressure, temperature, water vapor, and ozone profiles, and mixing ratio and number concentration profiles for any number of hydrometeor species. The pressure, temperature, water vapor and ozone profiles are passed to user provided molecular absorption routines (_i.e._ an interface for molecular absorption routines is provided, but actual algorithms are not). _SHDOMPPDA_ includes optical property calculations using scattering tables as a function of mass mean radius. Mass mean radius is chosen instead of effective radius because it can be directly calculated from the input hydrometeor mass mixing ratio and number concentration without further assumptions about the particle size distribution and shape. Molecular Rayleigh scattering is calculated for wavelengths shorter than 1 micron. The _SHDOMPPDA_ forward radiative transfer routine calculates the top of atmosphere monochromatic radiances (in reflectance units or brightness temperature) for a number of specified directions for a single input column of properties.

---

## Documentation
A journal article on the _SHDOMPP_ and _SHDOMPPDA_ algorithms has been published and is available as a [PDF file](http://nit.colorado.edu/~evans/shdomppda/Evans2007JAS_SHDOMPPDA.pdf) (870 kB, [doi:10.1175/2006JAS2047.1](https://doi.org/10.1175/2006JAS2047.1)). (Evans, K. F., 2007: _SHDOMPPDA_: A radiative transfer model for cloudy sky data assimilation. J. Atmos. Sci., 64, 3858-3868.) The [_SHDOMPP_ documentation](http://nit.colorado.edu/~evans/shdomppda/shdompp.doc) file has the input file format, list of input parameters, information about running _SHDOMPP_ and other programs in the distribution, and a list of files in the distribution.

The [_SHDOMPPDA_ documentation](http://nit.colorado.edu/~evans/shdomppda/shdomppda.txt) file has an overview of the *SHDOMPPDA_FWD_TL_ADJ* module, description of arguments for key routines, information about running the two programs to make scattering tables, and a list of files in the distribution.

---

## _SHDOMPP_ and _SHDOMPPDA_ Distributions
The models are distributed via this Web page. _SHDOMPP_ and the forward model of _SHDOMPPDA_ are freely distributed. The source files for the _SHDOMPPDA_ tangent linear (`shdomppda_ftls.f90`) and the adjoint (`shdomppda_ads.f90`) cannot be freely distributed due to licensing requirements by FASTOPT, the owner of the TAF compiler used to produce these files. Dummy routines are distributed instead. If the tangent linear or adjoint is needed, please request permission from Thomas.Kaminski@FastOpt.com, and then contact me ([evans@nit.colorado.edu](mailto:evans@nit.colorado.edu)) with evidence of this permission, so I can email the files to you.

_SHDOMPP_ and _SHDOMPPDA_ are distributed as a gzipped Unix tar files. Each distribution contains a documentation file, Fortran 90 source code and a makefile, and example scripts and data files that illustrate using the models. _SHDOMPPDA_ is meant to be used by calling the Fortran subroutines in the shdomppda_fwd_tl_adj module, but an example program is included.

---

## Original sources
- The [_SHDOMPP_ distribution](http://nit.colorado.edu/~evans/shdomppda/shdompp.tar.gz) is `275kB`.
- The [_SHDOMPPDA_ distribution](http://nit.colorado.edu/~evans/shdomppda/shdomppda.tar.gz) is `903 kB`.
- Get [`yang2005_ice_scatter.db`](http://nit.colorado.edu/~evans/shdomppda/yang2005_ice_scatter.db) (`16 MB`) for making _SHDOMPPDA_ ice crystal scattering tables.

---

[_SHDOMPPDA_](http://nit.colorado.edu/~evans/shdomppda/index.html) was developped by [Frank Evans](http://nit.colorado.edu/~evans/) from the [Dept. of Atmospheric and Oceanic Sciences](http://atoc.colorado.edu/) at the [University of Colorado, Boulder](http://www.colorado.edu/)
