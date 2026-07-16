# SOLERA

<p align="center">
  <img src="Logo.png" alt="My Logo" width="300"/>
</p>

SOLERA (soh-leh-rah) is a Diffusive Molecular Dynamics code designed to bridge atomistic and continuum time scales by transforming fast atomic vibrations into effective diffusive motions governed by stochastic dynamics. It enables the simulation of long-term phenomena such as phase transformations, defect migration, and hydrogen diffusion in metals while preserving atomic resolution and energy fidelity. Built on top of the PETSc library, SOLERA is fully parallelized using MPI and optimized for modern HPC architectures. Its modular design allows for flexible coupling between atomic interactions, continuum fields, and external driving forces, making it a powerful research tool for exploring microstructural evolution and multiscale transport phenomena in complex materials systems. This code is used as the basis for subsequent DMD developments by [TEP972](https://investigacion.us.es/sisius/sis_depgrupos.php?ct=&cs=&seltext=TEP-972&selfield=CodPAI).

## What means SOLERA?
A method of producing [sherry wine](https://en.m.wiktionary.org/wiki/sherry) in which small amounts of younger wines stored in an upper tier of casks are systematically blended with the more mature wine in the casks below. 

## Documentation

Please refer to [our wiki](https://github.com/migmolper/DMD/wiki) for information on compiling, and running the code.

For building **SOLERA as an installable library** and **linking it from another C++ project**, see [LIBRARY_BUILD.md](LIBRARY_BUILD.md) (copy under [docs/LIBRARY_BUILD.md](docs/LIBRARY_BUILD.md)).

For a **Spack-based** workflow (dependencies + install on clusters or laptops), see [SPACK.md](SPACK.md) (copy under [docs/SPACK.md](docs/SPACK.md)).

## Contributors 
* [Miguel Molinos](https://github.com/migmolper) (Asst. Prof.). Universidad Politécnica de Madrid, Spain. Main developer.
* [Jose Manuel Recio-López](https://github.com/jrecio1) (PhD. Student). Universidad de Sevilla, Spain. Al-Cu system.

Advised by
* [Pilar Ariza](https://github.com/mpariza) (Prof.). Universidad de Sevilla, Spain. 
* [Michael Ortiz](https://github.com/mortizcaltech) (Prof.). California Institute of Technology, USA.

## Additional context
- [Atomistic long-term simulation of heat and mass transport](https://www.sciencedirect.com/science/article/pii/S002250961400194X)
- [Acceleration of diffusive molecular dynamics simulations through mean field approximation and subcycling time integration](https://www.sciencedirect.com/science/article/pii/S0021999117306502)
- [Atomistic Simulation of Hydrogen Diffusion in Palladium Nanoparticles Using a Diffusive Molecular Dynamics Method](https://asmedigitalcollection.asme.org/IMECE/proceedings-abstract/IMECE2017/V009T12A026/261885)
- [Long-term atomistic simulation of hydrogen diffusion in metals](https://www.sciencedirect.com/science/article/pii/S0360319915001780)

## Useful links and software
* [Crystallography Open Database](http://www.crystallography.net/cod/search.html)
* [NIST Interatomic potentials](https://www.ctcms.nist.gov/potentials/)
* [Atomsk (preprocessor)](https://atomsk.univ-lille.fr)
* [Ovito (preprocessor/postprocessor)](https://www.ovito.org)
* [Paraview (postprocessor)](https://www.paraview.org)


[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/migmolper/DMD/blob/main/license.md)
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)](https://github.com/migmolper/DMD/wiki)
<!---
[![CircleCI](https://circleci.com/gh/geomechanics/mpm.svg?style=svg)](https://circleci.com/gh/geomechanics/mpm)
[![codecov](https://codecov.io/gh/geomechanics/mpm/branch/master/graph/badge.svg)](https://codecov.io/gh/geomechanics/mpm)
[![](https://img.shields.io/github/issues-raw/geomechanics/mpm.svg)](https://github.com/geomechanics/mpm/issues)
[![Coverity](https://scan.coverity.com/projects/14389/badge.svg)](https://scan.coverity.com/projects/14389/badge.svg)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/cb-geo/mpm.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/cb-geo/mpm/context:cpp)
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/geomechanics/mpm/projects/1)
-->

