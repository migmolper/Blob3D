# Blob3D

**Blob3D** is a C++ implementation of the *blob* (finite-particle) method for mass diffusion based on optimal transport. Mass is carried by finite-width particles whose motion follows from a time-discrete incremental variational principle that balances entropy and the Wasserstein cost of transport.

The formulation and algorithms implemented here follow:

- [Pandolfi, Stainier & Ortiz (2023)](https://arxiv.org/abs/2305.05315) — *An optimal-transport finite-particle method for mass diffusion* ([arXiv:2305.05315](https://arxiv.org/abs/2305.05315))
- [Pandolfi, Romero & Ortiz (2025)](https://doi.org/10.1016/j.cma.2025.118013) — *An optimal-transport finite-particle method for driven mass diffusion* ([CaltechAUTHORS](https://authors.library.caltech.edu/records/jyr2m-rtx05))

The method is meshless in any spatial dimension, can redistribute particles and track their evolution in time, and enforces general flux / mixed boundary conditions (including via an adsorption/depletion layer at the boundary). Built on [PETSc](https://petsc.org), Blob3D is parallelized with MPI and targets modern HPC systems.

## Documentation

For building **Blob3D as an installable library** and linking it from another C++ project, see [LIBRARY_BUILD.md](LIBRARY_BUILD.md).

Machine-specific configure scripts live under [`configs/`](configs/). An example application is available in [`examples/sphere/`](examples/sphere/).

## Contributors

* [Miguel Molinos](https://github.com/migmolper) (Asst. Prof.). Universidad Politécnica de Madrid, Spain. Main developer.

See also [AUTHORS.md](AUTHORS.md).

## References

- A. Pandolfi, L. Stainier, M. Ortiz. *An optimal-transport finite-particle method for mass diffusion*. [arXiv:2305.05315](https://arxiv.org/abs/2305.05315), 2023.
- A. Pandolfi, I. Romero, M. Ortiz. *An optimal-transport finite-particle method for driven mass diffusion*. *Computer Methods in Applied Mechanics and Engineering*, 442:118013, 2025. [doi:10.1016/j.cma.2025.118013](https://doi.org/10.1016/j.cma.2025.118013)

Related open-source finite-particle code by the authors of the 2025 paper: [gitlab.com/ignacio.romero/finite-particles](https://gitlab.com/ignacio.romero/finite-particles).

## Useful links and software

* [Ovito](https://www.ovito.org) (postprocessor)
* [ParaView](https://www.paraview.org) (postprocessor)
* [PETSc](https://petsc.org)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](license.md)
