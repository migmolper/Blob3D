/**
 * @file Simulation.hpp
 * @author Miguel Molinos ([migmolper](https://github.com/migmolper))
 * @brief RAII simulation facade for Blob3D particle systems.
 * @version 0.1
 * @date 2026-07-16
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef BLOB3D_SIMULATION_HPP
#define BLOB3D_SIMULATION_HPP

#include "Macros.hpp"
#include <utility>

/*******************************************************/

/**
 * @brief Thermodynamic environment driving the simulation
 *
 */
struct Environment {

  //! @param pressure: Environmental value of the pressure
  double pressure{0.0};

  //! @param temperature: Environmental value of the temperature
  double temperature{0.0};

  //! @param chemical_potential: Environmental chemical potential ({mu}).
  //! {mu} = {gamma}/{beta}
  double chemical_potential{0.0};
};

/*******************************************************/

/**
 * @brief Owns the PETSc DMSwarm (and its background mesh) plus site counts.
 *
 * Move-only. The background FE mesh attached via DMSwarmSetCellDM is destroyed
 * together with the swarm.
 */
class ParticleSwarm {
 public:
  ParticleSwarm() = default;
  ~ParticleSwarm();

  ParticleSwarm(ParticleSwarm&& other) noexcept;
  ParticleSwarm& operator=(ParticleSwarm&& other) noexcept;
  ParticleSwarm(const ParticleSwarm&) = delete;
  ParticleSwarm& operator=(const ParticleSwarm&) = delete;

  /**
   * @brief Access the owned DMSwarm
   *
   * @return DM Particle swarm (PIC)
   */
  DM dm() const { return dm_; }

  /**
   * @brief Access the dump-file to PETSc index mapping
   *
   * @return AO Application ordering between dump and PETSc numbering
   */
  AO dump_to_petsc() const { return dump2petsc_; }

  /**
   * @brief Number of particles in the global domain
   */
  PetscInt n_global() const { return n_global_; }
  PetscInt& n_global() { return n_global_; }

  /**
   * @brief Number of particles in the local domain (without ghost)
   */
  PetscInt n_local() const { return n_local_; }
  PetscInt& n_local() { return n_local_; }

  /**
   * @brief Number of ghost particles in the local domain
   */
  PetscInt n_ghost() const { return n_ghost_; }
  PetscInt& n_ghost() { return n_ghost_; }

  /**
   * @brief Destroy owned PETSc objects and reset counts
   */
  void release() noexcept;

  /**
   * @brief Take ownership of PETSc handles created during Simulation::initialize
   *
   * @param dm DMSwarm object
   * @param dump2petsc AO mapping from dump ordering to PETSc ordering
   * @param n_global Number of particles in the global domain
   * @param n_local Number of particles in the local domain (without ghost)
   */
  void adopt(DM dm, AO dump2petsc, PetscInt n_global, PetscInt n_local);

 private:
  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   @brief System information
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/

  //! @param dm_: Particle data (DMSwarm)
  DM dm_{nullptr};

  //! @param dump2petsc_: Track the index of each particle from .dump to PETSc
  //! ordering
  AO dump2petsc_{nullptr};

  //! @param n_global_: Number of particles in the global domain
  PetscInt n_global_{0};

  //! @param n_local_: Number of particles in the local domain (without ghost)
  PetscInt n_local_{0};

  //! @param n_ghost_: Number of ghost particles in the local domain
  PetscInt n_ghost_{0};
};

/*******************************************************/

/**
 * @brief Owns per-particle neighbor index sets (mechanical topology)
 *
 * Each entry is an IS listing the mechanical neighbors of a local (ghosted)
 * particle. Ownership of the IS array is transferred via adopt().
 */
class NeighborTopology {
 public:
  NeighborTopology() = default;
  ~NeighborTopology();

  NeighborTopology(NeighborTopology&& other) noexcept;
  NeighborTopology& operator=(NeighborTopology&& other) noexcept;
  NeighborTopology(const NeighborTopology&) = delete;
  NeighborTopology& operator=(const NeighborTopology&) = delete;

  /**
   * @brief Take ownership of an array of IS objects of length @p n
   *
   * @param neighs Array of IS (mechanical neighbors per site)
   * @param n Length of the array (typically local ghosted size)
   */
  void adopt(IS* neighs, PetscInt n);

  /**
   * @brief Destroy all IS objects and free the array
   *
   * @return PetscErrorCode
   */
  PetscErrorCode clear();

  /**
   * @brief Raw pointer to the IS array (may be nullptr)
   */
  IS* data() const { return neighs_; }
  IS*& data() { return neighs_; }

  /**
   * @brief Number of IS entries currently owned
   */
  PetscInt size() const { return n_; }

  /**
   * @brief True if no neighbor lists are stored
   */
  bool empty() const { return neighs_ == nullptr; }

 private:
  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   @brief Topological variables of each site
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/

  //! @param neighs_: Local indices of mechanical neighbors per site
  IS* neighs_{nullptr};

  //! @param n_: Length of neighs_
  PetscInt n_{0};
};

/*******************************************************/

/**
 * @brief Top-level Blob3D simulation: particles + topology + environment
 *
 * Move-only facade used by solvers and I/O. Construct an empty Simulation and
 * call initialize() (PETSc error-code style). Topology helpers mirror the
 * former DMSwarm*BlobsTopology API.
 */
class Simulation {
 public:
  Simulation() = default;
  ~Simulation();

  Simulation(Simulation&&) noexcept = default;
  Simulation& operator=(Simulation&&) noexcept = default;
  Simulation(const Simulation&) = delete;
  Simulation& operator=(const Simulation&) = delete;

  /**
   * @brief Create the background mesh, DMSwarm fields and dump→PETSc mapping
   *
   * @param dump Input dump_file with particle data
   * @param mesh_type Background mesh type (DMDA / DMPLEX)
   * @param r_cutoff Cutoff radius used to size the background mesh
   * @return PetscErrorCode
   */
  PetscErrorCode initialize(const dump_file& dump, BackgroundMeshType mesh_type,
                            double r_cutoff);

  /**
   * @brief Create ghost particles and mechanical neighbor lists
   *
   * @param buffer_width Search / ghost buffer width
   * @return PetscErrorCode
   */
  PetscErrorCode generate_topology(double buffer_width);

  /**
   * @brief Destroy topology, rebin particles and rebuild ghosts + neighbors
   *
   * @param buffer_width Search / ghost buffer width
   * @return PetscErrorCode
   */
  PetscErrorCode regenerate_topology(double buffer_width);

  /**
   * @brief Destroy neighbor lists and ghost particles
   *
   * @return PetscErrorCode
   */
  PetscErrorCode destroy_topology();

  /**
   * @brief Access the particle swarm sub-object
   */
  ParticleSwarm& particles() { return particles_; }
  const ParticleSwarm& particles() const { return particles_; }

  /**
   * @brief Access the neighbor topology sub-object
   */
  NeighborTopology& topology() { return topology_; }
  const NeighborTopology& topology() const { return topology_; }

  /**
   * @brief Access the thermodynamic environment
   */
  Environment& env() { return env_; }
  const Environment& env() const { return env_; }

  /*! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   @brief Convenience accessors used throughout the solvers
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !*/
  //! @brief Particle data (DMSwarm)
  DM dm() const { return particles_.dm(); }

  //! @brief Dump-file to PETSc index mapping
  AO dump_to_petsc() const { return particles_.dump_to_petsc(); }

  //! @param n_sites_global: Number of particles in the global domain
  PetscInt n_sites_global() const { return particles_.n_global(); }
  PetscInt& n_sites_global() { return particles_.n_global(); }

  //! @param n_sites_local: Number of particles in the local domain (without ghost)
  PetscInt n_sites_local() const { return particles_.n_local(); }
  PetscInt& n_sites_local() { return particles_.n_local(); }

  //! @param n_ghost: Number of ghost particles in the local domain
  PetscInt n_ghost() const { return particles_.n_ghost(); }
  PetscInt& n_ghost() { return particles_.n_ghost(); }

  //! @param mechanical_neighs_idx: Local indices of mechanical neighbors
  IS* mechanical_neighs_idx() const { return topology_.data(); }
  IS*& mechanical_neighs_idx() { return topology_.data(); }

 private:
  //! @param particles_: Particle swarm and site counts
  ParticleSwarm particles_;

  //! @param topology_: Mechanical neighbor lists
  NeighborTopology topology_;

  //! @param env_: Environmental thermo-chemo-mechanical variables
  Environment env_;
};

#endif /* BLOB3D_SIMULATION_HPP */
