/**
 * @file Blobs/Topology.hpp
 * @brief Topology helpers (thin wrappers around Simulation methods).
 */

#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "Simulation.hpp"

inline PetscErrorCode DMSwarmGenerateBlobsTopology(Simulation& simulation,
                                                   double r_cutoff) {
  return simulation.generate_topology(r_cutoff);
}

inline PetscErrorCode DMSwarmRegenerateBlobsTopology(Simulation& simulation,
                                                     double r_cutoff) {
  return simulation.regenerate_topology(r_cutoff);
}

inline PetscErrorCode DMSwarmDestroyBlobsTopology(Simulation& simulation) {
  return simulation.destroy_topology();
}

#endif /* TOPOLOGY_HPP */
