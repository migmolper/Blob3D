/**
 * @file Blobs/IO/dump-inputs.hpp
 * @author M.Molinos (@migmolper)
 * @brief
 * @version 0.1
 * @date 2023-12-12
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef BLOBS_DUMP_IO_HPP
#define BLOBS_DUMP_IO_HPP

#include "Macros.hpp"

dump_file DMSwarmBlobsReadDump(const char* SimulationFile);

void DMSwarmBlobsFreeDump(dump_file* Simulation_data);

#endif /* BLOBS_DUMP_IO_HPP */