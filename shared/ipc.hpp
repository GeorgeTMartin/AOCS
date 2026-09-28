// Interprocess Communication - shared memory struct
#pragma once

#include <Eigen/Dense>

// Standard C++ and POSIX headers for interprocess communication
#include <cstdint>
#include <fcntl.h>
#include <semaphore.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define SHM_NAME "/sim_shm_w0"               // location of share memory
#define SEM_SIM_TO_FSW "/sem_sim2fsw_w0"     // semaphore for simulation to flight software communication
#define SEM_FSW_TO_SIM "/sem_fsw2sim_w0"     // semaphore for flight software to simulation communication

struct sensorMessage {
    Eigen::Matrix<double, 3, 1> GPS_position_ECI;   // TODO: replace with ECEF
    Eigen::Matrix<double, 3, 1> GPS_velocity_ECI;   // TODO: replace with ECEF
    Eigen::Matrix<double, 3, 1> GPS_time;           // seconds elapsed since the GPS epoch
};

struct sharedMemoryMap {
    sensorMessage sensorData; // structure to hold sensor data
};