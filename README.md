# AOCS

This repository is intended to provide the basis of development for primitive AOCS capabilities, with an associated simulation framework.

## Git LFS

The large data binaries are stored with Git LFS, as they are relatively static and not expected to change as the repository evolves. Install and initialize Git
LFS before cloning:

```sh
sudo apt-get install git-lfs
git lfs install
git clone <repository-url>
```

For an existing clone, run `git lfs install` and then `git lfs pull` from the
repository root to download the large data files.

## Build

Requirements: CMake 3.20 or newer, a C++17 compiler, the JsonCpp and pybind11 development packages, and Python development headers. Eigen is included in the repository.

```sh
sudo apt-get update
sudo apt-get install cmake g++ libjsoncpp-dev pybind11-dev python3-dev ninja-build
```

From a fresh clone, configure and build with:

```sh
cmake -S . -B build -G Ninja
cmake --build build
```

Run the simulation from the repository root:

```sh
./build/sim/sim
```

## Simulation and GNC interface

The simulation initializes persistent GNC state once with
`aocs::initializeGnc()`, then publishes GPS position, velocity, time, and
measurement covariance to POSIX shared memory before calling
`aocs::runGncStep()` each simulation step. The GNC boundary maps shared-memory
messages into internal navigation inputs, evaluates the GNC state machine,
and calls the navigation step. Navigation owns the orbit-determination
instance, initializes from the first GPS measurement, and propagates the
estimate between measurements. The boundary maps the resulting ECI state and
6x6 covariance back into `sharedMemoryMap::navigationEstimate`; its `valid`
field indicates whether the estimate is usable. Guidance and control steps
are not implemented yet. The simulation currently supplies an identity GPS
measurement covariance as a placeholder; replace it with the sensor model's
covariance when available.
The GNC acceleration process-noise spectral density is configured internally
in [`gnc_config.hpp`](./gnc_lib/gnc_config.hpp), in radial, transverse, and
normal (RTN) order. The current `[1e-4, 1e-2, 1e-4] m²/s³` values are
provisional tuning values, not a sensor-certified specification.

The shared-memory transport is synchronous and in-process for now; a future
UDP transport can implement the same interface when GNC runs in a separate
container.

## Analysis

Requirements: Python 3.14 or newer

Create a virtual environment
```sh
# if necessary, install venv
apt install python3.14-venv
python3 -m venv .venv
source ./venv/bin/activate
```
Install requirements
```sh
pip install -r requirements.txt
```

The default simulation runs for two initial-condition orbital periods. It
propagates the true state at 100 Hz and runs the AOCS loop at 10 Hz, while GPS
measurements are sampled at 1 Hz and held between updates. The output CSV is
recorded at 10 Hz and includes the state estimate and full 6x6 estimate
covariance at each log step. Configure `true_state_rate_hz`, `aocs_rate_hz`,
`log_rate_hz`, and `gps_rate_hz` under `simulation` in `config.json`. The
true-state rate must be an integer multiple of the AOCS and log rates, and the
AOCS rate must be an integer multiple of the GPS rate.

After running the simulation, generate the interactive report from the
repository root:

```sh
python plotting/plot_results.py
```

The first HTML tab contains the 3D ECI orbit, position and velocity time
series in three columns, and the six osculating-element plots. The covariance
tab plots position and velocity variances and animates the 3σ position
covariance ellipsoid in a local offset frame. The plotter calls the
`aocs_orbit` pybind11 module, which invokes the C++
`IJKKeplerian::convert_ijk_to_keplerian` function for each state vector.
Angular elements are converted from radians to degrees.


## Roadmap
```text
Week 1: Inertial Nav & EKF MVP
├── High-Fidelity GPS
│   ├── GPS Sensor Model
│   └── Sensor Invalid Feature (General)
└── EKF MVP Architecture
    ├── Process Noise
    └── Increased STM fidelity - exponential drag, J2, Vallado Sun/Moon Ephem

Week 2: Hi-Fi Inertial Physics
├── Drag Model (MSIS) - Pull in src FORTRAN
├── Solar Pressure Model
└── DE440 Ephemeris

Week 3: Body Dynamics
├── Frame Library
│   ├── Sensor Frame
│   └── Body Frame
└── Bus Mass Simulation
    ├── CG Calculation - Wet/Dry Mass
    └── Inertia Dynamics - Tanks

Week 4: Attitude System
├── Attitude Dynamics
│   ├── Quaternion Attitude Representation
│   └── Attitude Propagation
└── Attitude Sensing
    ├── Sun Sensor
    ├── Magnetometer / IGRF or WMM
    ├── IMU
    └── Star Tracker

Week 5: Attitude Control I
├── MEKF - Markley
└── Attitude Control
    ├── Torque Rod
    └── B-Dot Control

Week 6: Hi-Fi Attitude Sim
├── Plate Model (SRP/Drag)
└── Gravity Gradient

Week 7: ADCS
├── Attitude Actuators
│   ├── Reaction Wheels
│   └── Thrusters
└── Attitude Control II
    ├── LQR Body Control
    ├── PD Momentum Control
    ├── L2 RW Allocation
    └── PID Torque Control

Week 8: Station Keeping Guidance
├── Hohmann Transfer
├── Bi-Elliptic Transfer
└── Inclination Change

Week 9: Hi-Fi Orbital Maneuvering
├── Thruster Allocation - Facet Method (?)
└── Mass Flow Dynamics      
```