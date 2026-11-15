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
`aocs::runGncStep()` each simulation step. The GNC entry point converts the
ECEF measurement to ECI, initializes its estimate from the first measurement,
and publishes the ECI state and 6x6 covariance in
`sharedMemoryMap::navigationEstimate`. Its `valid` field indicates whether
the published estimate is usable. The simulation currently supplies an
identity GPS measurement covariance as a placeholder; replace it with the
sensor model's covariance when available.

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
propagates and runs the AOCS loop at 10 Hz, while GPS measurements are sampled
at 1 Hz and held between updates. The output CSV includes the state estimate
and full 6x6 estimate covariance at every AOCS step.

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
