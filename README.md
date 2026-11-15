# AOCS

This repository is intended to provide the basis of development for primitive AOCS capabilities, with an associated simulation framework.

## Build

Requirements: CMake 3.20 or newer, a C++17 compiler, the JsonCpp and pybind11 development packages, and Python development headers. Eigen is included in the repository.

```sh
sudo apt-get update
sudo apt-get install cmake g++ libjsoncpp-dev pybind11-dev python3-dev
```

From a fresh clone, configure and build with:

```sh
cmake -S . -B build
cmake --build build
```

Run the simulation from the repository root:

```sh
./build/sim/sim
```

## Simulation and GNC interface

The simulation publishes its GPS position, GPS velocity, and GPS time inputs
to a POSIX shared-memory region, then calls `aocs::gnc::main` through the
`GncTransport` interface. The current GNC entry point is intentionally a
no-op scaffold. The shared-memory transport is synchronous and in-process for
now; a future UDP transport can implement the same interface when GNC runs in
a separate container.

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

After running the simulation, generate the orbital-element time-series plot
from the repository root:

```sh
python plotting/plot_results.py
```

The simulation CSV stores only time, position, and velocity. The plotter calls
the `aocs_orbit` pybind11 module, which invokes the C++
`IJKKeplerian::convert_ijk_to_keplerian` function for each state vector. It
plots the osculating semi-major axis, eccentricity, inclination, RAAN,
argument of periapsis, and true anomaly in six panels in `output.html`.
Angular values are converted from radians to degrees for plotting. The report
also includes a dark, retro-styled 3D ECI orbit view with an Earth sphere at
the origin.
