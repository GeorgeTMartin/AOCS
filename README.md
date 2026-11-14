# AOCS

This repository is intended to provide the basis of development for primitive AOCS capabilities, with an associated simulation framework.

## Build

Requirements: CMake 3.20 or newer, a C++17 compiler, and the JsonCpp development package. Eigen is included in the repository.

```sh
sudo apt-get update
sudo apt-get install cmake g++ libjsoncpp-dev
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

## Earth rotation and quaternions

`sim/include/environment/earth_rotation/iau2000b.hpp` reduces GCRS (ECI)
coordinates to ITRS (ECEF) using IAU 2000B nutation, IAU 2000 precession
and frame bias, Greenwich apparent sidereal time, and polar motion. The
returned matrix and quaternion are passive coordinate transforms; applying
either maps a GCRS vector into ITRS. The quaternion is scalar-first and uses
the Hamilton product convention.

Supply TT and UT1 as two-part Julian dates. TT drives precession and
nutation; UT1 drives Earth rotation. Polar motion `xp` and `yp` are in
radians and default to zero. The library does not convert UTC to TT/UT1:
callers must provide the appropriate time scales and Earth orientation
parameters. For elapsed seconds from J2000.0, use
`JulianDate::fromJ2000Seconds` for each time scale.

The simulation reads `simulation.initial_conditions.t_J2000_s` from
`sim/config.json` as its initial epoch and adds the elapsed simulation time
before each gravity evaluation. Since the sample configuration provides one
epoch only, it uses that epoch for both TT and UT1 and assumes zero polar
motion.

```cpp
using aocs::environment::EarthRotation;
using aocs::environment::JulianDate;

const auto reduction =
    EarthRotation(tt, ut1, xpRadians, ypRadians).reduce();
const Eigen::Vector3d rItrs = reduction.gcrsToItrs * rGcrs;
const Eigen::Vector3d rItrsFromQuaternion =
    reduction.gcrsToItrsQuaternion.rotatePassive(rGcrs);
```

The IAU 2000B model is approximately 1 milliarcsecond, and this reduction
does not include observed celestial-pole offsets (`dX`, `dY`) or sub-daily
Earth-orientation corrections. The equation of the equinoxes includes its
IAU 2000 complementary terms. The reusable quaternion type is in
`shared/math/quaternion.hpp`.


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
