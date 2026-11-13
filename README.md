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
