#pragma once

#include "ipc.hpp"
#include "navigation/orbit_determination.hpp"


namespace aocs {

int main(sharedMemoryMap& sharedMemory);

class Navigation {
public:
    Navigation() = default;
    ~Navigation() = default;

};
} // namespace aocs
