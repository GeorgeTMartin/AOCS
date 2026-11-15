#pragma once

#include <Eigen/Dense>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <new>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define SHM_NAME "/sim_shm_w0"

struct sensorMessage {
    Eigen::Vector3d GPS_position_ECEF = Eigen::Vector3d::Zero();
    Eigen::Vector3d GPS_velocity_ECEF = Eigen::Vector3d::Zero();
    double GPS_time_s = 0.0;
};

struct sharedMemoryMap {
    sensorMessage sensorData;
};

class SharedMemoryRegion {
public:
    explicit SharedMemoryRegion(const std::string& name = SHM_NAME) : name_(name) {
        descriptor_ = shm_open(name_.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
        if (descriptor_ == -1) {
            throw std::runtime_error(
                "Unable to create shared memory " + name_ + ": " +
                std::strerror(errno));
        }

        if (ftruncate(descriptor_, sizeof(sharedMemoryMap)) == -1) {
            const int error = errno;
            cleanup();
            throw std::runtime_error(
                "Unable to size shared memory " + name_ + ": " +
                std::strerror(error));
        }

        void* const mapping =
            mmap(nullptr, sizeof(sharedMemoryMap), PROT_READ | PROT_WRITE,
                 MAP_SHARED, descriptor_, 0);
        if (mapping == MAP_FAILED) {
            const int error = errno;
            cleanup();
            throw std::runtime_error(
                "Unable to map shared memory " + name_ + ": " +
                std::strerror(error));
        }

        try {
            data_ = ::new (mapping) sharedMemoryMap{};
        } catch (...) {
            munmap(mapping, sizeof(sharedMemoryMap));
            cleanup();
            throw;
        }
    }

    ~SharedMemoryRegion() {
        if (data_) {
            data_->~sharedMemoryMap();
            munmap(data_, sizeof(sharedMemoryMap));
        }
        cleanup();
    }

    SharedMemoryRegion(const SharedMemoryRegion&) = delete;
    SharedMemoryRegion& operator=(const SharedMemoryRegion&) = delete;
    SharedMemoryRegion(SharedMemoryRegion&&) = delete;
    SharedMemoryRegion& operator=(SharedMemoryRegion&&) = delete;

    sharedMemoryMap& data() {
        return *data_;
    }

private:
    void cleanup() {
        if (descriptor_ != -1) {
            close(descriptor_);
            descriptor_ = -1;
        }
        if (!name_.empty()) {
            shm_unlink(name_.c_str());
            name_.clear();
        }
    }

    std::string name_;
    int descriptor_ = -1;
    sharedMemoryMap* data_ = nullptr;
};

class GncTransport {
public:
    virtual ~GncTransport() = default;
    virtual int exchange(const sensorMessage& input) = 0;
};

class SharedMemoryGncTransport final : public GncTransport {
public:
    using GncEntryPoint = int (*)(sharedMemoryMap&);

    SharedMemoryGncTransport(sharedMemoryMap& memory, GncEntryPoint gncEntryPoint)
        : memory_(memory), gncEntryPoint_(gncEntryPoint) {
        if (!gncEntryPoint_) {
            throw std::invalid_argument("GNC entry point must not be null");
        }
    }

    int exchange(const sensorMessage& input) override {
        memory_.sensorData = input;
        return gncEntryPoint_(memory_);
    }

private:
    sharedMemoryMap& memory_;
    GncEntryPoint gncEntryPoint_;
};
