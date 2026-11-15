#pragma once

#include <array>
#include <cerrno>
#include <cstring>
#include <cstdint>
#include <fcntl.h>
#include <functional>
#include <new>
#include <stdexcept>
#include <string>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>

#define SHM_NAME "/sim_shm_w0"

namespace aocs::ipc {

inline constexpr double kGpsSecondsAtJ2000 = 630763200.0;
inline constexpr std::size_t kStateDimension = 6;

}

struct sensorMessage {
    std::array<double, 3> GPS_position_ECEF_m{};
    std::array<double, 3> GPS_velocity_ECEF_m_s{};
    double GPS_time_s = 0.0;
    double aocs_time_J2000_s = 0.0;
    std::uint32_t has_new_GPS_measurement = 0;
    std::array<double, aocs::ipc::kStateDimension * aocs::ipc::kStateDimension>
        GPS_covariance_row_major{};
};

struct navigationEstimateMessage {
    std::uint32_t valid = 0;
    double time_J2000_s = 0.0;
    std::array<double, 3> position_ECI_m{};
    std::array<double, 3> velocity_ECI_m_s{};
    std::array<double, aocs::ipc::kStateDimension * aocs::ipc::kStateDimension>
        covariance_row_major{};
};

struct sharedMemoryMap {
    sensorMessage sensorData;
    navigationEstimateMessage navigationEstimate;
};

class SharedMemoryRegion {
public:
    explicit SharedMemoryRegion(const std::string& name = SHM_NAME) : name_(name) {
        descriptor_ = shm_open(name_.c_str(), O_CREAT | O_RDWR, 0600);
        if (descriptor_ == -1) {
            throw std::runtime_error(
                "Unable to create shared memory " + name_ + ": " +
                std::strerror(errno));
        }

        if (flock(descriptor_, LOCK_EX | LOCK_NB) == -1) {
            const int error = errno;
            close(descriptor_);
            descriptor_ = -1;
            throw std::runtime_error(
                "Unable to lock shared memory " + name_ + ": " +
                std::strerror(error));
        }
        lockHeld_ = true;

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
        if (lockHeld_ && !name_.empty()) {
            shm_unlink(name_.c_str());
            name_.clear();
        }
        if (descriptor_ != -1) {
            if (lockHeld_) {
                flock(descriptor_, LOCK_UN);
                lockHeld_ = false;
            }
            close(descriptor_);
            descriptor_ = -1;
        }
    }

    std::string name_;
    int descriptor_ = -1;
    sharedMemoryMap* data_ = nullptr;
    bool lockHeld_ = false;
};

class GncTransport {
public:
    virtual ~GncTransport() = default;
    virtual int exchange(const sensorMessage& input) = 0;
};

class SharedMemoryGncTransport final : public GncTransport {
public:
    using GncEntryPoint = std::function<int(sharedMemoryMap&)>;

    SharedMemoryGncTransport(sharedMemoryMap& memory, GncEntryPoint gncEntryPoint)
        : memory_(memory), gncEntryPoint_(std::move(gncEntryPoint)) {
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
