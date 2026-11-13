#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

class EGM2008Coefficients {
public:
    struct Coefficient {
        std::uint16_t degree;
        std::uint16_t order;
        double C;
        double S;
        double sigmaC;
        double sigmaS;
    };

    explicit EGM2008Coefficients(const std::string& path) {
        std::ifstream input(path, std::ios::binary);
        if (!input) {
            throw std::runtime_error("Unable to open EGM2008 binary: " + path);
        }

        char magic[8];
        input.read(magic, sizeof(magic));
        const char expectedMagic[8] = {'E', 'G', 'M', '2', '0', '0', '8', '\0'};
        if (!input || std::memcmp(magic, expectedMagic, sizeof(magic)) != 0) {
            throw std::runtime_error("Invalid EGM2008 binary magic");
        }
        if (readU32(input) != 1) {
            throw std::runtime_error("Unsupported EGM2008 binary version");
        }

        earthGravityConstant = readDouble(input);
        referenceRadius = readDouble(input);
        maxDegree = readU32(input);
        const std::uint64_t recordCount = readU64(input);
        if (maxDegree < 2 || recordCount != expectedRecordCount(maxDegree)) {
            throw std::runtime_error("Invalid EGM2008 binary dimensions");
        }

        productType = readString(input);
        modelName = readString(input);
        errors = readString(input);
        norm = readString(input);
        tideSystem = readString(input);

        coefficients.reserve(static_cast<std::size_t>(recordCount));
        std::uint32_t degree = 0;
        std::uint32_t order = 0;
        for (std::uint64_t index = 0; index < recordCount; ++index) {
            while (degree == 1 || order > degree) {
                ++degree;
                order = 0;
            }

            Coefficient coefficient{
                readU16(input),
                readU16(input),
                readDouble(input),
                readDouble(input),
                readDouble(input),
                readDouble(input)
            };
            if (coefficient.degree != degree || coefficient.order != order) {
                throw std::runtime_error("Invalid EGM2008 coefficient sequence");
            }
            coefficients.push_back(coefficient);
            ++order;
        }
        if (input.peek() != std::char_traits<char>::eof()) {
            throw std::runtime_error("Unexpected trailing data in EGM2008 binary");
        }
    }

    const Coefficient& coefficient(std::uint16_t degree, std::uint16_t order) const {
        if (degree > maxDegree || degree == 1 || order > degree) {
            throw std::out_of_range("EGM2008 coefficient index is out of range");
        }
        const std::uint64_t index =
            degree == 0 ? 0 : static_cast<std::uint64_t>(degree) * (degree + 1) / 2 - 2 + order;
        return coefficients.at(static_cast<std::size_t>(index));
    }

    double earthGravityConstant;
    double referenceRadius;
    std::uint32_t maxDegree;
    std::string productType;
    std::string modelName;
    std::string errors;
    std::string norm;
    std::string tideSystem;

private:
    static std::uint64_t expectedRecordCount(std::uint32_t degree) {
        return static_cast<std::uint64_t>(degree + 1) * (degree + 2) / 2 - 2;
    }

    static std::uint16_t readU16(std::istream& input) {
        unsigned char bytes[2];
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        if (!input) {
            throw std::runtime_error("Truncated EGM2008 binary");
        }
        return static_cast<std::uint16_t>(bytes[0])
            | (static_cast<std::uint16_t>(bytes[1]) << 8);
    }

    static std::uint32_t readU32(std::istream& input) {
        unsigned char bytes[4];
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        if (!input) {
            throw std::runtime_error("Truncated EGM2008 binary");
        }
        return static_cast<std::uint32_t>(bytes[0])
            | (static_cast<std::uint32_t>(bytes[1]) << 8)
            | (static_cast<std::uint32_t>(bytes[2]) << 16)
            | (static_cast<std::uint32_t>(bytes[3]) << 24);
    }

    static std::uint64_t readU64(std::istream& input) {
        unsigned char bytes[8];
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        if (!input) {
            throw std::runtime_error("Truncated EGM2008 binary");
        }
        std::uint64_t value = 0;
        for (unsigned int index = 0; index < 8; ++index) {
            value |= static_cast<std::uint64_t>(bytes[index]) << (index * 8);
        }
        return value;
    }

    static double readDouble(std::istream& input) {
        const std::uint64_t bits = readU64(input);
        double value;
        std::memcpy(&value, &bits, sizeof(value));
        return value;
    }

    static std::string readString(std::istream& input) {
        const std::uint32_t length = readU32(input);
        if (length > 1024) {
            throw std::runtime_error("Invalid EGM2008 metadata length");
        }
        std::string value(length, '\0');
        input.read(value.data(), length);
        if (!input) {
            throw std::runtime_error("Truncated EGM2008 metadata");
        }
        return value;
    }

    std::vector<Coefficient> coefficients;
};
