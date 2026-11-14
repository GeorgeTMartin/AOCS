#include "astro/ijk_keplerian.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <stdexcept>

namespace py = pybind11;

py::array_t<double> convertIJKToKeplerianSeries(
    py::array_t<double, py::array::c_style | py::array::forcecast> positions,
    py::array_t<double, py::array::c_style | py::array::forcecast> velocities) {
    if (positions.ndim() != 2 || positions.shape(1) != 3 ||
        velocities.ndim() != 2 || velocities.shape(1) != 3) {
        throw std::invalid_argument(
            "positions and velocities must both have shape (N, 3)");
    }
    if (positions.shape(0) != velocities.shape(0)) {
        throw std::invalid_argument(
            "positions and velocities must contain the same number of states");
    }

    const auto positionValues = positions.unchecked<2>();
    const auto velocityValues = velocities.unchecked<2>();
    py::array_t<double> elements({positions.shape(0), py::ssize_t{6}});
    auto elementValues = elements.mutable_unchecked<2>();

    for (py::ssize_t row = 0; row < positions.shape(0); ++row) {
        const Eigen::Vector3d position(
            positionValues(row, 0), positionValues(row, 1), positionValues(row, 2));
        const Eigen::Vector3d velocity(
            velocityValues(row, 0), velocityValues(row, 1), velocityValues(row, 2));
        if (!position.allFinite() || !velocity.allFinite() ||
            !(position.norm() > 0.0) ||
            !(position.cross(velocity).norm() > 0.0)) {
            throw std::invalid_argument(
                "each state must have finite values, nonzero position, and nonzero angular momentum");
        }

        const orbitalElements state =
            IJKKeplerian::convert_ijk_to_keplerian(position, velocity);
        elementValues(row, 0) = state.sma;
        elementValues(row, 1) = state.e;
        elementValues(row, 2) = state.i;
        elementValues(row, 3) = state.omega;
        elementValues(row, 4) = state.argp;
        elementValues(row, 5) = state.nu;
    }

    return elements;
}

PYBIND11_MODULE(aocs_orbit, module) {
    module.doc() = "Bindings for AOCS orbital mechanics";
    module.def(
        "convert_ijk_to_keplerian_series",
        &convertIJKToKeplerianSeries,
        py::arg("positions_eci_m"),
        py::arg("velocities_eci_m_s"),
        "Convert matching arrays of ECI position and velocity vectors. "
        "Returns columns [sma_m, eccentricity, inclination_rad, raan_rad, "
        "argp_rad, true_anomaly_rad].");
}
