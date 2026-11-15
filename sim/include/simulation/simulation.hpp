#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace aocs::simulation {

struct SpacecraftState {
    Eigen::Vector3d position_eci_m;
    Eigen::Vector3d velocity_eci_m_s;
};

class Propagation {
public:
    using AccelerationFunction =
        std::function<Eigen::Vector3d(const SpacecraftState&, double)>;

    class RK4 {
    public:
        static void integrate(
            SpacecraftState& state,
            double time_s,
            double dt_sec,
            const AccelerationFunction& acceleration) {
            if (!std::isfinite(time_s) || !std::isfinite(dt_sec) ||
                dt_sec <= 0.0 || !acceleration) {
                throw std::invalid_argument(
                    "RK4 requires finite time, a positive finite timestep, and an acceleration function");
            }

            const StateDerivative k1 = evaluate(state, time_s, acceleration);
            const StateDerivative k2 = evaluate(
                offsetState(state, k1, 0.5 * dt_sec),
                time_s + 0.5 * dt_sec,
                acceleration);
            const StateDerivative k3 = evaluate(
                offsetState(state, k2, 0.5 * dt_sec),
                time_s + 0.5 * dt_sec,
                acceleration);
            const StateDerivative k4 = evaluate(
                offsetState(state, k3, dt_sec),
                time_s + dt_sec,
                acceleration);

            const double scale = dt_sec / 6.0;
            SpacecraftState nextState{
                state.position_eci_m +
                    scale * (k1.positionRate + 2.0 * k2.positionRate +
                             2.0 * k3.positionRate + k4.positionRate),
                state.velocity_eci_m_s +
                    scale * (k1.velocityRate + 2.0 * k2.velocityRate +
                             2.0 * k3.velocityRate + k4.velocityRate)};
            if (!nextState.position_eci_m.allFinite() ||
                !nextState.velocity_eci_m_s.allFinite()) {
                throw std::runtime_error("RK4 produced a non-finite spacecraft state");
            }
            state = nextState;
        }

    private:
        struct StateDerivative {
            Eigen::Vector3d positionRate;
            Eigen::Vector3d velocityRate;
        };

        static SpacecraftState offsetState(
            const SpacecraftState& state,
            const StateDerivative& derivative,
            double scale) {
            return {
                state.position_eci_m + scale * derivative.positionRate,
                state.velocity_eci_m_s + scale * derivative.velocityRate};
        }

        static StateDerivative evaluate(
            const SpacecraftState& state,
            double time_s,
            const AccelerationFunction& acceleration) {
            const Eigen::Vector3d velocityRate = acceleration(state, time_s);
            if (!velocityRate.allFinite()) {
                throw std::runtime_error("RK4 acceleration must be finite");
            }
            return {state.velocity_eci_m_s, velocityRate};
        }
    };

    void propagate(
        SpacecraftState& state,
        double& elapsedTime_s,
        double dt_sec,
        const AccelerationFunction& acceleration) const {
        RK4::integrate(state, elapsedTime_s, dt_sec, acceleration);
        elapsedTime_s += dt_sec;
    }
};

} // namespace aocs::simulation
