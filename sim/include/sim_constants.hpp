#pragma once

class SimConstants {
    /*
    @brief Constants for use in simulation.
    Convention:
    All scale factors assume multiplicative use as:
    out = kConstant * in;
    e.g.
    angleRadians = kDegreesToRadians * angleDegrees;
    All offsets assume additive use as:
    out = in + kConstant;
    e.g.
    julianDate = secondsSinceJ2000 / 86400.0 + kJulianDateJ2000toJulianDate;
    */
public:
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr double kJulianDateJ2000toJulianDate = 2451545.0;
    static constexpr double kSecondsPerDay = 86400.0;
    static constexpr double kRadtoDeg = 180.0 / kPi;
    static constexpr double kDegtoRad = kPi / 180.0;

    // Earth Constants
    static constexpr double kEarthRadius = 6371000.0; // meters
    static constexpr double kEarthMass = 5.972e24; // kilograms
    static constexpr double kGravitationalConstant = 6.67430e-11; // m^3 kg^-1 s^-2
    static constexpr double kEarthMu = kGravitationalConstant * kEarthMass; // m^3 s^-2
};
