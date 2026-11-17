#pragma once

namespace flight_constants {
    constexpr double kPi = 3.14159265358979323846;

    constexpr double kEarthMu = 3.986004418e14;
    constexpr double kEarthRadius = 6.371e6;
    constexpr double kEarthAngularRate = 7.2921150e-5;

    constexpr double kGPSSecondsAtJ2000 = 630763200.0;
    constexpr double kUT1toTTOffset_s = 69.184;
    constexpr double kTTtoUT1Offset_s = -kUT1toTTOffset_s;
    constexpr double kJ2000toJulianDate = 2451545.0;

    constexpr double kArcsecondsToRadians = kPi / 648000.0;
    constexpr double kDegreesToRadians = kPi / 180.0;
}