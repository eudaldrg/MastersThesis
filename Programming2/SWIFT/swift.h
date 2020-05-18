#pragma once

#include <functional>
#include <numeric>

#include "my_math.h"

namespace Swift {

inline double Rect(double x)
{
    double abs_x = std::abs(x);
    if (abs_x < 0.5)
        return 1;
    if (abs_x == 0.5) // TODO: A lot of double precision issues, but who cares, a single point doesn't affect the integral.
        return 0.5;
    //        if (abs_x > 0.5)
    return 0;
}

inline double Sinc(double x)
{
    return std::sin(PI * x) / (PI * x);
}

/// 2 / pi * (sum_{j=1}^{2^{J-1}} sin(pi * x * (2j - 1) / 2^J) / (2j - 1)
inline double SIApprox(double x, std::size_t J_bar)
{
    const uint64_t exp_2_J_minus_1 = 1UL << (J_bar - 1);
    const uint64_t exp_2_J = exp_2_J_minus_1 * 2;
    double result = 0.0;
    for (std::size_t j = 1; j <= exp_2_J_minus_1; ++j)
    {
        std::size_t two_j_minus_1 = 2 * j - 1;
        result += std::sin(two_j_minus_1 * Swift::PI * x / exp_2_J) / two_j_minus_1;
    }
    return result * 2 / Swift::PI;
}

/// Expression 14 in the paper. It computes the area under f through the trapezoid's rule
inline double AreaUnderCurve(std::size_t m, std::vector<double> c_m_k)
{
    return (std::accumulate(c_m_k.begin(), c_m_k.end(), 0.0) - c_m_k.front() * 0.5 - c_m_k.back() * 0.5)
           / std::pow(2, m * 0.5);
}

inline std::function<double(double)> GetTheta(std::size_t m, int k)
{
    return [m, k](double x)
    {
        return std::pow(2, m / 2) * Sinc(std::pow(2, m) * x - k);
    };
}

inline std::function<double(double)> GetPsi(std::size_t m, int k)
{
    return [m, k](double x)
    {
        return std::pow(2, m / 2) * (std::sin(PI * (std::pow(2, m) * x - k - 0.5)
                                              - std::sin(2 * PI * (std::pow(2, m) * x - k - 0.5))) / (PI * (std::pow(2, m) * x - k - 0.5)));
    };
}

inline double GetC_m_k(std::size_t m, int k, std::complex<double> IFT, std::size_t N)
{
    std::complex<double> complex_part = std::exp(1i * static_cast<double>(k) * PI / static_cast<double>(N)) * IFT;
    return std::pow(2.0, m * 0.5) * 2 * complex_part.real();
}

/// c_{m,k} \approx c^*_{m,k} = \dfrac{2^{m/2}}{2^{J-1} \sum_{j=1}^{2^{J-1}}} \R [\hat{f}(\dfrac{(2j - 1)
/// \pi 2^m}{2^J}) e^{\dfrac{ik\pi(2j-1)}{2^J}}])]}

//    double GetC_m_k_Error(std::size_t m, int k, double a, double e)
}