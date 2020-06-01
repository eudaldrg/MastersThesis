#pragma once

#include <functional>
#include <numeric>

#include "my_math.h"

class Distribution;
class OptionContract;

namespace Swift {

struct SwiftParameters
{
    SwiftParameters(size_t m, int k_1, int k_2, size_t iota, size_t iota_bar);

    std::size_t m_m;
    int m_k1;
    int m_k2;
    std::size_t m_iota_dens_coefs; // Iota for the C_m_k calculation.
    std::size_t m_iota_payoff_coefs; // Iota for the V_m_k calculation. (determines the truncation of the integration space).
};

class SwiftEvaluator
{
public:
    SwiftEvaluator(SwiftParameters const& params, Distribution const& distribution, OptionContract const& option_contract);

    [[nodiscard]] double GetPrice(double F, double K, double r, double q) const;

    SwiftParameters m_params;
    Distribution const& m_distribution;
    OptionContract const& m_option_contract;
};

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
        return std::pow(2, m / 2) * (std::sin(MY_PI * (std::pow(2, m) * x - k - 0.5)
                                              - std::sin(2 * MY_PI * (std::pow(2, m) * x - k - 0.5))) / (MY_PI * (std::pow(2, m) * x - k - 0.5)));
    };
}

}