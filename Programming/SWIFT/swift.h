#pragma once

#include <functional>
#include <numeric>

#include "my_math.h"

class Distribution;
class OptionContract;

namespace Swift {

class SwiftParameters
{
public:
    SwiftParameters(size_t m, int k_1, int k_2, size_t iota_density, size_t iota_payoff);

    int m_m;
    int m_two_to_the_m;
    double m_sqrt_two_to_the_m;
    int m_k1;
    double m_payoff_from;
    int m_k2;
    double m_payoff_to;
    int m_J_density;
    int m_N_density;
    int m_J_payoff;
    int m_N_payoff;
private:
    int m_iota_density; // Iota for the C_m_k calculation. Determines the number of intervals in which we subdivide the sinc integral.
    int m_iota_payoff; // Iota for the V_m_k calculation. (determines the truncation of the integration space).
};

class SwiftEvaluator
{
public:
    SwiftEvaluator(SwiftParameters const& params, Distribution const& distribution, OptionContract const& option_contract);

    [[nodiscard]] double GetPrice(double S, double K, double r, double q, bool is_call) const;
    [[nodiscard]] std::vector<double> GetGradient(double F, double K, double r, double q) const;

    mutable SwiftParameters m_params;
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