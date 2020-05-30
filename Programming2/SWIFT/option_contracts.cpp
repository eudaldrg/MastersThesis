#include "option_contracts.h"

double EuropeanOptionContract::GetPayoff(double K, std::size_t m, int k, int k1, int k2, std::size_t iota, bool is_call) const
{
    double return_value = 0.0;
    if ((is_call && k2 <= 0) || (!is_call && k1 >= 0))
    return return_value;

    int k1_bar = std::max(k1, 0);
    int k2_bar = std::min(k2, 0);

    std::size_t two_to_the_iota_minus_one = 1 << (iota - 1);
    double two_to_the_m = 1 << m;
    for (std::size_t j = 1; j <= two_to_the_iota_minus_one; ++j)
    {
        double I1 = is_call ? GetI1(k1_bar / two_to_the_m, k2 / two_to_the_m, j, iota, m, k) : GetI1(k1 / two_to_the_m, k2_bar / two_to_the_m, j, iota, m, k);
        double I2 = is_call ? GetI2(k1_bar / two_to_the_m, k2 / two_to_the_m, j, iota, m, k) : GetI2(k1 / two_to_the_m, k2_bar / two_to_the_m, j, iota, m, k);
        assert(std::isnormal(I1) || I1 == 0.0);
        assert(std::isnormal(I2) || I2 == 0.0);
        return_value += I1 - I2;
    }
    return_value *= K * std::sqrt(two_to_the_m) / two_to_the_iota_minus_one;
    if (!is_call)
    return_value *= -1.0;
    return return_value;
}

double EuropeanOptionContract::GetC_j(std::size_t j, std::size_t J) const
{
    return (2 * j - 1) * MY_PI / (1 << J);
}

double EuropeanOptionContract::GetI1(double a, double b, std::size_t j, std::size_t J, std::size_t m, int k) const
{
    double two_to_the_m = 1 << m;
    double C_j = GetC_j(j, J);
    double C_j_times_two_to_the_m = C_j * two_to_the_m;

    double C_j_b = C_j * (two_to_the_m * b - k);
    double C_j_a = C_j * (two_to_the_m * a - k);
    double common_factor = (C_j_times_two_to_the_m) / (1 + C_j_times_two_to_the_m * C_j_times_two_to_the_m);
    return common_factor * (std::exp(b) * std::sin(C_j_b) - std::exp(a) * std::sin(C_j_a) + (std::exp(b) * std::cos(C_j_b) - std::exp(a) * std::cos(C_j_a)) / C_j_times_two_to_the_m);
}

double EuropeanOptionContract::GetI2(double a, double b, std::size_t j, std::size_t J, std::size_t m, int k) const
{
    double two_to_the_m = 1 << m;
    double C_j = GetC_j(j, J);
    double C_j_times_two_to_the_m = C_j * two_to_the_m;

    double C_j_b = C_j * (two_to_the_m * b - k);
    double C_j_a = C_j * (two_to_the_m * a - k);

    double return_value = (std::sin(C_j_b) - std::sin(C_j_a)) / C_j_times_two_to_the_m;
    assert(std::isnormal(return_value) || return_value == 0.0);
    return return_value;
}

double CashOrNothingContract::GetPayoff(double /*K*/, std::size_t m, int k, int /*k1*/, int /*k2*/, std::size_t iota, bool is_call) const
{
//        const int sign = (k == 0) ? 0 : ((k > 0) ? 1 : -1);
//        return std::pow(2.0, m / -2.0) * (sign * SIApprox(std::abs(k), J) + 0.5);
    if (!is_call)
        throw std::runtime_error(std::string(__FUNCTION__) + ": Not implemented for puts.");
    return (Sign(k) * SIApprox(std::abs(k), iota) + 0.5) / std::pow(2.0, m / 2.0);
}