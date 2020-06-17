#include "swift.h"

#include <iostream>

#include "distributions.h"
#include "option_contracts.h"
#include "density_coefficients_calculators.h"
#include "payoff_coefficients_calculators.h"

namespace Swift {

SwiftParameters::SwiftParameters(size_t m, int k_1, int k_2, size_t iota_density, size_t iota_payoff) : m_m(m), m_k1(k_1), m_k2(k_2), m_iota_density(iota_density), m_iota_payoff(iota_payoff)
{
    // Iota bar is used for generating 2^(iota_payoff) frequency samples for the char_function Fourier transform.
    double log2_k_range = std::log2(m_k2 - m_k1);
    int ceil_log_2_k_range = std::lround(std::ceil(log2_k_range));
    if (m_iota_density < ceil_log_2_k_range)
    {
        std::cout << __FUNCTION__ << ": Wrong params. iota_density " << m_iota_density << " not big enough to apply FFT. Increasing to " << ceil_log_2_k_range << std::endl;
        m_iota_density = ceil_log_2_k_range;
    }

    // TODO: Introduce this here somehow.
//    // Cash-or-nothing
//    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
//    {
//        if (k == 0)
//            continue;
//        int required_iota = std::roundl(std::ceil(std::log2(MY_PI * std::abs(k))));
//        if (m_params.m_iota_payoff < required_iota)
//        {
//            std::cout << "SwiftEvaluator: iota_payoff " << m_params.m_iota_payoff << " is not sufficient for k1, k2 range. Updating to " << required_iota << std::endl;
//            m_params.m_iota_payoff = std::max(m_params.m_iota_payoff, required_iota);
//        }
//    }
//
//    // European
//    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
//    {
//        if (k == 0)
//            continue;
//        int required_iota = std::roundl(std::ceil(std::log2(MY_PI * std::max(std::abs(std::max(m_params.m_k1, 0) - k), std::abs(m_params.m_k2 - k)))));
//        if (m_params.m_iota_payoff < required_iota)
//        {
//            std::cout << "SwiftEvaluator: iota_payoff " << m_params.m_iota_payoff << " is not sufficient for k1, k2 range. Updating to " << required_iota << std::endl;
//            m_params.m_iota_payoff = std::max(m_params.m_iota_payoff, required_iota);
//        }
//    }

    m_two_to_the_m = std::pow(2, m_m);
    m_payoff_from = static_cast<double>(m_k1) / m_two_to_the_m;
    m_payoff_to = static_cast<double>(m_k2) / m_two_to_the_m;
    m_sqrt_two_to_the_m = std::sqrt(static_cast<double>(m_two_to_the_m));
    m_J_density = std::pow(2, m_iota_density);
    m_N_density = 2 * m_J_density;
    m_J_payoff = std::pow(2, m_iota_payoff);
    m_N_payoff = 2 * m_J_payoff;
}

SwiftEvaluator::SwiftEvaluator(SwiftParameters const& params, Distribution const& distribution, OptionContract const& option_contract)
    : m_params(params), m_distribution(distribution), m_option_contract(option_contract) { }

//double SwiftEvaluator::GetPrices(double S, std::vector<double> Ks, double r, double q, bool is_call) const
//{
//    double inverse_exp_interest_times_years = std::exp(-r * m_distribution.m_T);
//    if (is_call)
//    {
//        double put_price = GetPrice(S, K, r, q, !is_call);
//        return put_price + S - K * inverse_exp_interest_times_years;
//    }
//
//    double x = Distribution::GetXCompression(S, K, r, q, m_distribution.m_T);
//
//    FastVietaCalculator calculator(m_distribution, m_params);
////    ExplicitVietaCalculator calculator(m_distribution, m_params);
//    std::vector<double> c_m = calculator.GetCMCoefs(x);
//
////    double coefficients_sum = std::accumulate(c_m.begin(), c_m.end(), 0.0);
////    double theoretic_error = std::abs(1 - std::pow(2, -0.5 * m_params.m_m) * coefficients_sum);
////    std::cout << "Theoretic error should be under " << theoretic_error << std::endl;
//
//    double option_price = 0.0;
//    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
//    {
//        double payoff = m_option_contract.GetPayoff(K, m_params.m_m, k, m_params.m_k1, m_params.m_k2, m_params.m_iota_payoff, is_call);
//        option_price += c_m[k - m_params.m_k1] * payoff;
//    }
//    return option_price * inverse_exp_interest_times_years;
//}

double SwiftEvaluator::GetPrice(double S, double K, double r, double q, bool is_call) const
{
    double inverse_exp_interest_times_years = std::exp(-r * m_distribution.m_T);
//    if (is_call)
//    {
////////        double prev_k1 = m_params.m_k1;
////////        double prev_k2 = m_params.m_k2;
////////        m_params.m_k1 = -prev_k2;
////////        m_params.m_k2 = -prev_k1;
//        double put_price = GetPrice(S, K, r, q, !is_call);
////////        m_params.m_k1 = prev_k1;
////////        m_params.m_k2 = prev_k2;
//        return put_price + S - K * inverse_exp_interest_times_years;
//    }

    double x = Distribution::GetXCompression(S, K, r, q, m_distribution.m_T);

    FastVietaCalculator calculator(m_distribution, m_params);
//    ExplicitVietaCalculator calculator(m_distribution, m_params);
    std::vector<double> c_m = calculator.GetCMCoefs(x);

//    FFTPayoffCalculator payoff_calculator(m_option_contract, m_params);
    ExplicitPayoffCalculator payoff_calculator(m_option_contract, m_params);
    std::vector<double> non_k_payoff_coefficients = payoff_calculator.GetNonKPayoffCoefs(is_call);

//    double coefficients_sum = std::accumulate(c_m.begin(), c_m.end(), 0.0);
//    double theoretic_error = std::abs(1 - std::pow(2, -0.5 * m_params.m_m) * coefficients_sum);
//    std::cout << "Theoretic error should be under " << theoretic_error << std::endl;

    double option_price = 0.0;
    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
        option_price += c_m[k - m_params.m_k1] * non_k_payoff_coefficients[k - m_params.m_k1];
    return option_price * inverse_exp_interest_times_years * K;
}

std::vector<double> SwiftEvaluator::GetGradient(double F, double K, double r, double q) const
{
    double x = Distribution::GetXCompression(F, K, r, q, m_distribution.m_T);

    FastVietaCalculator fast_vieta_calculator(m_distribution, m_params);
    std::vector<std::vector<double>> c_m = fast_vieta_calculator.GetGradientCMCoefs(x);

    std::vector<double> option_gradient(m_distribution.GetNumberOfParameters(), 0.0);
    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
    {
        double payoff = m_option_contract.GetPayoff(K, k, m_params, true);
        for (std::size_t param = 0; param < m_distribution.GetNumberOfParameters(); ++param)
            option_gradient[param] += c_m[k - m_params.m_k1][param] * payoff;
    }
    for (std::size_t param = 0; param < m_distribution.GetNumberOfParameters(); ++param)
        option_gradient[param] *= std::exp(-r * m_distribution.m_T); // TODO: Add q or remove altogether.
    return option_gradient;
}

}