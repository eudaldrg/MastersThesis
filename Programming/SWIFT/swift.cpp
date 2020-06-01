#include "swift.h"

#include <iostream>

#include "distributions.h"
#include "option_contracts.h"
#include "c_m_calculators.h"

namespace Swift {

SwiftParameters::SwiftParameters(size_t m, int k_1, int k_2, size_t iota, size_t iota_bar) : m_m(m), m_k1(k_1), m_k2(k_2), m_iota_dens_coefs(iota), m_iota_payoff_coefs(iota_bar)
{
    // Iota bar is used for generating 2^(iota_bar) frequency samples for the char_function Fourier transform.
    double log2_k_range = std::log(m_k2 - m_k1) / std::log(2);
    std::size_t ceil_log_2_k_range = std::lround(std::ceil(log2_k_range));
    if (m_iota_payoff_coefs < ceil_log_2_k_range)
    {
        std::cout << __FUNCTION__ << ": Wrong params. iota_bar " << m_iota_payoff_coefs << " not big enough. Increasing to " << ceil_log_2_k_range << std::endl;
        m_iota_payoff_coefs = ceil_log_2_k_range;
    }
}

SwiftEvaluator::SwiftEvaluator(SwiftParameters const& params, Distribution const& distribution, OptionContract const& option_contract)
    : m_params(params), m_distribution(distribution), m_option_contract(option_contract)
{
    // Cash-or-nothing
    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
    {
        if (k == 0)
            continue;
        std::size_t required_iota = std::log(MY_PI * std::abs(k)) / std::log(2);
        if (m_params.m_iota_payoff_coefs < required_iota)
        {
            std::cout << "SwiftEvaluator: iota_bar " << m_params.m_iota_payoff_coefs << " is not sufficient for k1, k2 range. Updating to " << required_iota << std::endl;
            m_params.m_iota_payoff_coefs = std::max(m_params.m_iota_payoff_coefs, required_iota);
        }
    }

    // European
    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
    {
        if (k == 0)
            continue;
        std::size_t required_iota = std::log(MY_PI * std::max(std::abs(std::max(m_params.m_k1, 0) - k), std::abs(m_params.m_k2 - k))) / std::log(2);
        if (m_params.m_iota_payoff_coefs < required_iota)
        {
            std::cout << "SwiftEvaluator: iota_bar " << m_params.m_iota_payoff_coefs << " is not sufficient for k1, k2 range. Updating to " << required_iota << std::endl;
            m_params.m_iota_payoff_coefs = std::max(m_params.m_iota_payoff_coefs, required_iota);
        }
    }
}

double SwiftEvaluator::GetPrice(double F, double K, double r, double q) const
{
    double x = Distribution::GetXCompression(F, K, r, q, m_distribution.m_T);

    FastVietaCalculator fast_vieta_calculator(m_distribution, m_params);
    std::vector<double> c_m = fast_vieta_calculator.GetCMCoefs(x);

    double option_price = 0.0;
    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
    {
        double payoff = m_option_contract.GetPayoff(K, m_params.m_m, k, m_params.m_k1, m_params.m_k2, m_params.m_iota_payoff_coefs, true);
        option_price += c_m[k - m_params.m_k1] * payoff;
    }
    return option_price * std::exp(-r * m_distribution.m_T);
}

}