#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"
#include "SWIFT/payoff_coefficients_calculators.h"

BOOST_AUTO_TEST_SUITE(PayoffCalculatorsUT)

///c_m_k = 2^(m/2) * int_{-1/2}^{1/2} (Re(char(2^{m+1} pi t) * exp(i 2 pi k t) dt
BOOST_AUTO_TEST_CASE(FastCoeffsTest)
{
//    const double years = 0.1;
    EuropeanOptionContract contract;
    Swift::SwiftParameters params{
        4,
//        5,
        -7,
        8,
//        84,
        5,
        5
//        7
    };

//    Swift::DCTAndDSTPayoffCalculator dct_payoff_calculator(contract, params);
    Swift::FFTPayoffCalculator fft_payoff_calculator(contract, params);
    Swift::ExplicitVietaPayoffCalculator vieta_payoff_calculator(contract, params);
    Swift::ExplicitPayoffCalculator explicit_payoff_calculator(contract, params);

//    std::vector<double> dct_payoffs = dct_payoff_calculator.GetNonKPayoffCoefs(true);
    std::vector<double> fft_payoffs = fft_payoff_calculator.GetNonKPayoffCoefs(true);
    std::vector<double> explicit_payoffs = explicit_payoff_calculator.GetNonKPayoffCoefs(true);
    std::vector<double> explicit_vieta_payoffs = vieta_payoff_calculator.GetNonKPayoffCoefs(true);
    std::vector<double> original_payoffs;

    for (int k = params.m_k1; k <= params.m_k2; ++k)
    {
        std::size_t index = k - params.m_k1;
        original_payoffs.push_back(contract.GetPayoffNonKComponent(k, params, true));
//        BOOST_CHECK_CLOSE(dct_payoffs[index], explicit_payoffs[index], 1e-7);
        BOOST_CHECK_CLOSE(explicit_vieta_payoffs[index], explicit_payoffs[index], 1e-7);
        BOOST_CHECK_CLOSE(fft_payoffs[index], explicit_vieta_payoffs[index], 1e-7);
        BOOST_CHECK_CLOSE(original_payoffs[index], explicit_payoffs[index], 1e-7);
    }
}

BOOST_AUTO_TEST_SUITE_END()