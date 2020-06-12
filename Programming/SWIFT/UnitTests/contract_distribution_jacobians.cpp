#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/distributions.h"
#include "SWIFT/known_distribution_contract_combinations.h"
#include "SWIFT/my_math.h"
#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"

BOOST_AUTO_TEST_SUITE(GradientPricingUT)

BOOST_AUTO_TEST_CASE(PutCallParityBSUT)
{
    const double F = 100;
    const double years = 0.5;
    const double r = 0.3;
    const double q = 0.0;
//    const bool is_call = true;

    const double vol = 0.25;
    GBM underlying_geometry(vol, years);
    std::vector<double> strikes = { 80 };

    Swift::SwiftParameters swift_params{5, -31, 32, 8, 7};

    EuropeanOptionContract underlying_contract;
    for (double K : strikes)
    {
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        double heston_price_put = swift_evaluator.GetPrice(F, K, r, q, false);
        double heston_price_call = swift_evaluator.GetPrice(F, K, r, q, true);

        BOOST_CHECK_CLOSE(heston_price_call - heston_price_put, F - K * std::exp(-r * years), MATH_EPSILON * 100);
        std::cout << "Put " << heston_price_put << " Call " << heston_price_call << std::endl;
    }
}



BOOST_AUTO_TEST_CASE(PutCallParityHestonUT)
{
    const double F = 100;
    const double years = 0.5;
    const double r = 0.3;
    const double q = 0.0;
//    const bool is_call = true;

    HestonParameters parameters{ 1.0, 0.05, 0.2, -0.7, 0.04};
    HestonDistribution underlying_geometry(parameters, years);
    std::vector<double> strikes = { 80 };

    Swift::SwiftParameters swift_params{5, -31, 32, 10, 10};

    EuropeanOptionContract underlying_contract;
    for (double K : strikes)
    {
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        double heston_price_put = swift_evaluator.GetPrice(F, K, r, q, false);
        double heston_price_call = swift_evaluator.GetPrice(F, K, r, q, true);
        double real_heston_call = GetHestonEuropeanPrice(parameters, F, K, r, years, q);

        BOOST_CHECK_CLOSE(heston_price_call - heston_price_put, F - K * std::exp(-r * years), MATH_EPSILON * 100);
        std::cout << "Put " << heston_price_put << " Call " << heston_price_call << std::endl;
        std::cout << " Real Call " << real_heston_call << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(PutCallParityUT)
{
    const double F = 1;
    const double years = 0.5;
//    const double r = 0.03;
    const double r = 0.0;
    const double q = 0.0;
//    const bool is_call = true;

    HestonParameters parameters{ 1.0, 0.05, 0.2, -0.7, 0.04};
    HestonDistribution underlying_geometry(parameters, years);

    double error = 1e-5;
    std::vector<double> strikes = { 1.2 };
    std::size_t m = 5;
    // Wavelet scale determination.
    while (m < 10)
    {
        double estimated_error = std::numeric_limits<double>::lowest();
        double bad_strike = 0.0;
        for (double strike : strikes)
        {
            double current_error = underlying_geometry.GetErrorInTermsOfM(m, HestonDistribution::GetXCompression(F, strike, r, q, years));
            if (current_error > estimated_error)
            {
                estimated_error = current_error;
                bad_strike = strike;
            }
            if (estimated_error > error)
                break;
        }
        if (estimated_error > error)
        {
            ++m;
            std::cout << "Error " << estimated_error << " higher than " << error << " for strike " << bad_strike << " increasing m to " << m << std::endl;
        }
        else
        {
            std::cout << "m " << m << " fulfills " << error << " for all strikes (worst strike has error " << estimated_error << std::endl;
            break;
        }
    }
    double c = 5 * underlying_geometry.GetDomainTruncationInitialValue();
    int k = std::lround(std::ceil(std::pow(2, m) * c));
    std::size_t iota_density = std::lround(std::ceil(std::log2(k * MY_PI)));

    std::cout << "Params c " << c << " k " << k << " iota " << iota_density << std::endl;

    Swift::SwiftParameters swift_params{m, 1 - k, k, iota_density, iota_density};
//    Swift::SwiftParameters swift_params{5, -17, 32, 8, 7};

    EuropeanOptionContract underlying_contract;
    for (double K : strikes)
    {
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        double heston_price_put = swift_evaluator.GetPrice(F, K, r, q, false);
        double heston_price_call = swift_evaluator.GetPrice(F, K, r, q, true);

        double real_heston_price_call = GetHestonEuropeanPrice(parameters, F, K, r, years, q);
        double real_heston_price_put = real_heston_price_call - F + K * std::exp(-r * years);
        BOOST_CHECK_CLOSE(heston_price_call - heston_price_put, F - K * std::exp(-r * years), MATH_EPSILON * 100);
        std::cout << "Put " << heston_price_put << " Call " << heston_price_call << std::endl;
        std::cout << "Real Put " << real_heston_price_put << " Real Call " << real_heston_price_call << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(HestonGradientUT)
{
    const double F = 1;
    const double years = 0.5;
//    const double r = 0.03;
    const double r = 0.0;
    const double q = 0.0;
//    const bool is_call = true;

    HestonParameters parameters{ 1.0, 0.05, 0.2, -0.7, 0.04};
//    HestonParameters parameters{ 0.0, std::pow(0.25, 2), 1.0, 0.0, std::pow(0.25, 2)};
//    HestonParameters parameters{ 1.0, 1.0, 1.0, 1.0, 1.0};
    HestonDistribution underlying_geometry(parameters, years);

    double error = 1e-5;
    std::vector<double> strikes = { 1 };
    std::size_t m = 3;
    // Wavelet scale determination.
    while (m < 10)
    {
        double estimated_error = std::numeric_limits<double>::lowest();
        double bad_strike = 0.0;
        for (double strike : strikes)
        {
            double current_error = underlying_geometry.GetErrorInTermsOfM(m, HestonDistribution::GetXCompression(F, strike, r, q, years));
            if (current_error > estimated_error)
            {
                estimated_error = current_error;
                bad_strike = strike;
            }
            if (estimated_error > error)
                break;
        }
        if (estimated_error > error)
        {
            ++m;
            std::cout << "Error " << estimated_error << " higher than " << error << " for strike " << bad_strike << " increasing m to " << m << std::endl;
        }
        else
        {
            std::cout << "m " << m << " fulfills " << error << " for all strikes (worst strike has error " << estimated_error << std::endl;
            break;
        }
    }
    double c = underlying_geometry.GetDomainTruncationInitialValue();
    int k = std::lround(std::ceil(std::pow(2, m) * c));
    std::size_t iota_density = std::lround(std::ceil(std::log2(k * MY_PI)));

    std::cout << "Params c " << c << " k " << k << " iota " << iota_density << std::endl;

//    Swift::SwiftParameters swift_params{m, 1 - k, k, iota_density, iota_density};
    Swift::SwiftParameters swift_params{5, -31, 32, 10, 10};

    EuropeanOptionContract underlying_contract;
    for (double K : strikes)
    {
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        double heston_price = GetHestonEuropeanPrice(parameters, F, K, r, years, q);
        double approx_price = swift_evaluator.GetPrice(F, K, r, q, true);
        BOOST_CHECK_CLOSE(approx_price, heston_price, MATH_EPSILON * 100);

//        HestonPriceGradient heston_price_gradient = GetHestonEuropeanJacobian(parameters, F, K, r, years, q);
//        std::vector<double> real_gradient = {heston_price_gradient.m_k, heston_price_gradient.m_v_bar, heston_price_gradient.m_sigma, heston_price_gradient.m_rho, heston_price_gradient.m_v0};
//        std::vector<double> approx_european_heston = swift_evaluator.GetGradient(F, K, r, q);
//        for (std::size_t param = 0; param < underlying_geometry.GetNumberOfParameters(); ++param)
//            BOOST_CHECK_CLOSE(approx_european_heston[param], real_gradient[param], MATH_EPSILON * 100);
    }
}

BOOST_AUTO_TEST_SUITE_END()