#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/distributions.h"
#include "SWIFT/known_distribution_contract_combinations.h"
#include "SWIFT/my_math.h"
#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"

BOOST_AUTO_TEST_SUITE(OptionPricingUT)

BOOST_AUTO_TEST_CASE(TestMyChar)
{
    const double S = 100;
    const double T = 15;
    const double r = 0.03;
    const double q = 0.0;
//    const bool is_call = true;
    HestonParameters parameters{ 1.0, 0.05, 0.2, -0.7, 0.04};
    HestonDistribution distribution(parameters, T);
    EuropeanOptionContract underlying_contract;
    for (double K : std::vector<double>{ 80 })
    {
        Swift::SwiftParameters swift_params{3, -31, 32, 10, 10};
        Swift::SwiftEvaluator eval(swift_params, distribution, underlying_contract);

        double swift = eval.GetPrice(S, K, r, q, true);
        double real = 56.5072;
        double my_cui = GetHestonEuropeanPriceCuiMyChar(parameters, S, K, r, q, T);
        double cui = GetHestonEuropeanPrice(parameters, S, K, r, T, q);
        BOOST_CHECK_CLOSE(my_cui, cui, MATH_EPSILON * 100);
        BOOST_CHECK_CLOSE(swift, real, MATH_EPSILON * 100);
        std::cout << "real " << real << " my cui " << my_cui << " cui original " << cui << " swift " << swift << std::endl;
        std::cout << "factor " << real / my_cui << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(TestMyBSChar)
{
    const double S = 100;
    const double T = 15;
    const double r = 0.03;
    const double q = 0.0;
    const bool is_call = true;
    GBMParameters parameters{ 0.025 };
    GBM distribution(parameters.m_vol, T);
    EuropeanOptionContract contract;
    for (double K : std::vector<double>{ 80 })
    {
        Swift::SwiftParameters swift_params{5, -17, 32, 8, 7};
        Swift::SwiftEvaluator eval(swift_params, distribution, contract);

//        double inv_exp_r_times_t = std::exp(r * T);
        double my_cui = GetBSEuropeanPriceCuiMyChar(parameters, S, K, r, q, T);
        double real = GetBSEuropeanPrice(K, S, r, T, parameters.m_vol, true, q);
        double swift = eval.GetPrice(S, K, r, q, is_call);
        BOOST_CHECK_CLOSE(my_cui, real, MATH_EPSILON * 100);
        BOOST_CHECK_CLOSE(swift, real, MATH_EPSILON * 100);
        std::cout << my_cui - real << std::endl;
        std::cout << swift - real << std::endl;
        std::cout << "real " << real << " my cui " << my_cui << " swift " << swift << std::endl;
        std::cout << "factor " << real / my_cui << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(HestonUT)
{
    const double F = 80;
    const double years = 0.5;
    const double r = 0.03;
    const double q = 0.0;
    const bool is_call = true;
    HestonDistribution underlying_geometry({ 1.0, 0.05, 0.2, -0.7, 0.04}, years);
    EuropeanOptionContract underlying_contract;
    for (double K : std::vector<double>{ 80 })
    {
        Swift::SwiftParameters swift_params{5, -17, 32, 8, 7};
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        double european_price_heston = 5.1772;
        double approx_european_heston = swift_evaluator.GetPrice(F, K, r, q, is_call);
        BOOST_CHECK_CLOSE(approx_european_heston, european_price_heston, MATH_EPSILON * 100);
    }
}

BOOST_AUTO_TEST_CASE(CashOrNothingBSUTOrtizGraciaPaper)
{
    const double F = 100;
    const double years = 0.1;
    const double vol = 0.25;
    const double r = 0.1;
    const double q = 0.0;
    const bool is_call = true;
    GBM underlying_geometry(vol, years);
    CashOrNothingContract underlying_contract;

    struct TestParams{
        double K;
        Swift::SwiftParameters m_swift_parameters;
        double reference_value;
        double error;
    };

    std::vector<TestParams> tests = {
        {80, {1, -1, 2, 4, 3}, 0.988257979564030, 1.93e-1},
        {120, {2,-3, 2, 5, 4}, 0.0131034102155745, 4.42e-2},
        {120, {3, -7, 4, 6, 5}, 0.0131034102155745, 1.06e-2},
        {80, {4, -8, 16, 7, 6}, 0.988257979564030, 6.36e-6},
        {80, {5, -17, 32, 8, 7}, 0.988257979564030, 3.33e-16}
    };

    for (TestParams const& test_params : tests)
    {
        double K = test_params.K;
        double cash_or_nothing = GetBSCallCashOrNothingPrice(K, F, r, years, vol);
        Swift::SwiftEvaluator swift_evaluator(test_params.m_swift_parameters, underlying_geometry, underlying_contract);
        double approx_cash_or_nothing = swift_evaluator.GetPrice(F, K, r, q, is_call);
        BOOST_CHECK_CLOSE(test_params.reference_value, cash_or_nothing, MATH_EPSILON * 100);
//        std::cout << std::abs(approx_cash_or_nothing - cash_or_nothing) << std::endl;
        BOOST_CHECK_MESSAGE(std::abs(approx_cash_or_nothing - cash_or_nothing) < test_params.error, "Approx " << approx_cash_or_nothing << " Real " << cash_or_nothing
            << " diff " << approx_cash_or_nothing - cash_or_nothing << " expected " << test_params.error);
    }
}

BOOST_AUTO_TEST_CASE(BSUT)
{
    const double F = 100;
    const double years = 0.1;
    const double r = 0.1;
    const double q = 0.0;
    const bool is_call = true;
    const double vol = 0.25;
    GBM underlying_geometry(vol, years);
    GBMDerivative underlying_geometry_der(vol, years);
    EuropeanOptionContract underlying_contract;

    struct TestParams{
        double K;
        Swift::SwiftParameters m_parameters;
    };

    std::vector<TestParams> tests = {
//        {80, {1, -17, 32, 8, 7}},
//        {120, {2,-3, 2, 5, 4}},
//        {120, {3, -17, 18, 8, 7}},
//        {120, {4, -17, 18, 6, 5}},
//        {80, {3, -17, 18, 8, 7}},
        {80, {4, -17, 18, 8, 7}},
        {80, {5, -17, 18, 8, 7}},
        {80, {5, -38, 48, 7, 6}},
//        {120, {5, -20, 20, 8, 7}},
    };

    for (auto& [K, swift_params] : tests)
    {
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        double european_price_bs = GetBSEuropeanPrice(K, F, r, years, vol, is_call);
        double approx_european_bs = swift_evaluator.GetPrice(F, K, r, q, is_call);
        BOOST_CHECK_CLOSE(approx_european_bs, european_price_bs, MATH_EPSILON * 100);

        Swift::SwiftEvaluator swift_derivative(swift_params, underlying_geometry_der, underlying_contract);

        double european_vega_bs = GetBSEuropeanVega(K, F, r, years, vol, q);
        double approx_european_vega_bs = swift_derivative.GetPrice(F, K, r, q, is_call);
        BOOST_CHECK_CLOSE(european_vega_bs, approx_european_vega_bs, MATH_EPSILON * 100);
    }
}

BOOST_AUTO_TEST_SUITE_END()