#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/distributions.h"
#include "SWIFT/known_distribution_contract_combinations.h"
#include "SWIFT/my_math.h"
#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"

BOOST_AUTO_TEST_SUITE(OptionPricingUT)

BOOST_AUTO_TEST_CASE(HestonGradientUT)
{
    const double F = 80;
    const double years = 0.5;
    const double r = 0.03;
    const double q = 0.0;
//    const bool is_call = true;

    HestonParameters parameters{ 1.0, 0.05, 0.2, -0.7, 0.04};
    HestonDistribution underlying_geometry(parameters, years);
    EuropeanOptionContract underlying_contract;
    for (double K : std::vector<double>{ 80 })
    {
        Swift::SwiftParameters swift_params{5, -17, 32, 8, 7};
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        HestonPriceGradient heston_price_gradient = GetHestonEuropeanJacobian(parameters, F, K, r, years, q);
        std::vector<double> real_gradient = {heston_price_gradient.m_k, heston_price_gradient.m_v_bar, heston_price_gradient.m_sigma, heston_price_gradient.m_rho, heston_price_gradient.m_v0};
        std::vector<double> approx_european_heston = swift_evaluator.GetGradient(F, K, r, q);
        for (std::size_t param = 0; param < underlying_geometry.GetNumberOfParameters(); ++param)
            BOOST_CHECK_CLOSE(approx_european_heston[param], real_gradient[param], MATH_EPSILON * 100);
    }
}

BOOST_AUTO_TEST_SUITE_END()