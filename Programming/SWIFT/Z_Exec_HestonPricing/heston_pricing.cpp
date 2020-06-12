#include <iostream>
#include <SWIFT/option_contracts.h>
#include <SWIFT/known_distribution_contract_combinations.h>
#include <SWIFT/swift.h>
#include "vector"
#include "SWIFT/distributions.h"

int main()
{
    const double F = 80;
    const double years = 0.5;
    const double r = 0.03;
    const double q = 0.0;
//    const bool is_call = true;

    HestonParameters parameters{ 1.0, 0.05, 0.2, -0.7, 0.04};
    HestonDistribution underlying_geometry(parameters, years);
    EuropeanOptionContract underlying_contract;

    double error = 1e-5;
    std::vector<double> strikes = { 80 };
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



    Swift::SwiftParameters swift_params{m, -17, 32, 8, 7};
    Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);



    for (double K : strikes)
    {
        // PRICE
        double heston_price = GetHestonEuropeanPrice(parameters, F, K, r, years, q);
        double approx_price = swift_evaluator.GetPrice(F, K, r, q, true);
        if (!IsSame(approx_price, heston_price, MATH_EPSILON * 100))
            std::cout << "approx_european_heston " << approx_price << " different than real price " << heston_price << std::endl;

//        // GRADIENT
//        HestonPriceGradient heston_price_gradient = GetHestonEuropeanJacobian(parameters, F, K, r, years, q);
//        std::vector<double> real_gradient = {heston_price_gradient.m_k, heston_price_gradient.m_v_bar, heston_price_gradient.m_sigma, heston_price_gradient.m_rho, heston_price_gradient.m_v0};
//        std::vector<double> approx_european_heston = swift_evaluator.GetGradient(F, K, r, q);
//        for (std::size_t param = 0; param < underlying_geometry.GetNumberOfParameters(); ++param)
//            if (!IsSame(approx_european_heston[param], real_gradient[param], MATH_EPSILON * 100))
//                std::cout << "approx_european_heston " << approx_european_heston[param] << " for param " << param << " different than real_gradient " << real_gradient[param] << std::endl;
    }
}
