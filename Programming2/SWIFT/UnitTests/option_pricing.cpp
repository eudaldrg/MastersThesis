#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fstream>

#include "SWIFT/distributions.h"
#include "FFTW3/include_fftw3.h"
#include "SWIFT/known_distribution_contract_combinations.h"
#include "SWIFT/my_math.h"
#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"

BOOST_AUTO_TEST_SUITE(OptionPricingUT)

BOOST_AUTO_TEST_CASE(HestonUT)
{
//    const double F = 80;
//    const double years = 0.5;
//    const double r = 0.03;
//    const double q = 0.0;
//    const bool is_call = true;
//    HestonParameters heston_parameters(/*k*/ 1.0, /*v_bar*/ 0.05, /*sigma*/ 0.2, /*rho*/ -0.7, /*v_0*/ 0.04);
//    HestonDistribution underlying_geometry(heston_parameters);
//    for (double K : std::vector<double>{ 80 })
//    {
//        double heston_real = 5.1772;
//        std::cout << "Test K = " << K << " : (" << heston_real << ", " << 1 - heston_real << ")" << std::endl;
//
//        double x = GetXCompression(F, K, r, q, years);
//
//        //std::size_t m = 1;
//        std::size_t m = 4;
//        //int k1 = -1;
//        int k1 = -8;
//        //int k2 = 2;
//        int k2 = 16;
//        //std::size_t j = 4;
//        std::size_t j = 7;
//        std::size_t N = 1 << j;
//        //std::size_t j_bar = 3;
//        std::size_t j_bar = 6;
//
//        /// SLOW METHOD
//
//        std::vector<double> c_n;
//        const std::size_t two_to_the_j_minus_1 = 1 << (j - 1);
//        for (int k = k1; k <= k2; ++k)
//        {
//            double c_m_k = 0;
//            for (std::size_t jj = 1; jj <= two_to_the_j_minus_1; ++jj)
//            {
//                std::complex<double> f_hat = underlying_geometry.GetChar((2 * jj - 1) * MY_PI * (1 << m) / N, x, years);
//                std::complex<double> cv = f_hat
//                                          * std::exp(1i * static_cast<double>(k)* MY_PI* static_cast<double>(2 * jj - 1) / static_cast<double>(N));
//                c_m_k += cv.real();
//            }
//            c_n.push_back(c_m_k * std::pow(2.0, m / 2.0) / two_to_the_j_minus_1);
//        }
//
//        //fftw_complex* frequency_values, * time_values;
//        //fftw_plan p;
//        //frequency_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
//        //time_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
//        //p = fftw_plan_dft_1d(N, frequency_values, time_values, FFTW_BACKWARD, FFTW_ESTIMATE);
//        //std::vector<std::complex<double>> x(N);
//        //for (std::size_t i = 0; i < N; ++i)
//        //{
//        //	std::complex<double> value = underlying_geometry.GetChar((2 * j + 1) * MY_PI * (1 << m) / N);
//        //	frequency_values[i][0] = value.real();
//        //	frequency_values[i][1] = value.imag();
//        //}
//        //fftw_execute(p);
//        //std::vector<std::complex<double>> std_time_values;
//        //for (std::size_t i = 0; i < N; ++i)
//        //	std_time_values.emplace_back(time_values[i][0], time_values[i][1]);
//        double approx_heston = 0.0;
//        for (int k = k1; k <= k2; ++k)
//        {
//            double payoff = GetEuropeanPayoff(K, m, k, k1, k2, j_bar, is_call);
//            approx_heston += c_n[k - k1] * payoff;
//            std::cout << "k " << k << " payoff " << payoff << " shannon wavelet " << c_n[k - k1] << std::endl;
//        }
//        approx_heston *= std::exp(-r * years);
//        BOOST_CHECK_CLOSE(heston_real, approx_heston, MATH_EPSILON);
//    }
}

BOOST_AUTO_TEST_CASE(CashOrNothingBSUTOrtizGraciaPaper)
{
    const double F = 100;
    const double years = 0.1;
    const double vol = 0.25;
    const double r = 0.1;
    const double q = 0.0;
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
        double cash_or_nothing = GetBSCashOrNothingPrice(K, F, r, years, vol);
        Swift::SwiftEvaluator swift_evaluator(test_params.m_swift_parameters, underlying_geometry, underlying_contract);
        double approx_cash_or_nothing = swift_evaluator.GetPrice(F, K, r, q);
        BOOST_CHECK_CLOSE(test_params.reference_value, cash_or_nothing, MATH_EPSILON * 100);
//        std::cout << std::abs(approx_cash_or_nothing - cash_or_nothing) << std::endl;
        BOOST_CHECK(IsZero(approx_cash_or_nothing - cash_or_nothing, test_params.error));
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
    EuropeanOptionContract underlying_contract;

    struct TestParams{
        double K;
        Swift::SwiftParameters m_parameters;
    };

    std::vector<TestParams> tests = {
        {80, {1, -1, 2, 4, 3}},
        {120, {2,-3, 2, 5, 4}},
        {120, {3, -7, 4, 6, 5}},
        {80, {4, -8, 16, 7, 6}},
        {80, {5, -17, 32, 8, 7}},
//        {80, {5, -38, 48, 7, 6}},
//        {120, {5, -20, 20, 8, 7}},
    };

    for (auto& [K, swift_params] : tests)
    {
        Swift::SwiftEvaluator swift_evaluator(swift_params, underlying_geometry, underlying_contract);

        double european_price_bs = GetBSEuropeanPrice(K, F, r, years, vol, is_call);

        double approx_european_bs = swift_evaluator.GetPrice(F, K, r, q);
        BOOST_CHECK_CLOSE(approx_european_bs, european_price_bs, MATH_EPSILON * 100);
    }
}

BOOST_AUTO_TEST_SUITE_END()