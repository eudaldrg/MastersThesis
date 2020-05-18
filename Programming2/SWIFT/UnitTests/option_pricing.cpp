#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fstream>

#include "SWIFT/distributions.h"
#include "SWIFT/my_math.h"


BOOST_AUTO_TEST_SUITE(OptionPricingUT)

//BOOST_AUTO_TEST_CASE(HestonUT)
//{
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
//                std::complex<double> f_hat = underlying_geometry.GetChar((2 * jj - 1) * Swift::PI * (1 << m) / N, x, years);
//                std::complex<double> cv = f_hat
//                                          * std::exp(1i * static_cast<double>(k)* Swift::PI* static_cast<double>(2 * jj - 1) / static_cast<double>(N));
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
//        //	std::complex<double> value = underlying_geometry.GetChar((2 * j + 1) * Swift::PI * (1 << m) / N);
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
//        BOOST_CHECK_CLOSE(heston_real, approx_heston, Swift::MATH_EPSILON);
//    }
//}

BOOST_AUTO_TEST_CASE(CashOrNothingBSUT)
{
    const double F = 100;
    const double years = 0.1;
    const double vol = 0.25;
    const double r = 0.1;
    const double q = 0.0;
    GBM underlying_geometry(vol, years, r, 0);

    struct TestParams{
        std::size_t m;
        double K;
        int k1;
        int k2;
        int iota;
        int iota_bar;
    };

    std::vector<TestParams> tests = {
        {1, 80, -1, 2, 4, 3},
        {2,120, -3, 4, 5, 4},
        {3, 120, -7, 4, 6, 5},
        {4, 80, -8, 164, 7, 6},
        {5, 80, -17, 32, 8, 7}};

    for (auto const& [m, K, k1, k2, iota, iota_bar] : tests)
    {
        double cash_or_nothing = GetCashOrNothingPrice(K, F, r, years, vol);

//        double x = GetXCompression(F, K, r, q, years);
        (void)q;
        double x = std::log(F / K);

        std::size_t N = 1UL << iota;

        /// SLOW METHOD

        std::vector<double> c_n;
        const std::size_t two_to_the_j_minus_1 = 1UL << (iota - 1);
        for (int k = k1; k <= k2; ++k)
        {
            double c_m_k = 0;
            for (std::size_t jj = 1; jj <= two_to_the_j_minus_1; ++jj)
            {
                std::complex<double> f_hat = underlying_geometry.GetChar((2 * jj - 1) * Swift::PI * (1UL << m) / N, x);
                std::complex<double> cv = f_hat
                                          * std::exp(1i * static_cast<double>(k) * Swift::PI * static_cast<double>(2 * jj - 1) / static_cast<double>(N));
                c_m_k += cv.real();
            }
            c_n.push_back(c_m_k * std::pow(2.0, m / 2.0) / two_to_the_j_minus_1);
        }

        //fftw_complex* frequency_values, * time_values;
        //fftw_plan p;
        //frequency_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
        //time_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
        //p = fftw_plan_dft_1d(N, frequency_values, time_values, FFTW_BACKWARD, FFTW_ESTIMATE);
        //std::vector<std::complex<double>> x(N);
        //for (std::size_t i = 0; i < N; ++i)
        //{
        //	std::complex<double> value = underlying_geometry.GetChar((2 * j + 1) * Swift::PI * (1 << m) / N);
        //	frequency_values[i][0] = value.real();
        //	frequency_values[i][1] = value.imag();
        //}
        //fftw_execute(p);
        //std::vector<std::complex<double>> std_time_values;
        //for (std::size_t i = 0; i < N; ++i)
        //	std_time_values.emplace_back(time_values[i][0], time_values[i][1]);

        double approx_cash_or_nothing = 0.0;
        for (int k = k1; k <= k2; ++k)
        {
            double payoff = GetCashOrNothingCallPayOff(m, k, iota_bar);
            approx_cash_or_nothing += c_n[k - k1] * payoff;
//            std::cout << "k " << k << " payoff " << payoff << " shannon wavelet " << c_n[k - k1] << std::endl;
        }
        approx_cash_or_nothing *= std::exp(-r * years);
        BOOST_CHECK_CLOSE(cash_or_nothing, approx_cash_or_nothing, Swift::MATH_EPSILON);
    }
}

//BOOST_AUTO_TEST_CASE(BSUT)
//{
//    const double F = 100;
//    const double years = 1;
//    const double r = 0.1;
//    //const double q = 0.0;
//    const bool is_call = true;
//    const double vol = 0.25;
//    GBM underlying_geometry(vol, years, r, 0);
//    for (double K : std::vector<double>{100})
//    {
//        double real_value = 14.975790;
//        std::cout << "Test K = " << K << " : (" << real_value << ", " << 1 - real_value << ")" << std::endl;
//
//        double x = std::log(F / K);
//
//        std::size_t m = 3;
//        int k1 = -19;
//        int k2 = 20;
//        std::size_t j = 11;
//        std::size_t N = 1 << j;
//        std::size_t j_bar = 10;
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
//                std::complex<double> f_hat = underlying_geometry.GetChar((2 * jj - 1) * Swift::PI * (1 << m) / N, x);
//                std::complex<double> cv = f_hat
//                                          * std::exp(1i * static_cast<double>(k) * Swift::PI * static_cast<double>(2 * jj - 1) / static_cast<double>(N));
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
//        //	std::complex<double> value = underlying_geometry.GetChar((2 * j + 1) * Swift::PI * (1 << m) / N);
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
//        std::cout << "Original " << real_value << " Approx " << approx_heston << " Diff " << real_value - approx_heston << std::endl;
//    }
//}

BOOST_AUTO_TEST_SUITE_END()