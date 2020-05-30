#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fstream>

#include "SWIFT/distributions.h"



BOOST_AUTO_TEST_SUITE(DistributionsUT)

BOOST_AUTO_TEST_CASE(HestonCharUT)
{
    double S_0 = 1.0;
    double K = 1.1;
    double r = 0.02;
    double q = 0;
    double tau = 15;
    double F_t_T = S_0 * exp((r - q) * tau);
    double x = std::log(F_t_T / K);
//    HestonParameters heston_parameters(/*k*/ 1.0, /*v_bar*/ 0.05, /*sigma*/ 0.2, /*rho*/ -0.7, /*v_0*/ 0.04);
    HestonParameters heston_parameters( 3.0, 0.1, 0.25, -0.8, 0.08);
    HestonDistribution underlying_geometry(heston_parameters, tau);

    std::vector<double> u_tests;
    for (std::size_t i = 0; i <= 40; ++i)
        u_tests.push_back(i * 0.1);

    for (double u : u_tests)
    {
//        std::complex<double> char_position_factor = GetCharPositionFactor(x, u);
        std::complex<double> char_position_factor = x / x;
        std::complex<double> cui_char_val = underlying_geometry.GetCuiChar(u, tau) * char_position_factor;
        std::complex<double> gatheral_char_val = underlying_geometry.GetGatheralChar(u, tau) * char_position_factor;
        BOOST_CHECK_CLOSE(cui_char_val.real(), gatheral_char_val.real(), 1e-9);
        if (u < 1.0)
        {
            std::complex<double> heston_char_val = underlying_geometry.GetHestonChar(u, tau) * char_position_factor;
            double hes_gath_diff = heston_char_val.real() - gatheral_char_val.real();
            if (IsZero(u))
                BOOST_CHECK_MESSAGE(std::isnan(heston_char_val.real()),
                                    "U:" << u << " Heston Char:" << heston_char_val.real() << " Should be nan.");
            else
                BOOST_CHECK_MESSAGE(IsZero(hes_gath_diff),
                                    "U:" << u << " Heston Char:" << heston_char_val.real() << " Gatheral Char:" << gatheral_char_val.real() << " Diff:" << hes_gath_diff);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()