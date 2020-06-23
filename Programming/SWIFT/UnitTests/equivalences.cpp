#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <SWIFT/known_distribution_contract_combinations.h>
#include <SWIFT/swift.h>
#include <SWIFT/option_contracts.h>
#include <LEVMAR/levmar.h>

#include "SWIFT/distributions.h"

struct Data
{
    std::vector<double> Ks = {
        0.9371, 0.8603, 0.8112, 0.7760, 0.7470, 0.7216, 0.6699, 0.6137,
        0.9956, 0.9868, 0.9728, 0.9588, 0.9464, 0.9358, 0.9175, 0.9025,
        1.0427, 1.0463, 1.0499, 1.0530, 1.0562, 1.0593, 1.0663, 1.0766,
        1.2287, 1.2399, 1.2485, 1.2659, 1.2646, 1.2715, 1.2859, 1.3046,
        1.3939, 1.4102, 1.4291, 1.4456, 1.4603, 1.4736, 1.5005, 1.5328};
    std::vector<double> Ts = {
        0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
        0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
        0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
        0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
        0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143
    };
    std::vector<double> heston_matlab = {
        0.0691, 0.1500, 0.1995, 0.2390, 0.2691, 0.2981, 0.3556, 0.4157,
        0.0265, 0.0515, 0.0675, 0.0894, 0.1035, 0.1208, 0.1519, 0.1787,
        0.0079, 0.0234, 0.0287, 0.0396, 0.0439, 0.0527, 0.0686, 0.0814,
        0.0000, 0.0002, 0.0006, 0.0018, 0.0030, 0.0055, 0.0118, 0.0185,
        0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0002, 0.0010, 0.0024
    };
    double S = 1.0;
    double r = 0.02;
    double q = 0.0;
    HestonParameters hes_parameters{3.0, 0.10, 0.25, -0.8, 0.08};
    GBMParameters bs_parameters{0.25};
    Swift::SwiftParameters swift_params{3, -31, 32, 10, 10};
    EuropeanOptionContract european_contract;
};

BOOST_FIXTURE_TEST_SUITE( EquivalencesUT,  Data)

BOOST_AUTO_TEST_CASE(Matlab_EQUIVALENT_TO_ORIGINAL_CUI_UT)
{
    for (std::size_t i = 0; i < Ks.size(); ++i)
    {
        double matlab = heston_matlab[i];
        double cui = GetHestonEuropeanPrice(hes_parameters, S, Ks[i], r, Ts[i], 0.0);
        BOOST_CHECK_CLOSE(matlab, cui, MATH_EPSILON * 100);
        std::cout << "factor " << matlab / cui << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(MyCUI_EQUIVALENT_TO_ORIGINAL_CUI_UT)
{
    for (std::size_t i = 0; i < Ks.size(); ++i)
    {
        double my_cui = GetHestonEuropeanPriceCuiMyChar(hes_parameters, S, Ks[i], r, 0.0, Ts[i]);
        double cui = GetHestonEuropeanPrice(hes_parameters, S, Ks[i], r, Ts[i], 0.0);
        BOOST_CHECK_CLOSE(my_cui, cui, MATH_EPSILON * 100);
    }
}

BOOST_AUTO_TEST_CASE(CuiBSUT)
{
    for (std::size_t i = 0; i < Ks.size(); ++i)
    {
        double real = GetBSEuropeanPrice(Ks[i], S, r, Ts[i], bs_parameters.m_vol, true);
        double cui = GetBSEuropeanPriceCuiMyChar(bs_parameters, S, Ks[i], r, 0.0, Ts[i]);
        BOOST_CHECK_CLOSE(real, cui, MATH_EPSILON * 100);
        std::cout << "factor " << real / cui << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(SwiftBSUT)
{
    for (std::size_t i = 0; i < Ks.size(); ++i)
    {
        GBM distribution(bs_parameters.m_vol, Ts[i]);
        Swift::SwiftEvaluator eval(swift_params, distribution, european_contract);
        double real = GetBSEuropeanPrice(Ks[i], S, r, Ts[i], bs_parameters.m_vol, true);
        double swift = eval.GetPrice(S, Ks[i], r, 0.0, true);
        BOOST_CHECK_MESSAGE(std::abs(real - swift) < 1e-3, "Error " << real - swift);
//        std::cout << "factor " << real / swift << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(Jacobian_BS_UT)
{
    for (std::size_t i = 0; i < Ks.size(); ++i)
    {
        GBM distribution(bs_parameters.m_vol, Ts[i]);
        Swift::SwiftEvaluator eval(swift_params, distribution, european_contract);
        double real = GetBSEuropeanVega(Ks[i], S, r, Ts[i], bs_parameters.m_vol, q);
        double swift = eval.GetGradient(S, Ks[i], r, 0.0, true)[0];
        BOOST_CHECK_MESSAGE(std::abs(real - swift) < 1e-3, "Error " << real - swift);
        std::cout << "factor " << real / swift << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(Swift_Heston_UT)
{
    for (std::size_t i = 0; i < Ks.size(); ++i)
    {
        HestonDistribution distribution(hes_parameters, Ts[i]);
        Swift::SwiftEvaluator eval(swift_params, distribution, european_contract);
        double matlab = heston_matlab[i];
        double swift = eval.GetPrice(S, Ks[i], r, 0.0, true);
        BOOST_CHECK_CLOSE(matlab, swift, MATH_EPSILON * 100);
        std::cout << "factor " << matlab / swift << std::endl;
    }
}

BOOST_AUTO_TEST_SUITE_END()

