#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/distributions.h"
#include "SWIFT/c_m_calculators.h"
#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"

BOOST_AUTO_TEST_SUITE(CNCalculatorsUT)

BOOST_AUTO_TEST_CASE(EuropeanPayoffUT)
{
    bool is_call = true;
    Swift::SwiftParameters params
    {
        1,
        -1,
        2,
        3,
        3
    };
    EuropeanOptionContract european_option_contract;

    std::size_t J = 3;
    int k1 = -1;
    int k2 = 2;
    std::size_t m = 1;
    double K = 80;
    int k = -1;
    std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << european_option_contract.GetPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
    k = 0;
    std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << european_option_contract.GetPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
    k = 1;
    std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << european_option_contract.GetPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
    k = 2;
    std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << european_option_contract.GetPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()