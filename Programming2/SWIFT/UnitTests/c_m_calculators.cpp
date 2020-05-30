#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/distributions.h"
#include "SWIFT/swift.h"
#include "SWIFT/c_m_calculators.h"

BOOST_AUTO_TEST_SUITE(CNCalculatorsUT)

///c_m_k = 2^(m/2) * int_{-1/2}^{1/2} (Re(char(2^{m+1} pi t) * exp(i 2 pi k t) dt
BOOST_AUTO_TEST_CASE(FastCoeffsTest)
{
    const double years = 0.1;
    GBM underlying_geometry(0.25, years);
    Swift::SwiftParameters params{
        4,
//        5,
        -8,
        164,
//        84,
        12,
        12
//        7
    };

    Swift::ParsevalCalculator parseval_calculator(underlying_geometry, params);
    Swift::ExplicitVietaCalculator explicit_vieta_calculator(underlying_geometry, params);
    Swift::FastVietaCalculator fast_vieta_calculator(underlying_geometry, params);

    double x = Distribution::GetXCompression(100, 80, 0.1, 0.0, years);

    std::vector<double> c_m_parseval = parseval_calculator.GetCMCoefs(x);
    std::vector<double> c_m_vieta = explicit_vieta_calculator.GetCMCoefs(x);
    std::vector<double> c_m_vieta_dft = fast_vieta_calculator.GetCMCoefs(x);

    for (int k = params.m_k1; k <= params.m_k2; ++k)
    {
        std::size_t index = k - params.m_k1;
        BOOST_CHECK(IsZero(c_m_parseval[index] - c_m_vieta[index], 1e-7));
        BOOST_CHECK(IsZero(c_m_vieta[index] - c_m_vieta_dft[index], 1e-7));
    }
}

BOOST_AUTO_TEST_SUITE_END()