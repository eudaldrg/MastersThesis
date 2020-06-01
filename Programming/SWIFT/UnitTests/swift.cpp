#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fstream>

#include "SWIFT/distributions.h"

BOOST_AUTO_TEST_SUITE(SwiftUT)

BOOST_AUTO_TEST_CASE(SIApproxTes)
{
    std::size_t j_bar = 14;

    std::vector<std::pair<double, double>> test_results = {
        {1, 0.5894898722360836351},
        {2, 0.4514116667901403133},
        {3, 0.5330932376182719825},
        {-1, -0.589489872236083635},
        {-2, -0.451411666790140313},
        {-3, -0.533093237618271982}
    };

    for (auto const& [x, exact] : test_results)
        BOOST_CHECK_CLOSE(Sign(x) * SIApprox(std::abs(x), j_bar), exact, MATH_EPSILON);
}

BOOST_AUTO_TEST_SUITE_END()