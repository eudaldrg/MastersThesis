#include <boost/test/unit_test.hpp>

#include "SWIFT/my_math.h"

#include <iostream>

namespace std {
ostream& operator<<(ostream& ss, vector<complex<double>> const& vec)
{
    for (auto const& complex : vec)
        ss << complex << " ";
    return ss;
}
}

BOOST_AUTO_TEST_SUITE(MyMathUT)

BOOST_AUTO_TEST_CASE(BitReverse)
{
    BOOST_CHECK_EQUAL(bitReverse(0b1101, 4), 0b1011);
    BOOST_CHECK_EQUAL(bitReverse(0b110101, 6), 0b101011);
    BOOST_CHECK_EQUAL(bitReverse(0b110101, 4), 0b1010);
}

//BOOST_AUTO_TEST_CASE(OwnFFTUT)
//{
//    std::cout << "Begin FFTUT" << std::endl;
//    std::vector<std::complex<double>> x = { {1,0}, {2, -1}, {0, -1}, {-1, 2} };
//    std::vector<std::complex<double>> X = { {2,0}, {-2, -2}, {0, -2}, {4, 4} };
//    std::vector<std::complex<double>> real_X = x;
//
//    Swift::fft(real_X, false);
//
//    BOOST_CHECK_EQUAL(X, real_X);
//
//    std::vector<std::complex<double>> real_x = real_X;
//    Swift::fft(real_x, true);
//    BOOST_CHECK_EQUAL(x, real_x);
//}

BOOST_AUTO_TEST_SUITE_END()