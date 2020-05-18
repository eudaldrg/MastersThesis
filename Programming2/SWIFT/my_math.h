#pragma once

#include <boost/math/constants/constants.hpp>
#include <complex>
#include <cmath>
#include <iterator>

namespace Swift {
using namespace std::complex_literals;
const double PI = boost::math::constants::pi<double>();
const double MATH_EPSILON = 1e-5;

/// \param [in] x: positive integer to reverse.
/// \param [in] log2n: size of x in bits.
/// \return bitwise reverse of x. i.e. x = 0b1101 and log2n = 4 returns n = 0b1011
inline unsigned int bitReverse(std::size_t x, std::size_t log2n) {
    std::size_t n = 0;
    while (log2n > 0)
    {
        n <<= 1U;
        n |= (x & 1U);
        x >>= 1U;
        --log2n;
    }
    return n;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline bool IsZero(double val, double tol = MATH_EPSILON)
{
    return std::abs(val) <= tol;
}

inline bool IsSame(double lhs, double rhs, double abs_tol = MATH_EPSILON, double rel_tol = MATH_EPSILON)
{
    if (std::abs(lhs - rhs) > abs_tol)
        return false;
    if (!IsZero(lhs, abs_tol))
        return std::abs((lhs - rhs) / lhs) <= rel_tol;
    return std::abs((lhs - rhs) / rhs) <= rel_tol;
}

inline std::vector<std::complex<double>> DFT(std::vector<std::complex<double>> const& x)
{
    std::size_t N = x.size();
    std::vector<std::complex<double>> X;
    for (std::size_t k = 0; k < N; ++k)
    {
        std::complex<double> current_transform{0.0,0.0};
        for (std::size_t n = 0; n < N; ++n)
            current_transform += x[n] * std::exp(-1i * 2.0 * PI / static_cast<double>(N) * static_cast<double>(k) * static_cast<double>(n));
        X.push_back(current_transform);
    }
    return X;
}

inline std::vector<std::complex<double>> IDFT(std::vector<std::complex<double>> const& X)
{
    std::size_t N = X.size();
    std::vector<std::complex<double>> x;
    for (std::size_t n = 0; n < N; ++n)
    {
        std::complex<double> current_transform{0.0,0.0};
        for (std::size_t k = 0; k < N; ++k)
            current_transform += X[k] * std::exp(1i * 2.0 * PI / static_cast<double>(N) * static_cast<double>(k) * static_cast<double>(n));
        x.push_back(current_transform / static_cast<double>(N));
    }
    return x;
}

inline void fft(std::vector<std::complex<double>>& a, bool invert) {
    int n = a.size();
    if (n == 1)
        return;

    std::vector<std::complex<double>> a0(n / 2);
    std::vector<std::complex<double>> a1(n / 2);
    for (int i = 0; 2 * i < n; i++) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }
    fft(a0, invert);
    fft(a1, invert);

    double ang = 2 * PI / n * (invert ? -1 : 1);
    std::complex<double> w(1), wn(cos(ang), sin(ang));
    for (int i = 0; 2 * i < n; i++) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        if (invert) {
            a[i] /= 2;
            a[i + n / 2] /= 2;
        }
        w *= wn;
    }
}

inline std::vector<std::complex<double>> FFT(std::vector<std::complex<double>> a, std::size_t log2n)
{
    std::size_t n = 1UL << log2n;
    std::vector<std::complex<double>> b(n);
    for (std::size_t i = 0; i < n; ++i)
        b[bitReverse(i, log2n)] = a[i];

    for (std::size_t s = 1; s <= log2n; ++s) {
        std::size_t m = 1U << s;
        std::size_t m2 = m >> 1U;
        std::complex<double> w(1, 0);
        std::complex<double> wm = exp(-PI * 1i / static_cast<double>(m2));
        for (std::size_t j = 0; j < m2; ++j) {
            for (std::size_t k = j; k < n; k += m) {
                std::complex<double> t = w * b[k + m2];
                std::complex<double> u = b[k];
                b[k] = u + t;
                b[k + m2] = u - t;
            }
            w *= wm;
        }
    }
    return b;
}
}