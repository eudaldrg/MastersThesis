#include "c_m_calculators.h"

#include "SWIFT/distributions.h"
#include "FFTW3/include_fftw3.h"
#include "SWIFT/my_math.h"
#include "SWIFT/swift.h"

namespace Swift {

CMCalculator::CMCalculator(Distribution const& distribution, SwiftParameters const& params) : m_distribution(distribution), m_params(params)
{ }

std::vector<double> ParsevalCalculator::GetCMCoefs(double x) const
{
    std::vector<double> c_m;
    c_m.reserve(m_params.m_k2 - m_params.m_k1);
    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
    {
        auto f = [this, x, k](double t)
        {
            double u = std::pow(2.0, m_params.m_m + 1) * MY_PI * t;
            return (m_distribution.GetChar(u, x) * std::exp(1i * 2.0 * MY_PI * static_cast<double>(k) * t)).real();
        };
        c_m.push_back(TrapezoidInt(-0.5, 0.5, f, m_integral_buckets) * std::pow(2.0, m_params.m_m / 2.0));
    }
    return c_m;
}

ParsevalCalculator::ParsevalCalculator(Distribution const& distribution, SwiftParameters const& params, std::size_t integral_buckets)
: CMCalculator(distribution, params), m_integral_buckets(integral_buckets)
{ }

std::vector<double> ExplicitVietaCalculator::GetCMCoefs(double x) const
{
    const std::size_t N = 1UL << m_params.m_iota_dens_coefs;
    const std::size_t two_to_the_m = 1UL << m_params.m_m;
    const std::size_t two_to_the_j_minus_1 = std::pow(2.0, (m_params.m_iota_dens_coefs - 1));

    std::vector<double> c_m;
    c_m.reserve(m_params.m_k2 - m_params.m_k1);
    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
    {
        double c_m_k = 0;
        for (std::size_t i = 1; i <= two_to_the_j_minus_1; ++i)
        {
            double u = (2 * i - 1) * MY_PI * two_to_the_m / N;

            std::complex<double> f_hat = m_distribution.GetChar(u, x);
            std::complex<double> cv = f_hat * std::exp(1i * static_cast<double>(k) * MY_PI * (2.0 * i - 1.0) / static_cast<double>(N));
            c_m_k += cv.real();
        }
        c_m.push_back(c_m_k * std::sqrt(two_to_the_m) / two_to_the_j_minus_1);
    }
    return c_m;
}

ExplicitVietaCalculator::ExplicitVietaCalculator(Distribution const& distribution, SwiftParameters const& params) : CMCalculator(distribution, params)
{ }

std::vector<double> FastVietaCalculator::GetCMCoefs(double x) const
{
    std::size_t N = 1UL << m_params.m_iota_dens_coefs;
    std::size_t two_to_the_m = 1UL << m_params.m_m;
    const std::size_t two_to_the_j_minus_1 = std::pow(2.0, (m_params.m_iota_dens_coefs - 1));

    std::vector<double> c_m_dft;
    std::vector<double> c_m;
    std::vector<std::complex<double>> frequencies;

    for (std::size_t i = 0; i < N; ++i)
    {
        if (i < N / 2)
            frequencies.push_back(m_distribution.GetChar(static_cast<double>(2 * i + 1) * MY_PI * two_to_the_m / N, x));
        else
            frequencies.emplace_back(0, 0);
    }
    std::vector<std::complex<double>> times = MY_IDFT(frequencies);

    for (int k = m_params.m_k1; k <= m_params.m_k2; ++k)
    {
        std::size_t k_mod_N = Mod(k, static_cast<int>(N));
        std::complex<double> idf_part = times[k_mod_N];
        std::complex<double> c_m_k_complex_part_dft = std::exp(1i * static_cast<double>(k) * MY_PI / static_cast<double>(N)) * idf_part;
        c_m.push_back(c_m_k_complex_part_dft.real() * std::sqrt(two_to_the_m) / two_to_the_j_minus_1);
    }
    return c_m;
}

std::vector<std::complex<double>> FastVietaCalculator::MY_IDFT(std::vector<std::complex<double>> const& X, bool fast)
{
    if (!fast)
        return IDFT(X, false);

    fftw_complex* frequency_values, * time_values;
    fftw_plan p;
    frequency_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * X.size());
    time_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * X.size());
    p = fftw_plan_dft_1d(X.size(), frequency_values, time_values, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (std::size_t i = 0; i < X.size(); ++i)
    {
        frequency_values[i][0] = X[i].real();
        frequency_values[i][1] = X[i].imag();
    }

    fftw_execute(p);
    std::vector<std::complex<double>> time_values_std;
    for (std::size_t i = 0; i < X.size(); ++i)
        time_values_std.emplace_back(time_values[i][0], time_values[i][1]);
    return time_values_std;
}

FastVietaCalculator::FastVietaCalculator(Distribution const& distribution, SwiftParameters const& params) : CMCalculator(distribution, params)
{ }

}