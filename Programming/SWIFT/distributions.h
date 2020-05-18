#pragma once

#include <boost/math/distributions/normal.hpp>

#include "my_math.h"
#include "swift.h"

using namespace std::complex_literals;

inline double GetStandardCDF(double x)
{
	return boost::math::cdf(boost::math::normal(0, 1), x);
}
//inline double GetStandardPDF(double x)
//{
//	return boost::math::pdf(boost::math::normal(0, 1), x);
//}

inline double GetCashOrNothingCallPayOff(std::size_t m, int k, std::size_t J)
{
	const int sign = (k == 0) ? 0 : ((k > 0) ? 1 : -1);
	return std::pow(2.0, m / -2.0) * (sign * Swift::SIApprox(std::abs(k), J) + 0.5);
}

inline double GetC_j(std::size_t j, std::size_t J)
{
	return (2 * j - 1) * Swift::PI / (1 << J);
}

inline double GetI1(double a, double b, std::size_t j, std::size_t J, std::size_t m, int k)
{
	double two_to_the_m = 1 << m;
	double C_j = GetC_j(j, J);
	double C_j_times_two_to_the_m = C_j * two_to_the_m;
	return (C_j_times_two_to_the_m) / (1 + C_j_times_two_to_the_m * C_j_times_two_to_the_m)
		* (std::exp(b) * std::sin(C_j * (two_to_the_m * b - k)) - std::exp(a) * std::sin(C_j * (two_to_the_m * a - k))
			+ (std::exp(b) * std::cos(C_j * (two_to_the_m * b - k)) - std::exp(a) * std::cos(C_j * (two_to_the_m * a - k))) / C_j_times_two_to_the_m);
}

inline double GetI2(double a, double b, std::size_t j, std::size_t J, std::size_t m, int k)
{
	double two_to_the_m = 1 << m;
	double C_j = GetC_j(j, J);
	double C_j_times_two_to_the_m = C_j * two_to_the_m;
	double return_value = (std::sin(C_j * (two_to_the_m * b - k)) - std::sin(C_j * (two_to_the_m * a - k))) / C_j_times_two_to_the_m;
	assert(std::isnormal(return_value) || return_value == 0.0);
	return return_value;
}

inline double GetEuropeanPayoff(double /*strike_price*/ K, std::size_t m, int k, int k1, int k2, std::size_t J, bool is_call)
{
	double return_value = 0.0;
	if ((is_call && k2 <= 0) || (!is_call && k1 >= 0))
		return return_value;
	int k1_bar = std::max(k1, 0);
	int k2_bar = std::min(k2, 0);
	int two_to_the_j_minus_one = 1 << (J - 1);
	double two_to_the_m = 1 << m;
	for (std::size_t j = 1; j <= two_to_the_j_minus_one; ++j)
	{
		double I1 = is_call ? GetI1(k1_bar / two_to_the_m, k2 / two_to_the_m, j, J, m, k) : GetI1(k1 / two_to_the_m, k2_bar / two_to_the_m, j, J, m, k);
		double I2 = is_call ? GetI2(k1_bar / two_to_the_m, k2 / two_to_the_m, j, J, m, k) : GetI2(k1 / two_to_the_m, k2_bar / two_to_the_m, j, J, m, k);
		assert(std::isnormal(I1) || I1 == 0.0);
		assert(std::isnormal(I2) || I2 == 0.0);
		return_value += I1 - I2;
	}
	return_value *= K * std::sqrt(two_to_the_m) / two_to_the_j_minus_one;
	if (!is_call)
		return_value *= -1.0;
	return return_value;
}

inline double GetCashOrNothingPrice(double K, double F, double r, double years, double vol, double q = 0)
{
	double vol_times_sqrt_years = vol * sqrt(years);
	double d1 = (std::log(F / K) + (r - q + 0.5 * vol * vol) * years) / vol_times_sqrt_years;
	double d2 = d1 - vol_times_sqrt_years;
	return std::exp(-r * years) * GetStandardCDF(d2);
}


class GBM
{
public:
	GBM(double vol, double T, double r, double q = 0.0) : m_vol(vol), m_vol_2(m_vol * m_vol), m_T(T), m_r(r), m_q(q) {}

	[[nodiscard]] std::complex<double> GetChar(double w, double x) const
	{
		return std::exp(-w * m_T * (1i * (m_r - m_q - 0.5 * m_vol_2) + 0.5 * m_vol_2 * w)) * std::exp(-1i * x * w);
	}

	double m_vol;
	double m_vol_2;
	double m_T;
	double m_r;
	double m_q;
};

inline double GetBSPrice(double K, double F, double inverse_exp_times_interest, double years, double vol, bool is_call)
{
	double vol_times_sqrt_years = vol * sqrt(years);
	double d1 = std::log(F / K) / vol_times_sqrt_years + 0.5 * vol_times_sqrt_years;
	double d2 = d1 - vol_times_sqrt_years;
	return is_call ? F * GetStandardCDF(d1) - inverse_exp_times_interest * K * GetStandardCDF(d2)
	: GetStandardCDF(-d2) * K * inverse_exp_times_interest - GetStandardCDF(-d1) * F;
}


struct HestonParameters
{
public:
	HestonParameters(double k, double v_bar, double sigma, double rho, double v_0) : m_k(k), m_v_bar(v_bar), m_sigma(sigma), m_rho(rho), m_v_0(v_0) {}
	double m_k;
	double m_v_bar;
	double m_sigma;
	double m_rho;
	double m_v_0;
};

class HestonDistribution 
{
public:
	HestonDistribution(HestonParameters parameters, double r, double q = 0.0)
		: m_parameters(parameters), m_r(r), m_q(q)
	{
		m_sigma_squared = m_parameters.m_sigma * m_parameters.m_sigma;
	}

	struct HelperVariables
	{
		std::complex<double> e;
		std::complex<double> d;
		std::complex<double> A1;
		std::complex<double> A2;
		std::complex<double> A;
		std::complex<double> D;
	};

	HelperVariables GetHelperVariables(double u, double tau, HestonParameters const& heston_parameters)
	{
		auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
		std::complex<double> e = k - sigma * rho * u * 1i;
		std::complex<double> d = std::sqrt(e * e + m_sigma_squared * (u * u + u * 1i));
		std::complex<double> A1 = (u * u + 1i * u) * std::sinh(0.5 * d * tau);
		std::complex<double> A2 = d * std::cosh(0.5 * d * tau) / v_0 + e * std::sinh(d * tau * 0.5) / v_0;
		std::complex<double> A = A1 / A2;
		std::complex<double> D = std::log(d / v_0) + (k - d) * tau / 2.0 - std::log((d + e) / (2 * v_0) + (d - e) / (2 * v_0)
			* std::exp(-d * tau));

		HelperVariables hv{e, d, A1, A2, A, D};
	}

	[[nodiscard]] std::complex<double> GetChar(double u /*w in GBM*/, double x, double tau /*time to expiry*/)
	{
		//double forward_price = S_0 * std::exp((r - q) * tau);
		double forw_over_s0 = std::exp((m_r - m_q) * tau);
		auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
		auto const& [e, d, A1, A2, A, D] = GetHelperVariables(u, tau, m_parameters);
		return std::exp(1i * u * std::log(forw_over_s0) - (k * v_bar * rho * tau * 1i * u) / sigma - A + (2 * k * v_bar * D)
			/ (m_sigma_squared)) * std::exp(-1i * x * u);
	}
	[[nodiscard]] std::complex<double> GetHestonChar(double u, double tau)
	{
		double forw_over_s0 = std::exp((m_r - m_q) * tau);
		auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
		std::complex<double> e = k - sigma * rho * u * 1i;
		std::complex<double> d = std::sqrt(e * e + m_sigma_squared * (u * u + u * 1i));
		std::complex<double> g1 = (e + d) / (e - d);
		return std::exp(1i * u * std::log(forw_over_s0)
			+ k * v_bar / (sigma * sigma) * ((e + d) * tau - std::complex<double>(2) * std::log((std::complex<double>(1) - g1 * std::exp(d * tau)) / (std::complex<double>(1) - g1)))
			+ v_0 / (sigma * sigma) * (e + d) * ((std::complex<double>(1) - std::exp(d * tau)) / (std::complex<double>(1) - g1 * std::exp(d * tau))));
	}
	//[[nodiscard]] std::tuple<std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>> 
	//	GetCharDvol0(double u, double x, double tau)
	//{
	//	auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
	//	auto const& [e, d, A1, A2, A, D] = GetHelperVariables(u, tau, m_parameters);
	//	std::complex<double> char_val = GetChar(u, x, tau);
	//	auto dchar_dv_0 = -A / v_0 * char_val;
	//	auto dchar_dv_bar = (2.0 * k * D * char_val) / (sigma * sigma);
	//	/// Variables for dchar_drho
	//	auto dd_drho = -(e * sigma * 1i * u) / d;
	//	auto dA2_drho = -(sigma * 1i * (2.0 + e * tau) * (e * std::cosh(0.5 * d * tau) + d * std::sinh(0.5 * d * tau)) / (2 * d * v_0);
	//	auto dB_drho = std::exp(k * tau * 0.5) / v_0 * (dd_drho / A2 - dA2_drho / (A2 * A2));
	//	auto dA1_drho = -(1i * u * (u * u + 1i * u) * tau * e * sigma) / (2 * d) * std::cosh(0.5 * d * tau);
	//	auto dA_drho = dA1_drho / A2 - A / A2 * dA2_drho;
	//	auto dchar_drho = char_val * (-dA_drho + 2 * k * v_bar / (sigma * sigma * d) * (dd_drho - d / A2 * dA2_drho) - k * v_bar * tau * 1i * u / sigma);
	//	/// Variables for dchar_dk
	//	auto dchar_dk = ;
	//	/// Variable for dchar_dsigma
	//	auto dd_dsigma = (rho / sigma - 1.0 / e) * dd_drho + sigma * u * u / d;

	//	auto dA1_dsigma = ;
	//	auto dA2_dsigma = ;
	//	auto dA_dsigma = ;
	//	auto dchar_dsigma = ;
	//}
private:
	HestonParameters m_parameters;
	double m_r; //< Interest rate
	double m_q; //< dividents
	double m_sigma_squared;
};