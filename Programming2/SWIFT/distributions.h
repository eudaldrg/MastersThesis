#pragma once

#include "my_math.h"

class Distribution
{
public:
    Distribution(double T) : m_T(T)
    { }
    virtual ~Distribution() = default;
    [[nodiscard]] virtual std::complex<double> GetNonXChar(double u) const = 0;
    [[nodiscard]] std::complex<double> GetChar(double u, double x) const
    {
        return GetNonXChar(u) * GetCharPositionFactor(x, u);
    }

    // The overall char is regular char * this.
    // x = log(F_{t, T}/K) = log(S_t * exp((r - q) * tau) / K) = (r - q) * tau * log (S_t / K)
    static std::complex<double> GetCharPositionFactor(double x, double u)
    {
        return std::exp(-1i * x * u);
    }

    static double GetXCompression(double future, double strike, double r, double q, double tau)
    {
//    return std::log(future * std::exp((r - q) * tau) / strike);
        return (r - q) * tau + std::log(future / strike);
    }

    double m_T;
};

// AKA B&S
class GBM : public Distribution
{
public:
	GBM(double vol, double T) : Distribution(T), m_vol(vol), m_vol_2(m_vol * m_vol) {}

//    [[nodiscard]] std::complex<double> GetChar(double u, double x) const
//    {
//        return std::exp(-u * m_T * (1i * (m_r - m_q - 0.5 * m_vol_2) + 0.5 * m_vol_2 * u)) * std::exp(-1i * x * u);
//    }
    [[nodiscard]] virtual std::complex<double> GetNonXChar(double u) const override final
    {
        return std::exp(-u * m_T * 0.5 * m_vol_2 * (u - 1i));
    }

	double m_vol;
	double m_vol_2;
};

/// Parameters of the heston distribution assuming
/// dS = mu S dt TODO:
/// dv =
struct HestonParameters
{
public:
	HestonParameters(double k, double v_bar, double sigma, double rho, double v_0) : m_k(k), m_nu_bar(v_bar), m_sigma(sigma), m_rho(rho), m_nu(v_0) {}
	double m_k;
	double m_nu_bar;
	double m_sigma;
	double m_rho;
	double m_nu;
};

class HestonDistribution : public Distribution
{
public:
	HestonDistribution(HestonParameters parameters, double T)
		: Distribution(T), m_parameters(parameters)
	{
		m_sigma_squared = m_parameters.m_sigma * m_parameters.m_sigma;
	}

	struct HelperVariables
	{
		std::complex<double> xi;
		std::complex<double> d;
		std::complex<double> A1;
		std::complex<double> A2;
		std::complex<double> A;
		std::complex<double> D;
	};

	HelperVariables GetHelperVariables(double u, double tau) const
	{
		auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
		std::complex<double> xi = k - sigma * rho * u * 1i;
		std::complex<double> d = std::sqrt(xi * xi + m_sigma_squared * (u * u + u * 1i));
		std::complex<double> A1 = (u * u + 1i * u) * std::sinh(0.5 * d * tau);
		std::complex<double> A2 = d * std::cosh(0.5 * d * tau) / v_0 + xi * std::sinh(d * tau * 0.5) / v_0;
		std::complex<double> A = A1 / A2;
		std::complex<double> D = std::log(d / v_0) + (k - d) * tau / 2.0 - std::log((d + xi) / (2 * v_0) + (d - xi) / (2 * v_0)
                                                                                                           * std::exp(-d * tau));

		return{xi, d, A1, A2, A, D};
	}

	[[nodiscard]] std::complex<double> GetNonXChar(double u) const override final
    {
        return GetCuiChar(u, m_T);
    }

	[[nodiscard]] std::complex<double> GetCuiChar(double u, double tau) const
    {
        auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
        auto const& [xi, d, A1, A2, A, D] = GetHelperVariables(u, tau);
	    return std::exp(- (k * v_bar * rho * tau * 1i * u) / sigma - A + (2 * k * v_bar * D) / (m_sigma_squared));
    }

	// The original paper used:
	// beta instead of xi
	// lambda instead of kappa
	// eta instead of sigma
	[[nodiscard]] std::complex<double> GetGatheralChar(double u, double tau)
    {
        auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
        std::complex<double> xi = k - sigma * rho * u * 1i;
        std::complex<double> d = std::sqrt(xi * xi + m_sigma_squared * (u * u + u * 1i));

        std::complex<double> r_minus = (xi - d) / m_sigma_squared;
        std::complex<double> r_plus = (xi + d) / m_sigma_squared;
//        std::complex<double> g1 = r_plus / r_minus;
        std::complex<double> g2 = r_minus / r_plus;

        std::complex<double> exp_minus_d_times_tau = std::exp(-d * tau);
        std::complex<double> G2 = (1.0 - exp_minus_d_times_tau) / (1.0 - g2 * exp_minus_d_times_tau);

        std::complex<double> D = r_minus * G2;
        std::complex<double> C = k * (r_minus * tau - 2.0 / m_sigma_squared * std::log((1.0 - g2 * exp_minus_d_times_tau) / (1.0 - g2)));

        return std::exp(C * v_bar + D * v_0);
    }

	[[nodiscard]] std::complex<double> GetHestonChar(double u, double tau)
	{
		auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
		std::complex<double> xi = k - sigma * rho * u * 1i;
		std::complex<double> d = std::sqrt(xi * xi + m_sigma_squared * (u * u + u * 1i));

        std::complex<double> r_minus = (xi - d) / m_sigma_squared;
        std::complex<double> r_plus = (xi + d) / m_sigma_squared;
        std::complex<double> g1 = r_plus / r_minus;
//        std::complex<double> g2 = r_minus / r_plus;

        std::complex<double> exp_d_times_tau = std::exp(d * tau);
        std::complex<double> G1 = (1.0 - exp_d_times_tau) / (1.0 - g1 * exp_d_times_tau);

        std::complex<double> D = r_plus * G1;
        std::complex<double> C = k * (r_plus * tau - 2.0 / m_sigma_squared * std::log((1.0 - g1 * exp_d_times_tau) / (1.0 - g1)));
		return std::exp(C * v_bar + D * v_0);
	}
	//[[nodiscard]] std::tuple<std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>>
	//	GetCharDvol0(double u, double x, double tau)
	//{
	//	auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
	//	auto const& [xi, d, A1, A2, A, D] = GetHelperVariables(u, tau, m_parameters);
	//	std::complex<double> char_val = GetChar(u, x, tau);
	//	auto dchar_dv_0 = -A / v_0 * char_val;
	//	auto dchar_dv_bar = (2.0 * k * D * char_val) / (sigma * sigma);
	//	/// Variables for dchar_drho
	//	auto dd_drho = -(xi * sigma * 1i * u) / d;
	//	auto dA2_drho = -(sigma * 1i * (2.0 + xi * tau) * (xi * std::cosh(0.5 * d * tau) + d * std::sinh(0.5 * d * tau)) / (2 * d * v_0);
	//	auto dB_drho = std::exp(k * tau * 0.5) / v_0 * (dd_drho / A2 - dA2_drho / (A2 * A2));
	//	auto dA1_drho = -(1i * u * (u * u + 1i * u) * tau * xi * sigma) / (2 * d) * std::cosh(0.5 * d * tau);
	//	auto dA_drho = dA1_drho / A2 - A / A2 * dA2_drho;
	//	auto dchar_drho = char_val * (-dA_drho + 2 * k * v_bar / (sigma * sigma * d) * (dd_drho - d / A2 * dA2_drho) - k * v_bar * tau * 1i * u / sigma);
	//	/// Variables for dchar_dk
	//	auto dchar_dk = ;
	//	/// Variable for dchar_dsigma
	//	auto dd_dsigma = (rho / sigma - 1.0 / xi) * dd_drho + sigma * u * u / d;

	//	auto dA1_dsigma = ;
	//	auto dA2_dsigma = ;
	//	auto dA_dsigma = ;
	//	auto dchar_dsigma = ;
	//}
private:
	HestonParameters m_parameters;
	double m_sigma_squared;
};