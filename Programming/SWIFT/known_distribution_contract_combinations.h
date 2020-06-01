#pragma once

#include <boost/math/distributions/normal.hpp>
#include "distributions.h"

inline double GetStandardCDF(double x)
{
    return boost::math::cdf(boost::math::normal(0, 1), x);
}
//inline double GetStandardPDF(double x)
//{
//	return boost::math::pdf(boost::math::normal(0, 1), x);
//}

/// https://quant.stackexchange.com/questions/15489/derivation-of-the-formulas-for-the-values-of-european-asset-or-nothing-and-cash
inline double GetBSCashOrNothingPrice(double K, double F, double r, double years, double vol, double q = 0)
{
    double vol_times_sqrt_years = vol * sqrt(years);
    double d1 = (std::log(F / K) + (r - q + 0.5 * vol * vol) * years) / vol_times_sqrt_years;
    double d2 = d1 - vol_times_sqrt_years;
    return std::exp(-r * years) * GetStandardCDF(d2);
}

inline double GetBSEuropeanPrice(double K, double F, double r, double years, double vol, bool is_call, double q = 0)
{
    double inverse_exp_t_r = std::exp(-r * years);
    double inverse_exp_t_q = std::exp(-q * years);
    double vol_times_sqrt_years = vol * sqrt(years);
    double d1 = (std::log(F / K) + (r - q + 0.5 * vol * vol) * years) / vol_times_sqrt_years;
    double d2 = d1 - vol_times_sqrt_years;
    return is_call ? F * inverse_exp_t_q * GetStandardCDF(d1) - inverse_exp_t_r * K * GetStandardCDF(d2)
                   : K * inverse_exp_t_r * GetStandardCDF(-d2) - F * inverse_exp_t_q * GetStandardCDF(-d1);
}

inline double GetHestonEuropeanPrice(HestonParameters const& parameters, double S, double K, double r, double T, double q)
{
    double inverse_exponential_r_times_T = std::exp(-r * T);
    double inverse_exponential_q_times_T = std::exp(-q * T);
    double K_over_S = K / S;

    // retrieve model parameters
    double sigma_squared = parameters.m_sigma * parameters.m_sigma;

    // Heston Pricing Formula: C = S * e^(-qT) * P1 - K * e^(-rT) * P2
    // P1 = 0.5 + 1/pi * Int_0^inf(Re(e^(-i * u * log (K / S_0)) / iu * Char(u - i) / Char(i)) = 0.5 + Qv1 / (S_t * pi)
    // P2 = 0.5 + 1/pi * Int_0^inf(Re(e^(-i * u * log (K / S_0)) / iu * Char(u)) = 0.5 + Qv2 / pi

    // If we compute S * C / S, we can avoid some computations in each of the quadratures by pre-computing K / S
    // Heston Pricing Final Formula: C = S * 0.5 * (e^(-qT) - K / S * e^(-rT)) + e^(-rT) / PI *
    //  Int_0^inf
    //  [
    //    Re( e^(-iu * log K / S) / iu (Char(u - i) -  K / S * Char(u - i))) du)
    //  ]

    auto F = [K_over_S, r, T, &parameters, sigma_squared](double u)
    {
        auto const& [kappa, v_bar, sigma, rho, v0] = parameters;
        // We need to evaluate the Char function at (u - i) and at u and do
        std::complex<double> ui = 1i * u;
        std::complex<double> u_minus_i_times_i = 1i * (u - 1i);
        // xi = kappa-1i*sigma*rho*u1;
        double sigma_rho = sigma * rho;
        std::complex<double> xi_M1 = kappa - sigma_rho * u_minus_i_times_i;
        std::complex<double> xi_M2 = xi_M1 + sigma_rho;

        // m = 1i*u1 + pow(u1,2);
        std::complex<double> m_M1 = ui + 1.0 + pow(u - 1i, 2); // m_M1 = (PQ_M-1i)*1i + pow(PQ_M-1i, 2);
        std::complex<double> m_M2 = ui + pow(u - 0.0 * 1i, 2);

        // d = sqrt(pow(kes,2) + m*pow(sigma,2));
        std::complex<double> d_M1 = sqrt(pow(xi_M1, 2) + m_M1 * sigma_squared);
        std::complex<double> d_M2 = sqrt(pow(xi_M2, 2) + m_M2 * sigma_squared);

        // g = exp(-kappa*v_bar*rho*T*u1*1i/sigma);
        double tmp1 = -kappa * v_bar * rho * T / sigma;
        double g = exp(tmp1);
        std::complex<double> g_M2 = exp(tmp1 * ui);
        std::complex<double> g_M1 = g_M2*g;

        // alp, calp, salp
        double half_years_to_expiry = 0.5 * T;
        std::complex<double> alpha = d_M1*half_years_to_expiry;
        std::complex<double> calp_M1 = cosh(alpha);
        std::complex<double> salp_M1 = sinh(alpha);

        alpha = d_M2*half_years_to_expiry;
        std::complex<double> calp_M2 = cosh(alpha);
        std::complex<double> salp_M2 = sinh(alpha);

        // A2 = d*calp + kes*salp;
        std::complex<double> A2_M1 = d_M1*calp_M1 + xi_M1 * salp_M1;
        std::complex<double> A2_M2 = d_M2*calp_M2 + xi_M2 * salp_M2;

        // A1 = m*salp;
        std::complex<double> A1_M1 = m_M1*salp_M1;
        std::complex<double> A1_M2 = m_M2*salp_M2;

        // A = A1/A2;
        std::complex<double> A_M1 = A1_M1/A2_M1;
        std::complex<double> A_M2 = A1_M2/A2_M2;

        double twice_kappa_v_bar_over_sigma_squared = 2 * kappa * v_bar / sigma_squared;
        std::complex<double> D_M1 = log(d_M1) + (kappa - d_M1) * half_years_to_expiry - log((d_M1 + xi_M1) * 0.5 + (d_M1 - xi_M1) * 0.5 * exp(-d_M1 * T));
        std::complex<double> D_M2 = log(d_M2) + (kappa - d_M2) * half_years_to_expiry - log((d_M2 + xi_M2) * 0.5 + (d_M2 - xi_M2) * 0.5 * exp(-d_M2 * T));

        // F = S * e^(r * T);
        double log_F_over_S = r * T;
        // characteristic function: Char(u) = exp(ui * log(F / S)) * exp(-v0*A) * g * exp(2*kappa*v_bar/pow(sigma,2)*D) = exp(ui * log(F / S) - v0 * A + 2 * kappa * v_bar / sigma^2 * D) * g
        std::complex<double> char_u_minus_i = exp(u_minus_i_times_i * log_F_over_S - v0 * A_M1 + twice_kappa_v_bar_over_sigma_squared * D_M1) * g_M1;
        std::complex<double> char_u = exp(ui * log_F_over_S - v0 * A_M2 + twice_kappa_v_bar_over_sigma_squared * D_M2) * g_M2;

        // h = e^(-iu log(K / S)) / ui = (K / S) ^ (-ui) / ui
        std::complex<double> h = std::pow(K_over_S, -ui) / ui;
        return std::real(h * (char_u_minus_i - K_over_S * char_u));
    };

    GaussLegendreIntegrator gauss_legendre_integrator(GaussLegendreIntegrator::GaussLegendreMode::Fast);
    double quadrature_value = gauss_legendre_integrator.GetIntegral(F, 0, 200);
    return S * 0.5 * (inverse_exponential_q_times_T - K_over_S * inverse_exponential_r_times_T) + inverse_exponential_r_times_T / MY_PI * quadrature_value;
}
