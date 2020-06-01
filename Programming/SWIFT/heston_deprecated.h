pragma once

#include <complex>
#include <cmath>

namespace Heston 
{
/// Model parameters
/// k: 
/// v_bar: 
/// sigma: 
/// rho: 
/// v_0: 

/// Market parameters
///
/// S_t: Spot price
/// S_0: future price
/// K: strike
/// tau: time to expiry
/// q: dividents

using namespace std::complex_literals;

/// u is the frequency domain variable associated to what K means in time domain.
std::complex<double> EvaluateCharacteristicFunction(HestonParameters model_parameters, std::complex<double> u, double tau,
	double S_t, double K, double r, double q)
{
	double F = S_t * std::exp((r - q) * tau);
	auto const& [k, v_bar, sigma, rho, v_0] = model_parameters;
	std::complex<double> e = k - sigma * u * 1i;
	std::complex<double> d = std::sqrt(e * e + sigma * sigma * (u * u + u * 1i));
	std::complex<double> A1 = (u * u + 1i * u) * std::sinh(0.5 * d * tau);
	std::complex<double> A2 = d * std::cosh(0.5 * d * tau) / v_0 + e * std::sinh(d * tau * 0.5) / v_0;
	std::complex<double> A = A1 / A2;
	std::complex<double> D = std::log(d / v_0) + (k - d) * tau / 2.0 - std::log((d + e) / (2 * v_0) + (d - e) / (2 * v_0) * std::exp(-d * tau));
	return std::exp(1i * u * std::log(F / S_0) - (k * v_bar * rho * tau * 1i * u) / sigma - A + (2 * k * v_bar * D) / (sigma * sigma));
}
}