#pragma once

#include <boost/math/distributions/normal.hpp>

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