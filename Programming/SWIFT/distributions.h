#pragma once

#include "my_math.h"

class Distribution
{
public:
    explicit Distribution(double T) : m_T(T)
    { }
    virtual ~Distribution() = default;
    [[nodiscard]] virtual std::complex<double> GetNonXChar(double u) const = 0;
    [[nodiscard]] std::complex<double> GetChar(double u, double x) const
    {
        return GetNonXChar(u) * GetCharPositionFactor(x, u);
    }
    [[nodiscard]] virtual std::vector<std::complex<double>> GetNonXCharGradient(double u) const = 0;

    [[nodiscard]] std::vector<std::complex<double>> GetCharGradient(double u, double x) const
    {
        std::vector<std::complex<double>> result_components = GetNonXCharGradient(u);
        std::complex<double> x_factor = GetCharPositionFactor(x, u);
        for (auto& val : result_components)
            val *= x_factor;
        return result_components;
    }

    // The overall char is regular char * this.
    // x = log(F_{t, T}/K) = log(S_t * exp((r - q) * tau) / K) = (r - q) * tau * log (S_t / K)
    static std::complex<double> GetCharPositionFactor(double x, double u)
    {
        return std::exp(-1i * x * u);
    }

    virtual std::size_t GetNumberOfParameters() const = 0;

    static double GetXCompression(double future, double strike, double r, double q, double tau)
    {
//    return std::log(future * std::exp((r - q) * tau) / strike);
        return (r - q) * tau + std::log(future / strike);
    }

    double m_T;
};

struct GBMParameters
{
    double m_vol;
    static std::size_t GetNumberOfParameters()
    {
        return 1;
    }
};


// AKA B&S
class GBM : public Distribution
{
public:
	GBM(double vol, double T) : Distribution(T), m_parameters{vol}, m_vol_2(m_parameters.m_vol * m_parameters.m_vol)
	{ }

//    [[nodiscard]] std::complex<double> GetChar(double u, double x) const
//    {
//        return std::exp(-u * m_T * (1i * (m_r - m_q - 0.5 * m_vol_2) + 0.5 * m_vol_2 * u)) * std::exp(-1i * x * u);
//    }
    [[nodiscard]] std::complex<double> GetNonXChar(double u) const final
    {
        return std::exp(-u * m_T * 0.5 * m_vol_2 * (u - 1i));
    }

    [[nodiscard]] std::vector<std::complex<double>> GetNonXCharGradient(double /*u*/) const final
    {
	    throw std::runtime_error("Not implemented.");
    }
    [[nodiscard]] std::size_t GetNumberOfParameters() const final
    {
	    return GBMParameters::GetNumberOfParameters();
    }

    GBMParameters m_parameters;
    double m_vol_2;
};

/// Parameters of the heston distribution
struct HestonParameters
{
public:
    HestonParameters(double k, double v_bar, double sigma, double rho, double v_0) : m_k(k), m_v_bar(v_bar), m_sigma(sigma), m_rho(rho), m_v0(v_0) {}

    static std::size_t GetNumberOfParameters()
    {
        return 5;
    }

	double m_k;         // mean reversion rate
	double m_v_bar;     // long term variance
	double m_sigma;     // variance of volatility
	double m_rho;       // correlation between spot and volatility
	double m_v0;        // initial variance
};

std::ostream& operator<<(std::ostream& out, HestonParameters const& heston_parameters);

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

	[[nodiscard]] std::complex<double> GetNonXChar(double u) const final
    {
        return GetCuiChar(u, m_T);
    }

    [[nodiscard]] std::vector<std::complex<double>> GetNonXCharGradient(double u) const final
    {
        return GetCuiGradient(u, m_T);
    }

	[[nodiscard]] std::complex<double> GetCuiChar(double u, double tau) const
    {
        auto const& [k, v_bar, sigma, rho, v_0] = m_parameters;
        auto const& [xi, d, A1, A2, A, D] = GetHelperVariables(u, tau);
	    return std::exp(- (k * v_bar * rho * tau * 1i * u) / sigma - A + (2 * k * v_bar * D) / (m_sigma_squared));
    }

    [[nodiscard]] std::vector<std::complex<double>> GetCuiGradient(double u, double T) const
    {
        auto const& [kappa, v_bar, sigma, rho, v0] = m_parameters;

        std::complex<double> ui = 1i * u;
        std::complex<double> u_squared = u * u;

        // We need to evaluate everything at u1 = u - i and u2 = u.
        double sigma_times_rho = sigma * rho; // TAG: PRECOMPUTE
        std::complex<double> xi = kappa - sigma_times_rho * ui; // xi = kappa - sigma * rho * u * i
        std::complex<double> m = ui + u_squared; // m = u * i + u^2;
        std::complex<double> d = sqrt(pow(xi, 2) + m * m_sigma_squared); // d = sqrt(pow(xi,2) + m*pow(sigma,2));

        // g = exp(-kappa * b * rho * T * u1 * i / sigma);
        double kappa_v_bar_rho_T = kappa * v_bar * rho * T;  // TAG: PRECOMPUTE
        double minus_kappa_v_bar_rho_T_over_sigma = -kappa_v_bar_rho_T / sigma;  // TAG: PRECOMPUTE
        std::complex<double> g = exp(minus_kappa_v_bar_rho_T_over_sigma * ui);

        // alp, calp, salp
        double halft = 0.5 * T;
        std::complex<double> alpha = d * halft;
        std::complex<double> cosh_alpha = cosh(alpha);
        std::complex<double> sinh_alpha = sinh(alpha);
        std::complex<double> A2_times_v0 = d * cosh_alpha + xi * sinh_alpha;
        std::complex<double> A1 = m * sinh_alpha;
        std::complex<double> A_over_v0 = A1 / A2_times_v0;

        // B = d * exp(kappa * T / 2) / (A2 * v0);
        double exp_kappa_times_half_T = exp(kappa * halft); // exp(kappa * T / 2)
        std::complex<double> B = d * exp_kappa_times_half_T / A2_times_v0;

        double two_kappa_v_bar_over_sigma_squared = 2 * kappa * v_bar / m_sigma_squared;
        std::complex<double> D = log(d) + (kappa - d) * halft - log((d + xi) * 0.5 + (d - xi) * 0.5 * exp(-d * T));

        std::complex<double> H = xi * cosh_alpha + d * sinh_alpha;

        // lnB = log(B);
        std::complex<double> lnB = D;

        // partial b: y3 = y1*(2*kappa*lnB/pow(sigma,2)-kappa*rho*T*u1*i/sigma);
        double two_kappa_over_sigma_squared = two_kappa_v_bar_over_sigma_squared / v_bar;
        double minus_kappa_rho_T_over_sigma = minus_kappa_v_bar_rho_T_over_sigma / v_bar;

        std::complex<double> h_v_bar = two_kappa_over_sigma_squared * lnB + minus_kappa_rho_T_over_sigma * ui;

        // partial rho:
        double minus_kappa_v_bar_t_over_sigma = minus_kappa_v_bar_rho_T_over_sigma / rho; //-kappa * v_bar * T/sigma;

        std::complex<double> sigma_ui_over_d = sigma * ui / d;
        std::complex<double> pd_prho = -xi * sigma_ui_over_d;
        std::complex<double> pA1_prho = m * cosh_alpha * halft * pd_prho;
        std::complex<double> pA2_prho = -sigma_ui_over_d * H * (1.0 + xi * halft);
        std::complex<double> pA_prho = (pA1_prho - A_over_v0 * pA2_prho) / A2_times_v0;
        std::complex<double> pd_phrho_minus_pA2_prho_times_d_over_A2 = pd_prho - pA2_prho * d / A2_times_v0;
        std::complex<double> pB_prho = exp_kappa_times_half_T / A2_times_v0 * pd_phrho_minus_pA2_prho_times_d_over_A2;
        std::complex<double> h_rho = -v0 * pA_prho + two_kappa_v_bar_over_sigma_squared * pd_phrho_minus_pA2_prho_times_d_over_A2 / d + minus_kappa_v_bar_t_over_sigma * ui;

        // partial kappa:
        double v_bar_rho_T_over_sigma = v_bar * rho * T / sigma;
        double two_v_bar_over_sigma_squared = two_kappa_v_bar_over_sigma_squared / kappa; // 2 * v_bar / sigma_squared;

        std::complex<double> minus_one_over_sigma_ui = -1.0 / (sigma * ui);
        std::complex<double> pB_pa = minus_one_over_sigma_ui * pB_prho + B * halft;
        std::complex<double> h_kappa = -v0 * pA_prho * minus_one_over_sigma_ui + two_v_bar_over_sigma_squared * lnB + kappa * two_v_bar_over_sigma_squared * pB_pa / B - v_bar_rho_T_over_sigma * ui;

        // partial sigma:
        double rho_over_sigma = rho / sigma;
        double four_kappa_v_bar_over_sigma_cubed = 4 * kappa * v_bar / pow(sigma, 3);
        double kappa_v_bar_rho_T_over_sigma_squared = kappa_v_bar_rho_T / m_sigma_squared;
        std::complex<double> pd_pc = (rho_over_sigma - 1.0 / xi) * pd_prho + sigma * u_squared / d;
        std::complex<double> pA1_pc = m * cosh_alpha * halft * pd_pc;
        std::complex<double> pA2_pc = rho_over_sigma * pA2_prho - 1.0 / ui * (2.0 / (T * xi) + 1.0) * pA1_prho + sigma * halft * A1;
        std::complex<double> pA_pc = pA1_pc / A2_times_v0 - A_over_v0 / A2_times_v0 * pA2_pc;
        std::complex<double> h_sigma = -v0 * pA_pc - four_kappa_v_bar_over_sigma_cubed * lnB + two_kappa_v_bar_over_sigma_squared / d * (pd_pc - d / A2_times_v0 * pA2_pc) + kappa_v_bar_rho_T_over_sigma_squared * ui;

        // F = S * e^((r - q) * T)
        // characteristic function: y1 = exp(i * log(F / S) * u) * exp(-A + 2 * kappa * b / pow(sigma, 2) * D) * g
        // But we only care about the second exponential, the rest depends only on market parameters and will be computed separately.
        std::complex<double> char_u = exp(-v0 * A_over_v0 + two_kappa_v_bar_over_sigma_squared * D) * g;

        return {char_u * h_kappa, char_u * h_v_bar, char_u * h_sigma, char_u * h_rho, char_u * A_over_v0};
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


    [[nodiscard]] std::size_t GetNumberOfParameters() const final
    {
        return HestonParameters::GetNumberOfParameters();
    }
private:
	double m_sigma_squared;
	HestonParameters m_parameters;
};