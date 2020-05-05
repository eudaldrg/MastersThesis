#include <iostream>

//#include <fftw3.h>
#include "fftw3.h"
#include "distributions.h"
#include "swift.h"

//C:\Program Files (x86)\mingw-w64\i686-8.1.0-posix-dwarf-rt_v6-rev0\mingw32\bin>g++ D:\GoogleDrive\TFM\Exec\main.cpp -o D:\GoogleDrive\TFM\Exec\a.exe -I "C:\Program Files (x86)\fftw" -I "C:\Program Files\boost\boost_1_71_0" -L "C:\Program Files (x86)\fftw" -l fftw3-3

void TestSI()
{
	std::size_t j_bar = 5;
	
	std::cout << "TEST SI APPROX WITH JBAR EQUAL " << j_bar << std::endl;
	std::cout << "x = 1. Approx " << Swift::SIApprox(1, j_bar) << " Expected " << 0.58949 << std::endl;
	std::cout << "x = 2. Approx " << Swift::SIApprox(2, j_bar) << " Expected " << 0.45141<< std::endl;
	std::cout << "x = 3. Approx " << Swift::SIApprox(3, j_bar) << " Expected " << 0.53309 << std::endl;
}

void TestFFTW3()
{
	std::cout << "Begin TESTFFTW3" << std::endl;
	fftw_complex* in, * out;
	fftw_plan p;
	int N = 4;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	std::vector<std::complex<double>> x = { {1,0}, {2, -1}, {0, -1}, {-1, 2} };
	for (std::size_t i = 0; i < x.size(); ++i)
	{
		in[i][0] = x[i].real();
		in[i][1] = x[i].imag();
	}
	fftw_execute(p); /* repeat as needed */
	std::vector<std::complex<double>> X = { {2,0}, {-2, -2}, {0, -2}, {4, 4} };
	std::cout << "Expected result";
	for (std::complex<double> val : X)
		std::cout << " " << val;
	std::cout << std::endl;
	std::cout << "Obtained result";
	for (std::size_t i = 0; i < N; ++i)
		std::cout << " (" << out[i][0] << "," << out[i][1] << ")";
	std::cout << std::endl;
	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
}

void BitReverseUT()
{
	std::cout << "Begin BITREVERSEUT" << std::endl;
	{
		std::size_t x = 0b1101;
		std::size_t log2n = 4;
		std::size_t n = 0b1011;
		std::cout << "Reverse of x: " << Swift::bitReverse(x, log2n) << " expected " << n << std::endl;
	}
	{
		std::size_t x = 0b110101;
		std::size_t log2n = 6;
		std::size_t n = 0b101011;
		std::cout << "Reverse of x: " << Swift::bitReverse(x, log2n) << " expected " << n << std::endl;
	}
	{
		std::size_t x = 0b110101;
		std::size_t log2n = 4;
		std::size_t n = 0b1010;
		std::cout << "Reverse of x: " << Swift::bitReverse(x, log2n) << " expected " << n << std::endl;
	}
}

void FFTUT()
{
	std::cout << "Begin FFTUT" << std::endl;
	std::vector<std::complex<double>> x = { {1,0}, {2, -1}, {0, -1}, {-1, 2} };
	std::vector<std::complex<double>> X = { {2,0}, {-2, -2}, {0, -2}, {4, 4} };
	std::vector<std::complex<double>> real_X = x;
	Swift::fft(real_X, false);
	std::cout << "Expected result";
	for (std::complex<double> val : X)
		std::cout << " " << val;
	std::cout << std::endl;
	std::cout << "Obtained result";
	for (std::complex<double> val : real_X)
		std::cout << " " << val;
	std::cout << std::endl;

	std::cout << "Original vec";
	for (std::complex<double> val : x)
		std::cout << " " << val;
	std::cout << std::endl;
	std::vector<std::complex<double>> real_x = real_X;
	Swift::fft(real_x, true);
	std::cout << "Recovered vec";
	for (std::complex<double> val : real_x)
		std::cout << " " << val;
	std::cout << std::endl;
}

void CashOrNothingTest()
{
	const double F = 100;
	const double years = 0.1;
	const double vol = 0.25;
	const double r = 0.1;
	GBM underlying_geometry(vol, years, r, 0);
	for (double K : std::vector<double>{80, 100, 120})
	{
		double cash_or_nothing = GetCashOrNothingPrice(K, F, r, years, vol);
		std::cout << "Test K = " << K << " : (" << cash_or_nothing << ", " << 1 - cash_or_nothing << ")" << std::endl;

		double x = std::log(F / K);
		
		//std::size_t m = 1;
		std::size_t m = 4;
		//int k1 = -1;
		int k1 = -8;
		//int k2 = 2;
		int k2 = 16;
		//std::size_t j = 4;
		std::size_t j = 7;
		std::size_t N = 1UL << j;
		//std::size_t j_bar = 3;
		std::size_t j_bar = 6;

		/// SLOW METHOD

		std::vector<double> c_n;
		const std::size_t two_to_the_j_minus_1 = 1UL << (j - 1);
		for (int k = k1; k <= k2; ++k)
		{
			double c_m_k = 0;
			for (std::size_t jj = 1; jj <= two_to_the_j_minus_1; ++jj)
			{
				std::complex<double> f_hat = underlying_geometry.GetChar((2 * jj - 1) * Swift::PI * (1UL << m) / N, x);
				std::complex<double> cv = f_hat
				    * std::exp(1i * static_cast<double>(k) * Swift::PI * static_cast<double>(2 * jj - 1) / static_cast<double>(N));
				c_m_k += cv.real();
			}
			c_n.push_back(c_m_k * std::pow(2.0, m / 2.0) / two_to_the_j_minus_1);
		}
		
		//fftw_complex* frequency_values, * time_values;
		//fftw_plan p;
		//frequency_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
		//time_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
		//p = fftw_plan_dft_1d(N, frequency_values, time_values, FFTW_BACKWARD, FFTW_ESTIMATE);
		//std::vector<std::complex<double>> x(N);
		//for (std::size_t i = 0; i < N; ++i)
		//{
		//	std::complex<double> value = underlying_geometry.GetChar((2 * j + 1) * Swift::PI * (1 << m) / N);
		//	frequency_values[i][0] = value.real();
		//	frequency_values[i][1] = value.imag();
		//}
		//fftw_execute(p);
		//std::vector<std::complex<double>> std_time_values;
		//for (std::size_t i = 0; i < N; ++i)
		//	std_time_values.emplace_back(time_values[i][0], time_values[i][1]);
		double approx_cash_or_nothing = 0.0;
		for (int k = k1; k <= k2; ++k)
		{
			double payoff = GetCashOrNothingCallPayOff(m, k, j_bar);
			approx_cash_or_nothing += c_n[k - k1] * payoff;
			std::cout << "k " << k << " payoff " << payoff << " shannon wavelet " << c_n[k - k1] << std::endl;
		}
		approx_cash_or_nothing *= std::exp(-r * years);
		std::cout << "Original " << cash_or_nothing << " Approx " << approx_cash_or_nothing << " Diff " << cash_or_nothing - approx_cash_or_nothing << std::endl;
	}
}

void HestonTest()
{
	const double F = 80;
	const double years = 0.5;
	const double r = 0.03;
	//const double q = 0.0;
	const bool is_call = true;
	HestonParameters heston_parameters(/*k*/ 1.0, /*v_bar*/ 0.05, /*sigma*/ 0.2, /*rho*/ -0.7, /*v_0*/ 0.04);
	HestonDistribution underlying_geometry(heston_parameters, r);
	for (double K : std::vector<double>{ 80 })
	{
		double heston_real = 5.1772;
		std::cout << "Test K = " << K << " : (" << heston_real << ", " << 1 - heston_real << ")" << std::endl;

		double x = std::log(F / K);

		//std::size_t m = 1;
		std::size_t m = 4;
		//int k1 = -1;
		int k1 = -8;
		//int k2 = 2;
		int k2 = 16;
		//std::size_t j = 4;
		std::size_t j = 7;
		std::size_t N = 1 << j;
		//std::size_t j_bar = 3;
		std::size_t j_bar = 6;

		/// SLOW METHOD

		std::vector<double> c_n;
		const std::size_t two_to_the_j_minus_1 = 1 << (j - 1);
		for (int k = k1; k <= k2; ++k)
		{
			double c_m_k = 0;
			for (std::size_t jj = 1; jj <= two_to_the_j_minus_1; ++jj)
			{
				std::complex<double> f_hat = underlying_geometry.GetChar((2 * jj - 1) * Swift::PI * (1 << m) / N, x, years);
				std::complex<double> cv = f_hat
					* std::exp(1i * static_cast<double>(k)* Swift::PI* static_cast<double>(2 * jj - 1) / static_cast<double>(N));
				c_m_k += cv.real();
			}
			c_n.push_back(c_m_k * std::pow(2.0, m / 2.0) / two_to_the_j_minus_1);
		}

		//fftw_complex* frequency_values, * time_values;
		//fftw_plan p;
		//frequency_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
		//time_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
		//p = fftw_plan_dft_1d(N, frequency_values, time_values, FFTW_BACKWARD, FFTW_ESTIMATE);
		//std::vector<std::complex<double>> x(N);
		//for (std::size_t i = 0; i < N; ++i)
		//{
		//	std::complex<double> value = underlying_geometry.GetChar((2 * j + 1) * Swift::PI * (1 << m) / N);
		//	frequency_values[i][0] = value.real();
		//	frequency_values[i][1] = value.imag();
		//}
		//fftw_execute(p);
		//std::vector<std::complex<double>> std_time_values;
		//for (std::size_t i = 0; i < N; ++i)
		//	std_time_values.emplace_back(time_values[i][0], time_values[i][1]);
		double approx_heston = 0.0;
		for (int k = k1; k <= k2; ++k)
		{
			double payoff = GetEuropeanPayoff(K, m, k, k1, k2, j_bar, is_call);
			approx_heston += c_n[k - k1] * payoff;
			std::cout << "k " << k << " payoff " << payoff << " shannon wavelet " << c_n[k - k1] << std::endl;
		}
		approx_heston *= std::exp(-r * years);
		std::cout << "Original " << heston_real << " Approx " << approx_heston << " Diff " << heston_real - approx_heston << std::endl;
	}
}
	
void BlackAndScholesTest()
{
	const double F = 100;
	const double years = 1;
	const double r = 0.1;
	//const double q = 0.0;
	const bool is_call = true;
	const double vol = 0.25;
	GBM underlying_geometry(vol, years, r, 0);
	for (double K : std::vector<double>{ 100 })
	{
		double real_value = 14.975790;
		std::cout << "Test K = " << K << " : (" << real_value << ", " << 1 - real_value << ")" << std::endl;

		double x = std::log(F / K);

		std::size_t m = 3;
		int k1 = -19;
		int k2 = 20;
		std::size_t j = 11;
		std::size_t N = 1 << j;
		std::size_t j_bar = 10;

		/// SLOW METHOD

		std::vector<double> c_n;
		const std::size_t two_to_the_j_minus_1 = 1 << (j - 1);
		for (int k = k1; k <= k2; ++k)
		{
			double c_m_k = 0;
			for (std::size_t jj = 1; jj <= two_to_the_j_minus_1; ++jj)
			{
				std::complex<double> f_hat = underlying_geometry.GetChar((2 * jj - 1) * Swift::PI * (1 << m) / N, x);
				std::complex<double> cv = f_hat
					* std::exp(1i * static_cast<double>(k)* Swift::PI* static_cast<double>(2 * jj - 1) / static_cast<double>(N));
				c_m_k += cv.real();
			}
			c_n.push_back(c_m_k * std::pow(2.0, m / 2.0) / two_to_the_j_minus_1);
		}

		//fftw_complex* frequency_values, * time_values;
		//fftw_plan p;
		//frequency_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
		//time_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
		//p = fftw_plan_dft_1d(N, frequency_values, time_values, FFTW_BACKWARD, FFTW_ESTIMATE);
		//std::vector<std::complex<double>> x(N);
		//for (std::size_t i = 0; i < N; ++i)
		//{
		//	std::complex<double> value = underlying_geometry.GetChar((2 * j + 1) * Swift::PI * (1 << m) / N);
		//	frequency_values[i][0] = value.real();
		//	frequency_values[i][1] = value.imag();
		//}
		//fftw_execute(p);
		//std::vector<std::complex<double>> std_time_values;
		//for (std::size_t i = 0; i < N; ++i)
		//	std_time_values.emplace_back(time_values[i][0], time_values[i][1]);
		double approx_heston = 0.0;
		for (int k = k1; k <= k2; ++k)
		{
			double payoff = GetEuropeanPayoff(K, m, k, k1, k2, j_bar, is_call);
			approx_heston += c_n[k - k1] * payoff;
			std::cout << "k " << k << " payoff " << payoff << " shannon wavelet " << c_n[k - k1] << std::endl;
		}
		approx_heston *= std::exp(-r * years);
		std::cout << "Original " << real_value << " Approx " << approx_heston << " Diff " << real_value - approx_heston << std::endl;
	}
}

void EuropeanPayoffTest()
{
	bool is_call = true;
	std::size_t J = 3;
	int k1 = -1;
	int k2 = 2;
	std::size_t m = 1;
	double K = 80;
	int k = -1;
	std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << GetEuropeanPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
	k = 0;
	std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << GetEuropeanPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
	k = 1;
	std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << GetEuropeanPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
	k = 2;
	std::cout << "(k1,k,k2): " << k1 << ", " << k << ", " << k2 << " Payoff " << GetEuropeanPayoff(K, m, k, k1, k2, J, is_call) << std::endl;
}

void HestonCharTest()
{
	double r = 0.02;
	HestonParameters heston_parameters(/*k*/ 1.0, /*v_bar*/ 0.05, /*sigma*/ 0.2, /*rho*/ -0.7, /*v_0*/ 0.04);
	HestonDistribution underlying_geometry(heston_parameters, r);
	for (double u : std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0})
	{
		std::complex<double> char_val = underlying_geometry.GetChar(u, 0, 15);
		std::complex<double> heston_char_val = underlying_geometry.GetHestonChar(u, 15);
		std::cout << "u = " << u << " Re(char) = " << char_val.real() << " Heston char " << heston_char_val.real() << std::endl;
	}
}

int main() {
	//TestSI();
	//TestFFTW3();
	//FFTUT();
	//BitReverseUT();
	//CashOrNothingTest();
	//EuropeanPayoffTest();
	//HestonTest();
	HestonCharTest();
	//BlackAndScholesTest();
}
