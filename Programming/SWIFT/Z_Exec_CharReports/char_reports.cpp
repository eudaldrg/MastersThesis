#include <iostream>
#include "fstream"
#include "vector"
#include "SWIFT/distributions.h"

int main()
{
    double tau = 15;
//    HestonParameters heston_parameters(/*k*/ 1.0, /*v_bar*/ 0.05, /*sigma*/ 0.2, /*rho*/ -0.7, /*v_0*/ 0.04);
    HestonParameters heston_parameters( 3.0, 0.1, 0.25, -0.8, 0.08);
    HestonDistribution underlying_geometry(heston_parameters, tau);

    double S = 1.0, K = 1.1, r = 0.02, q=0.0;
    double x = HestonDistribution::GetXCompression(S, K, r,q, tau);

    std::vector<double> u_tests;
    for (std::size_t i = 0; i <= 40; ++i)
        u_tests.push_back(i * 0.1);

    std::ofstream report;
    report.open("/home/eudald/Projects/MastersThesis/Programming/SWIFT/Reports/cui_report.csv");
    if (!report.is_open())
        std::cout << "Couldn't open cui" << std::endl;
    for (double u : u_tests)
    {
        std::complex<double> complex_num = underlying_geometry.GetCuiChar(u, tau) / HestonDistribution::GetCharPositionFactor(x,u);
        report << u << ", " << complex_num.real() << ", " << complex_num.imag() << ", " << std::norm(complex_num) << ", " << std::arg(complex_num) << std::endl;
    }
    report.close();
    report.open("/home/eudald/Projects/MastersThesis/Programming/SWIFT/Reports/heston_report.csv");
    if (!report.is_open())
        std::cout << "Couldn't open heston" << std::endl;
    for (double u : u_tests)
    {
        std::complex<double> complex_num = underlying_geometry.GetHestonChar(u, tau) / HestonDistribution::GetCharPositionFactor(x,u);
        report << u << "," << complex_num.real() << "," << complex_num.imag() << "," << std::norm(complex_num) << "," << std::arg(complex_num) << std::endl;
    }
    report.close();
}
