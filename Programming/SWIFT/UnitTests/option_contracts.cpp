#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/distributions.h"
#include "SWIFT/density_coefficients_calculators.h"
#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"

BOOST_AUTO_TEST_SUITE(PayoffUT)

/* Payoff coefficients */
//Estamos suponiendo k1<=0 por eso el limite inferior de I1 e I2 es 0

#define PI M_PI /* PI to machine precision, defined in math.h */

struct LuisCodeSimulator
{
    int n = 2;
    double pw2n = pow(2, n);
    double pw2n2 = sqrt(pw2n);

    const double r = 0.1;
    const double L = 10.;
    const double S0 = 100.;
    const double reference = 6.6383090775296700; const double K=120; const double T=1; const double sigma=0.25;
    double x0 = log(S0/K);

    //Intervalo a priori y k mínima y máxima
    double c1 = x0+(r-0.5*pow(sigma,2))*T;
    double c2 = pow(sigma,2)*T;
    double ext_inf = c1-L*sqrt(c2);
    double ext_sup = c1+L*sqrt(c2);

    //ext_inf=-8.;
    //ext_sup=12.;

    int min_k=ceil(pw2n*ext_inf);
    int max_k=floor(pw2n*ext_sup);
    int k_2 = max_k;

    //J usado en la densidad
    double ad = std::max(std::fabs(ext_inf),std::fabs(ext_sup));
    double Mm =std::max(std::fabs(pw2n*ad-min_k),std::fabs(pw2n*ad+max_k));
    int Jd = std::roundl(std::ceil(log2(MY_PI*Mm)));
    double pwJd = pow(2,Jd-1);

    double maxabsk = abs(max_k-min_k);
    double Jp = ceil(log2(MY_PI * maxabsk));
    double pw2Jp = pow(2,Jp-1);

    double V(int k)
    {
        int j;
        double sup=k_2/pw2n;
        double sum=0;
        for (j=1;j<=pw2Jp;j++)
            sum+=I1(0,sup,k,j)-I2(0,sup,k,j);
        return ((1./pw2Jp)*sum);
    }

    double I1(double a,double b,int k,int j)
    {
        double C=((2*j-1.)/(2.*pw2Jp))*PI;
        double C2=(C*pw2n)/(1.+pow(C*pw2n,2));
        double C3=exp(b)*sin(C*(pw2n*b-k))-exp(a)*sin(C*(pw2n*a-k));
        double C4=(1./(C*pw2n))*(exp(b)*cos(C*(pw2n*b-k))-exp(a)*cos(C*(pw2n*a-k)));
        return(C2*(C3+C4));
    }

    double I2(double a,double b,int k,int j)
    {
        double C=((2*j-1.)/(2.*pw2Jp))*PI;
        double C2=(1./(C*pw2n))*(sin(C*(pw2n*b-k))-sin(C*(pw2n*a-k)));
        return(C2);
    }
};

BOOST_FIXTURE_TEST_CASE(EuropeanPayoffComparisonWithLuisUT, LuisCodeSimulator)
{
    bool is_call = true; // Luis only implemented calls
    Swift::SwiftParameters params{
        static_cast<std::size_t>(n),
        1 - k_2,
        k_2,
        static_cast<std::size_t>(Jd - 1),
        static_cast<std::size_t>(Jp - 1)
    };
    EuropeanOptionContract european_option_contract;

    for (int k = params.m_k1; k < params.m_k2; ++k)
    {
        double my_payoff_old = european_option_contract.GetPayoffNonKComponentOldPaper(k, params, is_call);
        double my_payoff_new = european_option_contract.GetPayoffNonKComponentOldPaper(k, params, is_call);
        double my_luis_payoff = european_option_contract.V(k, params, is_call);
        double luis_payoff = V(k);
        std::cout << "k " << k << " MyOld " << my_payoff_old << " MyNew " << my_payoff_new << " MyLuis_payoff " << my_luis_payoff << " luis_payoff " << luis_payoff << " diff_old " << my_payoff_old - luis_payoff<< " diff_new " << my_payoff_new - luis_payoff << std::endl;
        BOOST_CHECK_CLOSE(my_luis_payoff, luis_payoff,1e-3);
    }
}

BOOST_AUTO_TEST_CASE(EuropeanPayoffComparisonUT)
{
    bool is_call = false;
    Swift::SwiftParameters params{
        4,
        -8,
        164,
        12,
        12
    };
    EuropeanOptionContract european_option_contract;

    for (int k = params.m_k1; k < params.m_k2; ++k)
    {
        double old_payoff = european_option_contract.GetPayoffNonKComponentOldPaper(k, params, is_call);
        double new_payoff = european_option_contract.GetPayoffNonKComponentNewPaper(k, params, is_call);
        std::cout << "old " << old_payoff << " new " << new_payoff << " diff " << old_payoff - new_payoff << std::endl;
        BOOST_CHECK_CLOSE(old_payoff, new_payoff,1e-3);
    }
}

BOOST_AUTO_TEST_SUITE_END()