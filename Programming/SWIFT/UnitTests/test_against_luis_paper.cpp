#include <iostream>
#include <boost/test/unit_test.hpp>

#include "SWIFT/distributions.h"
#include "SWIFT/density_coefficients_calculators.h"
#include "SWIFT/known_distribution_contract_combinations.h"
#include "SWIFT/option_contracts.h"
#include "SWIFT/swift.h"

BOOST_AUTO_TEST_SUITE(ORT16UT)


#define PI M_PI /* PI to machine precision, defined in math.h */

/* Payoff coefficients */
//Estamos suponiendo k1<=0 por eso el limite inferior de I1 e I2 es 0

struct LuisCodeSimulator
{
    int n = 2;
    double r = 0.1;
    double L = 10.;
    double S0 = 100.;
    double real_price = 6.6383090775296700; const double reference_without_changes = 6.732761537207680; const double K=120; const double T=1; const double sigma=0.25;
    double x0 = log(S0/K);


    double pw2n;
    double pw2n2;
    double c1;
    double c2;
    double ext_inf;
    double ext_sup;
    int min_k;
    int max_k;
    int k_2;
    double Mm;
    int Jd;
    double pwJd;
    double maxabsk;
    double Jp;
    double pw2Jp;

    std::vector<double> ueval;
    std::vector<std::complex<double>> f_hat;
    std::vector<double> coefic;

    LuisCodeSimulator()
    {
        pw2n=pow(2,n);
        pw2n2=sqrt(pw2n);

        //Intervalo a priori y k mínima y máxima
        c1=x0+(r-0.5*pow(sigma,2))*T;
        c2=pow(sigma,2)*T;
        // TODO: This effectively makes the payoff K and the density k different. I may want to add that.
        ext_inf=c1-L*sqrt(c2);
        ext_sup=c1+L*sqrt(c2);

        //ext_inf=-8.;
        //ext_sup=12.;

        min_k=ceil(pw2n*ext_inf);
        max_k=floor(pw2n*ext_sup);
        k_2=max_k;

        //J usado en la densidad
        double a=std::max(std::fabs(ext_inf),std::fabs(ext_sup));
        Mm=std::max(std::fabs(pw2n*a-min_k),std::fabs(pw2n*a+max_k));
        Jd=ceil(log2(PI*Mm));

        //Jd=10;
        pwJd=pow(2,Jd-1);

        /* Fourier transfom evaluation */
        ueval.reserve(pwJd);
        f_hat.reserve(pwJd);
        for(int i=0;i<pwJd;i++)
        {
            ueval.push_back(((2.*i+1.)*PI)/(2.*pwJd));
            f_hat.push_back(CharacteristicFunction(pw2n * ueval[i]));
        }

        /*Coefficients*/
        coefic.resize(max_k - min_k + 1);
        for(int k=min_k;k<=max_k;k++) {
            coefic[k-min_k]=coef(k);
        }

        //Payoff
        //Suponemos K1<=0 y por tanto k1 barra=0, de esta forma el maximo Nk=k2-k1
        maxabsk=abs(max_k-min_k);
        Jp=ceil(log2(PI*maxabsk));
        pw2Jp=pow(2,Jp-1);
    }


    double PriceCall()
    {
        /*Call price*/
        double call_shannon=0;

        for(int k=min_k;k<=max_k;k++) {
            call_shannon+=coefic[k-min_k]*V(k);
            //printf("V(%d)=%lf\n",k,V(k));
        }
        call_shannon*=K*exp(-r*T)*pw2n2;
        return call_shannon;
    }

/*Funciones de escala*/
    double phi(int loc_n, int k, double x)
    {
        double pwj,pwj2,fun;

        pwj=pow(2, loc_n);
        pwj2=sqrt(pwj);
        fun=pwj2*sin(PI*(pwj*x-k))/(PI*(pwj*x-k));

        return(fun);
    }

/* Calculo de los coeficientes de la aproximacion */
    double coef(int k)
    {
        double sum=0;
        for (int j = 0; j < pwJd; j++)
            sum+=(f_hat[j]*std::exp(k*ueval[j]*1i)).real();

        return(pw2n2 * sum / (pwJd));
    }

    std::complex<double> CharacteristicFunction(double w)
    {
        std::complex<double> p1;
        std::complex<double> p2;

        p1=std::exp(-1.0i*w*(r-0.5*std::pow(sigma,2))*T);
        p2=std::exp(-0.5*std::pow(sigma,2)*std::pow(w,2)*T);

        return(std::exp(-w*x0*1i)*p1*p2);
    }

/* Payoff coefficients */
//Estamos suponiendo k1<=0 por eso el limite inferior de I1 e I2 es 0
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
        min_k,
        k_2,
        static_cast<std::size_t>(Jd - 1),
        static_cast<std::size_t>(Jp - 1)
    };
    EuropeanOptionContract european_option_contract;

    for (int k = params.m_k1; k < params.m_k2; ++k)
    {
        double my_payoff_old = european_option_contract.GetPayoffNonKComponentOldPaper(k, params, is_call);
        double my_payoff_new = european_option_contract.GetPayoffNonKComponentNewPaper(k, params, is_call);
        double my_luis_payoff = european_option_contract.V(k, params, is_call);
        double luis_payoff = pw2n2 * V(k);
//        std::cout << "k " << k << " MyOld " << my_payoff_old << " MyNew " << my_payoff_new << " MyLuis_payoff " << my_luis_payoff << " luis_payoff " << luis_payoff << " diff_old " << my_payoff_old - luis_payoff<< " diff_new " << my_payoff_new - luis_payoff << std::endl;
        BOOST_CHECK_CLOSE(my_payoff_old, my_payoff_new,1e-3);
        BOOST_CHECK_CLOSE(my_payoff_old, my_luis_payoff,1e-3);
        BOOST_CHECK_CLOSE(my_payoff_old, luis_payoff,1e-3);
    }
}

BOOST_FIXTURE_TEST_CASE(EuropeanPriceVsLuisPrice, LuisCodeSimulator)
{
    bool is_call = true; // Luis only implemented calls
    Swift::SwiftParameters params{
        static_cast<std::size_t>(n),
        min_k,
        k_2,
        static_cast<std::size_t>(Jd - 1),
        static_cast<std::size_t>(Jp - 1)
    };
    GBM distribution(sigma, T);
    EuropeanOptionContract contract;

    double luis_price = PriceCall();
    Swift::SwiftEvaluator eval(params, distribution, contract);
    double my_price = eval.GetPrice(S0, K, r, 0.0, is_call);
    double option_price = GetBSEuropeanPrice(K, S0, r, T, sigma, true, 0.0);

    std::cout << "My price " << my_price << " luis price " << luis_price << " before changes " << reference_without_changes << " real_price " << option_price << " diff " << my_price - luis_price << std::endl;
    BOOST_CHECK_CLOSE(my_price, luis_price,1e-3);
}

BOOST_AUTO_TEST_SUITE_END()