#include <iostream>
#include <boost/test/unit_test.hpp>
#include <SWIFT/quick_callibration_swift.h>

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
    int n = 5;
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
    std::vector<double> payoffcoefic;

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

        /* Density Coefficients*/
        coefic.resize(max_k - min_k + 1);
        for(int k=min_k;k<=max_k;k++)
        {
            coefic[k-min_k]=coef(k);
        }

        //Payoff
        //Suponemos K1<=0 y por tanto k1 barra=0, de esta forma el maximo Nk=k2-k1
        maxabsk=abs(max_k-min_k);
        Jp=ceil(log2(PI*maxabsk));
        pw2Jp=pow(2,Jp-1);
        /*Payoff with FFT*/
        payoffcoefic.resize(max_k-min_k+1);

        FFT_ForPayoffCoefficients(payoffcoefic);
    }


    double PriceCall()
    {
        /*Call price*/
        double call_shannon=0;

        for(int k=min_k;k<=max_k;k++) {
//            call_shannon+=coefic[k-min_k]*V(k);
            call_shannon+=coefic[k-min_k]*payoffcoefic[k-min_k];
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

    void FFT_ForPayoffCoefficients(std::vector<double>& FFT_payoffcoefVector)
    {
        int i,j;
        double sup=k_2/pw2n; //Suponemos que estamos en el escenario de intervalo: [0,k2/2^m]
        double inf=0;
        double A,B,Co,sa,sb,ca,cb,I11,I12,I21,I22,ea,eb;
        double *d,*e;

        int Npoints = pw2Jp;
        double *in,*out,*in2,*out2;
        fftw_plan plan1;
        fftw_plan plan2;

        ea=exp(inf);
        eb=exp(sup);

        d = (double*) fftw_malloc(sizeof(double)*pw2Jp);
        e = (double*) fftw_malloc(sizeof(double)*pw2Jp);

        in = (double*) fftw_malloc(sizeof(double)*Npoints);
        out = (double*) fftw_malloc(sizeof(double)*Npoints);
        in2 = (double*) fftw_malloc(sizeof(double)*Npoints);
        out2 = (double*) fftw_malloc(sizeof(double)*Npoints);

        for(j=0;j<pw2Jp;j++)
        {
            Co=((2*j+1.)/(2.*pw2Jp))*PI;
            B=(1./(Co*pw2n));
            A=(Co*pw2n)/(1.+pow(Co*pw2n,2));
            sb=sin(Co*pw2n*sup);
            sa=sin(Co*pw2n*inf);
            cb=cos(Co*pw2n*sup);
            ca=cos(Co*pw2n*inf);
            I11=eb*sb-ea*sa+B*eb*cb-B*ea*ca;
            I12=-eb*cb+ea*ca+B*eb*sb-B*ea*sa;
            I21=sb-sa;
            I22=ca-cb;
            d[j]=A*I11-B*I21;
            e[j]=A*I12-B*I22;
        }

        //Calculem amb FFT

        plan1 = fftw_plan_r2r_1d(Npoints, d, out, FFTW_REDFT10, FFTW_ESTIMATE);       //Here we set which kind of transformation we want to perform
        fftw_execute(plan1);                                                             //Execution of FFT

        plan2 = fftw_plan_r2r_1d(Npoints, e, out2, FFTW_RODFT10, FFTW_ESTIMATE);       //Here we set which kind of transformation we want to perform
        fftw_execute(plan2);                                                             //Execution of FFT

        //Para k=0
        FFT_payoffcoefVector[-min_k]=(1./pw2Jp)*(out[0]/2.);
        //printf("res[%d]=%lf\n",-min_k,FFT_payoffcoefVector[-min_k]);

        //Para k>0

        for(i=1;i<=max_k;i++)
        {
            FFT_payoffcoefVector[-min_k+i]=(1./pw2Jp)*((out[i]+out2[i-1])/2.); //out viene multiplicado por 2 !!!
            //printf("res[%d]=%lf\n",-min_k+i,FFT_payoffcoefVector[-min_k+i]);
        }

        //Para k<0

        for(i=1;i<=-min_k;i++)
        {
            FFT_payoffcoefVector[-min_k-i]=(1./pw2Jp)*((out[i]-out2[i-1])/2.); //out viene multiplicado por 2 !!!
            //printf("res[%d]=%lf\n",-min_k-i,FFT_payoffcoefVector[-min_k-i]);
        }

        fftw_destroy_plan(plan1);                                            //Destroy plan
        fftw_destroy_plan(plan2);
        free(d);
        free(e);
        free(in);                                                            //Free memory
        free(out);                                                           //Free memory
        free(in2);                                                           //Free memory
        free(out2);

        return;
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

//    std::cout << "My price " << my_price << " luis price " << luis_price << " before changes " << reference_without_changes << " real_price " << option_price << " diff " << my_price - luis_price << std::endl;
    BOOST_CHECK_CLOSE(my_price, luis_price,1e-3);
    BOOST_CHECK_CLOSE(my_price, option_price,1e-3);
}

BOOST_FIXTURE_TEST_CASE(EuropeanVega, LuisCodeSimulator)
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

    Swift::SwiftEvaluator eval(params, distribution, contract);
    std::vector<double> my_vega = eval.GetGradient(S0, K, r, 0.0, is_call);
    double option_vega = GetBSEuropeanVega(K, S0, r, T, sigma, 0.0);

    std::cout << "My vega " << my_vega[0] << " real_vega " << option_vega << " diff " << my_vega[0] - option_vega << std::endl;
    BOOST_CHECK_CLOSE(my_vega[0], option_vega,1e-3);
}

BOOST_FIXTURE_TEST_CASE(Fast, LuisCodeSimulator)
{
    bool is_call = true; // Luis only implemented calls
    Swift::SwiftParameters params{3, -7, 8, 8, 7};
    double locT = 0.119047619047619;
    GBM distribution(0.25, locT);
    EuropeanOptionContract contract;

    double locr = 0.02;
    double locK = 0.9371;
    double locS = 1.0;

    Swift::SwiftEvaluator eval(params, distribution, contract);
    double my_price = eval.GetPrice(locS, locK, locr, 0.0, is_call);
    Swift::SwiftInvariantData invariant_data(params, locT, contract, is_call, locS, {locK}, locr, 0.0);
    Swift::QuickCallibrationSwiftEvaluator quick_eval(invariant_data, params, distribution, contract);
    double my_quick_price = quick_eval.GetPrice(locS, {locK}, locr, 0.0)[0];

    std::cout << "my_price " << my_price << " my_quick_price " << my_quick_price << " diff " << my_price - my_quick_price << " factor " << my_quick_price / my_price << std::endl;
    BOOST_CHECK_CLOSE(my_price, my_quick_price,1e-3);
}

BOOST_AUTO_TEST_CASE(HestonGradient)
{
    bool is_call = true; // Luis only implemented calls
    Swift::SwiftParameters params{5, -51, 43, 7, 7};
    double T = 0.119048;
    HestonParameters heston_parameters{3.00, 0.10, 0.25, -0.8, 0.08};
    HestonDistribution distribution(heston_parameters, T);
    EuropeanOptionContract contract;

    double r = 0.02;
    double K = 0.9956;
    double S = 1.0;

    // K 0.9371 T 0.119048 jacobian 0.000134411 0.0280474 0.00330929 -0.00115206 0.151462
    // K 0.9956 T 0.119048 jacobian 0.0280474 0.00330929 -0.00115206 0.151462 0.000226969
    // K 1.0427 T 0.119048 jacobian 0.00330929 -0.00115206 0.151462 0.000226969 0.0375877
    // K 1.2287 T 0.119048 jacobian -0.00115206 0.151462 0.000226969 0.0375877 -0.000662676
    // K 1.3939 T 0.119048 jacobian 0.151462 0.000226969 0.0375877 -0.000662676 -6.89813e-05

    Swift::SwiftEvaluator eval(params, distribution, contract);
    std::vector<double> my_gradient = eval.GetGradient(S, K, r, 0.0, is_call);
    Swift::SwiftInvariantData invariant_data(params, T, contract, is_call, S, {K}, r, 0.0);
    Swift::QuickCallibrationSwiftEvaluator quick_eval(invariant_data, params, distribution, contract);
    std::vector<double> my_quick_gradient = quick_eval.GetGradient(S, {K}, r, 0.0)[0];
    std::vector<double> real_gradient{-0.00209303, 0.122985, 9.70534e-05, 0.0648615, 0.00843795};
    std::vector<double> my_cui_gradient = GetHestonEuropeanGradientCuiMyChar(heston_parameters, S, K, r, 0.0, T);

    std::cout << "my_gradient " << my_gradient[0] << ", " << my_gradient[1] << ", " << my_gradient[2] << ", " << my_gradient[3] << ", " << my_gradient[4] << std::endl
              << " my_quick_gradient " << my_quick_gradient[0] << ", " << my_quick_gradient[1] << ", " << my_quick_gradient[2] << ", " << my_quick_gradient[3] << ", " << my_quick_gradient[4] << std::endl
              << " my_cui_gradient " << my_cui_gradient[0] << ", " << my_cui_gradient[1] << ", " << my_cui_gradient[2] << ", " << my_cui_gradient[3] << ", " << my_cui_gradient[4] << std::endl
              << " real " << real_gradient[0] << ", " << real_gradient[1] << ", " << real_gradient[2] << ", " << real_gradient[3] << ", " << real_gradient[4] << std::endl;
//    BOOST_CHECK_CLOSE(my_price, my_quick_price,1e-3);
}

BOOST_AUTO_TEST_SUITE_END()