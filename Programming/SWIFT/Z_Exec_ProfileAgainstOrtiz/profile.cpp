/*Version 1: la f_hat se evalua 1 sola vez sobre la malla de puntos para todos los coeficientes */
/*La formula de Viete se usa para calcular la densidad*/
/*Aplicamos FFT a coeficientes densidad y a coeficientes pay-off*/

// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <ccomplex>
#include "FFTW3/include_fftw3.h"

#include <SWIFT/swift.h>
#include <SWIFT/distributions.h>
#include <SWIFT/option_contracts.h>
#include <SWIFT/known_distribution_contract_combinations.h>

using namespace std::complex_literals;

#define PI M_PI /* PI to machine precision, defined in math.h */
#define TWOPI (2.0*PI)
#define max(x,y) (x>y?x:y)

int number_of_repetitions = 100;

#define n 5

#define S0 100.
//#define K 100.
//#define K 50.
#define K 200.
#define r 0.0
#define q 0.
//#define T 45.
#define T 0.04

#define  mu r
#define kappa 1.5768 // kappa
#define vmean 0.0398
#define sigma 0.5751 // Sigma
#define rho -0.5711
#define v0 0.0175

//#define reference 5.007054273487142e+01 /*K=50 L=10 COS 50000 terms*/
//#define reference 5.785155700563814e+00 /*K=100 L=10 COS 50000 terms*/

//#define reference 3.887893512658775e+01 /* T=30 */
#define reference 4.691153171658516e+01 /* T=45 */

#define L 12.
#define x0 log(S0/K)

double V(int);
double phi(int,int,double);
double coef(int);
std::complex<double> chf(double);
double I1(double,double,int,int);
double I2(double,double,int,int);
void FFT_ForPDFCoefficients(double []);
void FFT_ForPayoffCoefficients(double []);

double *ueval,pw2n,pw2n2;
int m;
std::complex<double> *f_hat; double trap_payoff;
int k_2,Jd,pwJd;
int Jp,pw2Jp,maxabsk;
int min_k,max_k;

//double ext_inf,ext_sup;
//void InitData()
//{
//    pw2n=pow(2,n);
//    pw2n2=sqrt(pw2n);
//
//    double p1,p2,p3,p4,p5;
//    double c1,c2;
//    //Intervalo
//    c1=mu*T+(1.-exp(-kappa*T))*(vmean-v0)/(2.*kappa)-0.5*vmean*T;
//    p1=sigma*T*kappa*exp(-kappa*T)*(v0-vmean)*(8.*kappa*rho-4.*sigma);
//    p2=kappa*rho*sigma*(1.-exp(-kappa*T))*(16.*vmean-8.*v0);
//    p3=2.*vmean*kappa*T*(-4.*kappa*rho*sigma+pow(sigma,2)+4.*pow(kappa,2));
//    p4=pow(sigma,2)*((vmean-2.*v0)*exp(-2.*kappa*T)+vmean*(6.*exp(-kappa*T)-7.)+2.*v0);
//    p5=8.*pow(kappa,2)*(v0-vmean)*(1.-exp(-kappa*T));
//    c2=(1./(8.*pow(kappa,3)))*(p1+p2+p3+p4+p5);
//
//    ext_inf=x0+c1-L*sqrt(std::fabs(c2));
//    ext_sup=x0+c1+L*sqrt(std::fabs(c2));
//    //printf("ext_inf=%lf\text_sup=%lf\n",ext_inf,ext_sup);
//
//    //Rango de k's para la densidad
//    min_k=ceil(pw2n*ext_inf); printf("min_k=%d\n",min_k);
//    max_k=floor(pw2n*ext_sup); printf("max_k=%d\n",max_k);
//    k_2=max_k;
//
//    //J usado en la densidad
//    double a,Mm;
//    a=max(std::fabs(ext_inf),std::fabs(ext_sup)); //printf("a=%lf\n",a);
//    Mm=max(std::fabs(pw2n*a-min_k),std::fabs(pw2n*a+max_k)); //printf("Mm=%lf\n",Mm);
//    Jd=ceil(log2(PI*Mm)); //printf("Jd=%d\n",Jd);
//
//    //printf("Jd>=%d\t\tJdnoredondeo=%lf\n",Jd,log2(PI*Mm));
//    //Jd=10;
//    pwJd=pow(2,Jd-1); //printf("pwJd=%d\n",pwJd);
//}


void run_ortiz()
{
    printf("Ortiz");
    double xi,sum,ext_inf,ext_sup,c1,c2;
    double p1,p2,p3,p4,p5;
    int i,k;
    FILE *out;
//    ,*outpayoff;
    double *coefic;
    double call_shannon;
    clock_t inicio,parada;
//    double max_error;
    double a,Mm;
//    ,D;
    double A;
    double *payoffcoefic;

    int cont;

    inicio=clock();
    for(cont=0; cont < number_of_repetitions;cont++){
        pw2n=pow(2,n);
        pw2n2=sqrt(pw2n);

        //Intervalo
        c1=mu*T + (1.-exp(-kappa * T)) * (vmean - v0) / (2. * kappa) - 0.5 * vmean * T;
        p1= sigma * T * kappa * exp(-kappa * T) * (v0 - vmean) * (8. * kappa * rho - 4. * sigma);
        p2= kappa * rho * sigma * (1. - exp(-kappa * T)) * (16. * vmean - 8. * v0);
        p3= 2. * vmean * kappa * T * (-4. * kappa * rho * sigma + pow(sigma, 2) + 4. * pow(kappa, 2));
        p4= pow(sigma, 2) * ((vmean - 2. * v0) * exp(-2. * kappa * T) + vmean * (6. * exp(-kappa * T) - 7.) + 2. * v0);
        p5= 8. * pow(kappa, 2) * (v0 - vmean) * (1. - exp(-kappa * T));
        c2= (1./(8.*pow(kappa, 3))) * (p1 + p2 + p3 + p4 + p5);

        ext_inf=x0+c1-L*sqrt(fabs(c2));
        ext_sup=x0+c1+L*sqrt(fabs(c2));
//        printf("ext_inf=%lf\text_sup=%lf\n",ext_inf,ext_sup);

        //Rango de k's para la densidad
        min_k=ceil(pw2n*ext_inf);
//        printf("min_k=%d\n",min_k);
        max_k=floor(pw2n*ext_sup);
//        printf("max_k=%d\n",max_k);
        k_2=max_k;

        //J usado en la densidad
        a=max(fabs(ext_inf),fabs(ext_sup));
//        printf("a=%lf\n",a);
        Mm=max(fabs(pw2n*a-min_k),fabs(pw2n*a+max_k));
//        printf("Mm=%lf\n",Mm);
        Jd=ceil(log2(PI*Mm)); //printf("Jd=%d\n",Jd);

        //printf("Jd>=%d\t\tJdnoredondeo=%lf\n",Jd,log2(PI*Mm));
        //Jd=10;
        pwJd=pow(2,Jd-1); //printf("pwJd=%d\n",pwJd);

        /* Fourier transfom evaluation */
        ueval=static_cast<double*>(malloc(pwJd*sizeof(double)));
        //f_hat=malloc(pwJd*sizeof(std::complex<double>));

        for(i=0;i<pwJd;i++)
        {
            ueval[i]=((2.*i+1.)*PI)/(2.*pwJd);
            //f_hat[i]=chf(pw2n*ueval[i]);
        }

        /*Coefficients*/
        coefic=static_cast<double*>(malloc((max_k-min_k+1)*sizeof(double)));

        /*Coefficients without FFT*/
        //for(k=min_k;k<=max_k;k++) coefic[k-min_k]=coef(k);

        /*Coefficients with FFT*/
        FFT_ForPDFCoefficients(coefic);

        //Payoff
        //Suponemos K1<=0 y por tanto k1 barra=0, de esta forma el maximo Nk=k2-k1

        maxabsk=abs(max_k-min_k);
        Jp=ceil(log2(PI*maxabsk)); //printf("Jp=%d\n",Jp);
        pw2Jp=pow(2,Jp-1);
//        std::cout << std::string("m, min_k, k_2, Jd, Jp ") << m << ", " << min_k << ", "<< k_2 << ", " << Jd << ", " << Jp << std::endl;

        /*Payoff with FFT*/
        payoffcoefic=static_cast<double*>(malloc((max_k-min_k+1)*sizeof(double)));

        FFT_ForPayoffCoefficients(payoffcoefic);

        /*Call price*/
        //Sin FFT
        call_shannon=0;
        //for(k=min_k;k<=max_k;k++) call_shannon+=coefic[k-min_k]*V(k);

        //Con FFT
        for(k=min_k;k<=max_k;k++) call_shannon+=coefic[k-min_k]*payoffcoefic[k-min_k];
        //for(k=min_k;k<=5;k++) call_shannon+=coefic[k-min_k]*payoffcoefic[k-min_k];

        call_shannon*=K*exp(-r*T)*pw2n2;

    }
    parada=clock();
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

    printf("call_shannon=%.15lf\terror=%.2e\n",call_shannon,fabs(call_shannon-reference));
    printf("Tiempo ejecucion: %.3f\n",((double)(parada-inicio))/CLOCKS_PER_SEC);
    std::cout << std::endl;

    /*Area bajo la funcion de densidad=1?*/
    A=0.5*(coefic[0]+coefic[max_k-min_k]);
    for(k=min_k+1;k<max_k;k++) A+=coefic[k-min_k];
    A*=1./pw2n2;
    printf("A=%.15e\t%.2e\n",A,fabs(A-1));

    /*Aproximacion a la funcion de densidad*/
    out=fopen("density_Shannon_viete_density.dat","w");
    for(xi=ext_inf+0.01;xi<ext_sup;xi+=0.01)
    {
        sum=0;
        for(k=min_k;k<=max_k;k++)
        {
            sum+=coefic[k-min_k]*phi(n,k,xi);
        }
        fprintf(out,"%lf %e\n",xi,sum);
    }
    fclose(out);

    free(ueval);
    //free(f_hat);
    free(coefic);
    free(payoffcoefic);
#pragma GCC diagnostic pop
//    return 0;
}

/*Funciones de escala*/
double phi(int j,int k,double x)
{
    double pwj,pwj2,fun;

    pwj=pow(2,j);
    pwj2=sqrt(pwj);
    fun=pwj2*sin(PI*(pwj*x-k))/(PI*(pwj*x-k));

    return(fun);
}

/* Calculo de los coeficientes de la aproximacion */
double coef(int k)
{
    int j;
    double sum;

    sum=0;
    for(j=0;j<pwJd;j++) sum+=(f_hat[j]*std::exp(k*ueval[j]*1i)).real();

    return(pw2n2*sum/(pwJd));
}

std::complex<double> chf(double w)
{
    std::complex<double> D=std::sqrt(std::pow(kappa + rho * sigma * w * 1i, 2) + (std::pow(w, 2) - w * 1i) * std::pow(sigma, 2));
    std::complex<double> G= (kappa + rho * sigma * w * 1i - D) / (kappa + rho * sigma * w * 1i + D);
    std::complex<double> e1=-mu*w*T*1i+ (v0/pow(sigma, 2)) * ((1. - std::exp(-D * T)) / (1. - G * std::exp(-D * T))) * (kappa + rho * sigma * w * 1i - D);
    std::complex<double> e2= ((kappa * vmean) / pow(sigma, 2)) * (T * (kappa + rho * sigma * w * 1i - D) - 2. * std::log((1. - G * std::exp(-D * T)) / (1. - G)));
    return(std::exp(-w*x0*1i)*std::exp(e1)*std::exp(e2));
}

/* Payoff coefficients */
//Estamos suponiendo k1<=0 por eso el limite inferior de I1 e I2 es 0
double V(int k)
{
    int j;
    double sum,sup;

    sup=k_2/pw2n;
    sum=0;
    for(j=1;j<=pw2Jp;j++) sum+=I1(0,sup,k,j)-I2(0,sup,k,j);

    return ((1./pw2Jp)*sum);
}

double I1(double a,double b,int k,int j)
{
    double Co,C2,C3,C4;

    Co=((2*j-1.)/(2.*pw2Jp))*PI;
    C2=(Co*pw2n)/(1.+pow(Co*pw2n,2));
    C3=exp(b)*sin(Co*(pw2n*b-k))-exp(a)*sin(Co*(pw2n*a-k));
    C4=(1./(Co*pw2n))*(exp(b)*cos(Co*(pw2n*b-k))-exp(a)*cos(Co*(pw2n*a-k)));

    return(C2*(C3+C4));
}

double I2(double a,double b,int k,int j)
{
    double Co,C2;

    Co=((2*j-1.)/(2.*pw2Jp))*PI;
    C2=(1./(Co*pw2n))*(sin(Co*(pw2n*b-k))-sin(Co*(pw2n*a-k)));

    return(C2);
}

void FFT_ForPDFCoefficients(double FFT_CoefVector[])
{
    std::complex<double> sum1;

    int i;
    int Npoints = 2*pwJd;
    fftw_complex *in, *out1;
    fftw_plan plan1;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);                 //allocating memory
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);               //allocating memory

    for(i=0; i<pwJd;i++)
    {
        std::complex<double> current_f_hat=chf(ueval[i]*pw2n);
        in[i][0] = current_f_hat.real();
        in[i][1] = current_f_hat.imag();
    }
    for(i=pwJd; i < Npoints;i++)
    {
        in[i][0] = 0.;
        in[i][1] = 0.;
    }

    //Pels coeficients positius
    plan1 = fftw_plan_dft_1d(Npoints, in, out1, FFTW_BACKWARD, FFTW_ESTIMATE);       //Here we set which kind of transformation we want to perform
    fftw_execute(plan1);                                                             //Execution of FFT
    for(i=1;i<=max_k;i++)
    {
        sum1=std::exp(M_PI*i*1i/ static_cast<double >(Npoints))*std::complex<double>(out1[i][0], out1[i][1]);
        FFT_CoefVector[-min_k+i]=pw2n2*(sum1.real())/pwJd;
    }
    //Pels coeficients negatius
    plan1 = fftw_plan_dft_1d(Npoints, in, out1, FFTW_FORWARD, FFTW_ESTIMATE);        //Here we set which kind of transformation we want to perform
    fftw_execute(plan1);                                                              //Execution of FFT
    for(i=0;i<=-min_k;i++)
    {
        sum1=std::exp(-M_PI*i*1i/ static_cast<double>(Npoints))*std::complex<double>(out1[i][0], out1[i][1]);
        FFT_CoefVector[-min_k-i]=pw2n2*(sum1.real())/pwJd;
    }
    fftw_destroy_plan(plan1);                                                 //Destroy plan
    fftw_free(in);                                                            //Free memory
    fftw_free(out1);                                                          //Free memor
    return;
}

void FFT_ForPayoffCoefficients(double FFT_payoffcoefVector[])
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

void run_romo()
{
    printf("Romo");
    volatile double call_shannon;
    double loc_x0 = Distribution::GetXCompression(S0, K, r, q, T);
    HestonParameters heston_parameters(kappa, vmean, sigma, rho, v0);
    HestonDistribution distribution(heston_parameters, T);
//    Swift::SwiftParameters params(n, min_k, max_k, Jd - 1, Jp - 1);
    Swift::SwiftParameters params(n, distribution, loc_x0, loc_x0);
    EuropeanOptionContract contract;
    Swift::SwiftEvaluator eval(params, distribution, contract);

    clock_t inicio=clock();
    for (int cont=1;cont<=number_of_repetitions;cont++){
        call_shannon = eval.GetPrice(S0, K, r, q, true);
    }
    clock_t parada=clock();
    printf("call_shannon=%.15lf\terror=%.2e\n", call_shannon, std::fabs(call_shannon - reference));
    printf("Tiempo ejecucion: %.3f\n", ((double) (parada - inicio)) / CLOCKS_PER_SEC);
    std::cout << std::endl;
}

void run_cui()
{
    printf("Cui");
    clock_t inicio=clock();
    volatile double call_shannon;
    int cont=1;
    for (; cont<=number_of_repetitions; ++cont) {
        HestonParameters heston_parameters(kappa, vmean, sigma, rho, v0);
        call_shannon = GetHestonEuropeanPriceCuiMyChar(heston_parameters, S0, K, r, q, T, 200);
    }
    clock_t parada=clock();
    printf("call_shannon=%.15lf\terror=%.2e\n", call_shannon, std::fabs(call_shannon - reference));
    printf("Tiempo ejecucion: %.3f\n", ((double) (parada - inicio)) / CLOCKS_PER_SEC);
    std::cout << std::endl;
}

int main()
{
//    run_ortiz();
    run_romo();
    run_cui();
}