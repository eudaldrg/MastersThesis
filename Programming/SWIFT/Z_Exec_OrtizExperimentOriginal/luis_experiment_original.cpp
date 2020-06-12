/*Version 1: la f_hat se evalua 1 sola vez sobre la malla de puntos para todos los coeficientes */
/*La formula de Viete se usa para calcular la densidad*/

// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <complex>

#define PI M_PI /* PI to machine precision, defined in math.h */
#define TWOPI (2.0*PI)

const int n = 2;
const double S0 = 100.;
//const double K = 120.;
const double r = 0.1;
//const double T = 50;
//const double sigma = 0.25;
const double L = 10.;

//const double reference = 14.9757907783113000; const double K=100; const double T=1; const double sigma=0.25;
//const double reference = 0.044577814073288; const double K=120; const double T=0.1; const double sigma=0.25;
const double reference = 6.6383090775296700; const double K=120; const double T=1; const double sigma=0.25;
//const double reference = 99.2025928525532000; const double K=120; const double T=50; const double sigma=0.25;
//#define reference 99.9945609694213000 /*K=120 T=100 sigma=0.25*/
//#define reference 2.4690088235654300 /*sigma=0.01 challenging BENCHOP*/

const double x0 = log(S0/K);

double V(int);
double phi(int,int,double);
double coef(int);
std::complex<double> CharacteristicFunction(double);
double I1(double,double,int,int);
double I2(double,double,int,int);

double *ueval,pw2n,pw2n2;
int m;
std::complex<double> *f_hat; double trap_payoff;
int k_2,Jd,pwJd;
int Jp,pw2Jp,maxabsk;

using namespace std::complex_literals;

int main()
{
    double xi,sum,ext_inf,ext_sup,c1,c2;
    int i,k,min_k,max_k;
    FILE *out;
//    FILE *outpayoff;
    double *coefic;
    double call_shannon;
    clock_t inicio,parada;
    double exact,max_error;
    double a,Mm,A;
//    double D;

//    int cont;

    inicio=clock();
//for(cont=1;cont<=100;cont++){ 
    pw2n=pow(2,n);
    pw2n2=sqrt(pw2n);

    //Intervalo a priori y k mínima y máxima
    c1=x0+(r-0.5*pow(sigma,2))*T;
    c2=pow(sigma,2)*T;
    ext_inf=c1-L*sqrt(c2); printf("ext_inf=%lf\n",ext_inf);
    ext_sup=c1+L*sqrt(c2); printf("ext_sup=%lf\n",ext_sup);

    //ext_inf=-8.;
    //ext_sup=12.;

    min_k=ceil(pw2n*ext_inf); printf("min_k=%d\n",min_k);
    max_k=floor(pw2n*ext_sup); printf("max_k=%d\n",max_k);
    k_2=max_k;

    //J usado en la densidad
    a=std::max(std::fabs(ext_inf),std::fabs(ext_sup));
    Mm=std::max(std::fabs(pw2n*a-min_k),std::fabs(pw2n*a+max_k));
    Jd=ceil(log2(PI*Mm));

    printf("Jd>=%d\t\tJdnoredondeo=%lf\n",Jd,log2(PI*Mm));
    //Jd=10;
    pwJd=pow(2,Jd-1);

    /* Fourier transfom evaluation */
    ueval= static_cast<double*>(malloc(pwJd * sizeof(double)));
    f_hat= static_cast<std::complex<double>*>(malloc(pwJd * sizeof(std::complex<double>)));
    for(i=0;i<pwJd;i++)
    {
        ueval[i]=((2.*i+1.)*PI)/(2.*pwJd);
        f_hat[i]= CharacteristicFunction(pw2n * ueval[i]);
    }

    /*Coefficients*/
    coefic= static_cast<double*>(malloc((max_k - min_k + 1) * sizeof(double)));

    for(k=min_k;k<=max_k;k++) {
        coefic[k-min_k]=coef(k);
        //printf("%e\n",coefic[k-min_k]);
    }
    printf("f(k_1/2^m)=%.15e\tf(k_2/2^m)=%.15e\n",pw2n2*coefic[0],pw2n2*coefic[max_k-min_k]);

    //Payoff
    //Suponemos K1<=0 y por tanto k1 barra=0, de esta forma el maximo Nk=k2-k1
    maxabsk=abs(max_k-min_k);
    Jp=ceil(log2(PI*maxabsk));
    printf("Jp>=%d\t\tJpnoredondeo=%lf\n",Jp,log2(PI*maxabsk));
    //Jp=5;
    pw2Jp=pow(2,Jp-1);

    /*Call price*/
    call_shannon=0;

    for(k=min_k;k<=max_k;k++) {
        call_shannon+=coefic[k-min_k]*V(k);
        //printf("V(%d)=%lf\n",k,V(k));
    }
    call_shannon*=K*exp(-r*T)*pw2n2;

//} 
    parada=clock();

    printf("call_shannon=%.15lf\terror=%.2e\n",call_shannon,std::fabs(call_shannon-reference));
    printf("Tiempo ejecucion: %.3f\n",((double)(parada-inicio))/CLOCKS_PER_SEC);

    /*Area bajo la funcion de densidad=1?*/
    A=0.5*(coefic[0]+coefic[max_k-min_k]);
    for(k=min_k+1;k<max_k;k++) A+=coefic[k-min_k];
    A*=1./pw2n2;
    printf("A=%.15e\tA_error=%.2e\n",A,std::fabs(A-1.));

    /*Aproximacion a la funcion de densidad*/
    max_error=1.e-20;
    out=fopen("density_Shannon_viete_density.dat","w");
    for(xi=ext_inf+0.01;xi<ext_sup;xi+=0.01)
    {
        sum=0;
        for(k=min_k;k<=max_k;k++)
        {
            sum+=coefic[k-min_k]*phi(n,k,xi);
        }
        exact=(1./(sqrt(T)*sigma*sqrt(2*PI)))*exp(-pow((xi-log(S0/K)-(r-0.5*pow(sigma,2))*T),2)/(2*pow(sigma,2)*T));
        fprintf(out,"%lf %e %e %e\n",xi,sum,exact,std::fabs(sum-exact));
        if(std::fabs(sum-exact)>max_error) max_error=std::fabs(sum-exact);
    }
    fclose(out);
    //printf("max_error_density=%.2e\n",max_error);

    free(ueval);
    free(f_hat);
    free(coefic);

    return 0;
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