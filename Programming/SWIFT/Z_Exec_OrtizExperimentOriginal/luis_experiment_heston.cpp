/*Version 1: la f_hat se evalua 1 sola vez sobre la malla de puntos para todos los coeficientes */
/*La formula de Viete se usa para calcular la densidad*/
/*Aplicamos FFT a coeficientes densidad y a coeficientes pay-off*/

// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>

#define PI M_PI /* PI to machine precision, defined in math.h */
#define TWOPI (2.0*PI)
#define max(x,y) (x>y?x:y)

#define n 1

#define S0 100.
#define K 100.
#define r 0.
#define q 0.
#define T 45.

#define  mu r
#define lambda 1.5768
#define eta 0.5751
#define vmean 0.0398
#define v0 0.0175
#define rho -0.5711

//#define reference 5.007054273487142e+01 /*K=50 L=10 COS 50000 terms*/
//#define reference 5.785155700563814e+00 /*K=100 L=10 COS 50000 terms*/

//#define reference 3.887893512658775e+01 /* T=30 */
#define reference 4.691153171658516e+01 /* T=45 */

#define L 12.
#define x0 log(S0/K)

double V(int);
double phi(int,int,double);
double coef(int);
double complex chf(double);
double I1(double,double,int,int);
double I2(double,double,int,int);
void FFT_ForPDFCoefficients(double []);
void FFT_ForPayoffCoefficients(double []);

double *ueval,pw2n,pw2n2;
int m;
double complex *f_hat; double trap_payoff;
int k_2,Jd,pwJd;
int Jp,pw2Jp,maxabsk;
int min_k,max_k;

int main()
{
    double xi,sum,ext_inf,ext_sup,c1,c2;
    double p1,p2,p3,p4,p5;
    int i,k;
    FILE *out,*outpayoff;
    double *coefic;
    double call_shannon;
    clock_t inicio,parada;
    double max_error;
    double a,Mm,D;
    double A;
    double *payoffcoefic;

    int cont;

    inicio=clock();
    for(cont=1;cont<=1;cont++){
        pw2n=pow(2,n);
        pw2n2=sqrt(pw2n);

        //Intervalo
        c1=mu*T+(1.-exp(-lambda*T))*(vmean-v0)/(2.*lambda)-0.5*vmean*T;
        p1=eta*T*lambda*exp(-lambda*T)*(v0-vmean)*(8.*lambda*rho-4.*eta);
        p2=lambda*rho*eta*(1.-exp(-lambda*T))*(16.*vmean-8.*v0);
        p3=2.*vmean*lambda*T*(-4.*lambda*rho*eta+pow(eta,2)+4.*pow(lambda,2));
        p4=pow(eta,2)*((vmean-2.*v0)*exp(-2.*lambda*T)+vmean*(6.*exp(-lambda*T)-7.)+2.*v0);
        p5=8.*pow(lambda,2)*(v0-vmean)*(1.-exp(-lambda*T));
        c2=(1./(8.*pow(lambda,3)))*(p1+p2+p3+p4+p5);

        ext_inf=x0+c1-L*sqrt(fabs(c2));
        ext_sup=x0+c1+L*sqrt(fabs(c2));
        //printf("ext_inf=%lf\text_sup=%lf\n",ext_inf,ext_sup);

        //Rango de k's para la densidad
        min_k=ceil(pw2n*ext_inf); printf("min_k=%d\n",min_k);
        max_k=floor(pw2n*ext_sup); printf("max_k=%d\n",max_k);
        k_2=max_k;

        //J usado en la densidad
        a=max(fabs(ext_inf),fabs(ext_sup)); //printf("a=%lf\n",a);
        Mm=max(fabs(pw2n*a-min_k),fabs(pw2n*a+max_k)); //printf("Mm=%lf\n",Mm);
        Jd=ceil(log2(PI*Mm)); //printf("Jd=%d\n",Jd);

        //printf("Jd>=%d\t\tJdnoredondeo=%lf\n",Jd,log2(PI*Mm));
        //Jd=10;
        pwJd=pow(2,Jd-1); //printf("pwJd=%d\n",pwJd);

        /* Fourier transfom evaluation */
        ueval=malloc(pwJd*sizeof(double));
        //f_hat=malloc(pwJd*sizeof(double complex));

        for(i=0;i<pwJd;i++)
        {
            ueval[i]=((2.*i+1.)*PI)/(2.*pwJd);
            //f_hat[i]=chf(pw2n*ueval[i]);
        }

        /*Coefficients*/
        coefic=malloc((max_k-min_k+1)*sizeof(double));

        /*Coefficients without FFT*/
        //for(k=min_k;k<=max_k;k++) coefic[k-min_k]=coef(k);

        /*Coefficients with FFT*/
        FFT_ForPDFCoefficients(coefic);

        //Payoff
        //Suponemos K1<=0 y por tanto k1 barra=0, de esta forma el maximo Nk=k2-k1

        maxabsk=abs(max_k-min_k);
        Jp=ceil(log2(PI*maxabsk)); //printf("Jp=%d\n",Jp);
        pw2Jp=pow(2,Jp-1);

        /*Payoff with FFT*/
        payoffcoefic=malloc((max_k-min_k+1)*sizeof(double));

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

    printf("call_shannon=%.15lf\terror=%.2e\n",call_shannon,fabs(call_shannon-reference));
    printf("Tiempo ejecucion: %.3f\n",((double)(parada-inicio))/CLOCKS_PER_SEC);

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
    int j;
    double sum;

    sum=0;
    for(j=0;j<pwJd;j++) sum+=creal(f_hat[j]*cexp(k*ueval[j]*I));

    return(pw2n2*sum/(pwJd));
}

double complex chf(double w)
{
    double complex e1,e2,D,G;
    D=csqrt(cpow(lambda+rho*eta*w*I,2)+(cpow(w,2)-w*I)*cpow(eta,2));
    G=(lambda+rho*eta*w*I-D)/(lambda+rho*eta*w*I+D);
    e1=-mu*w*T*I+(v0/pow(eta,2))*((1.-cexp(-D*T))/(1.-G*cexp(-D*T)))*(lambda+rho*eta*w*I-D);
    e2=((lambda*vmean)/pow(eta,2))*(T*(lambda+rho*eta*w*I-D)-2.*clog((1.-G*cexp(-D*T))/(1.-G)));

    return(cexp(-w*x0*I)*cexp(e1)*cexp(e2));
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
    complex double sum1;

    int i;
    int Npoints = 2*pwJd;
    fftw_complex *in, *out1;
    fftw_plan plan1;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);                 //allocating memory
    out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npoints);               //allocating memory

    for(i=0;i<pwJd;i++)
        in[i]=chf(ueval[i]*pw2n);
    for(i=pwJd;i<Npoints;i++)
        in[i] = 0.;

    //Pels coeficients positius
    plan1 = fftw_plan_dft_1d(Npoints, in, out1, FFTW_BACKWARD, FFTW_ESTIMATE);       //Here we set which kind of transformation we want to perform
    fftw_execute(plan1);                                                             //Execution of FFT
    for(i=1;i<=max_k;i++)
    {
        sum1=cexp(M_PI*i*I/Npoints)*(out1[i]);
        FFT_CoefVector[-min_k+i]=pw2n2*(creal(sum1))/pwJd;
    }

    //Pels coeficients negatius
    plan1 = fftw_plan_dft_1d(Npoints, in, out1, FFTW_FORWARD, FFTW_ESTIMATE);        //Here we set which kind of transformation we want to perform
    fftw_execute(plan1);                                                              //Execution of FFT
    for(i=0;i<=-min_k;i++)
    {
        sum1=cexp(-M_PI*i*I/Npoints)*(out1[i]);
        FFT_CoefVector[-min_k-i]=pw2n2*(creal(sum1))/pwJd;
    }

    fftw_destroy_plan(plan1);                                                 //Destroy plan
    fftw_free(in);                                                            //Free memory
    fftw_free(out1);                                                          //Free memor

    return;

}

void FFT_ForPayoffCoefficients(double FFT_payoffcoefVector[])
{
    int i,j,k;
    double sup=k_2/pw2n; //Suponemos que estamos en el escenario de intervalo: [0,k2/2^m]
    double inf=0;
    double A,B,Co,sa,sb,ca,cb,I11,I12,I21,I22,ea,eb;
    double *d,*e;

    double res;

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
    free(out);                                                           //Free memor
    free(in2);                                                            //Free memory
    free(out2);

    return;
}