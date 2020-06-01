#include <iostream>
#include <SWIFT/known_distribution_contract_combinations.h>
#include "fstream"
#include "vector"
#include "SWIFT/distributions.h"

// market parameters: you may change the number of observations by modifying the size of T and K
struct MarketParameters{
    double S;
    double r;
    double q = 0;
    double T[40];
    double K[40];
};

int main()
{
//    int m = 5;  // # of parameters
    int n_observations = 40; // # of observations (consistent with the struct MarketParameters)

    MarketParameters market_parameters;


    //// INIT MARKET PARAMETERS ////
    // array of strikes
    std::vector<double> K_over_S = {
        0.9371, 0.8603, 0.8112, 0.7760, 0.7470, 0.7216, 0.6699, 0.6137,
        0.9956, 0.9868, 0.9728, 0.9588, 0.9464, 0.9358, 0.9175, 0.9025,
        1.0427, 1.0463, 1.0499, 1.0530, 1.0562, 1.0593, 1.0663, 1.0766,
        1.2287, 1.2399, 1.2485, 1.2659, 1.2646, 1.2715, 1.2859, 1.3046,
        1.3939, 1.4102, 1.4291, 1.4456, 1.4603, 1.4736, 1.5005, 1.5328};

    // array of expiries
    std::vector<double> expiries = {0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
                         0.119047619047619	, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714	, 1.07142857142857, 1.42857142857143	,
                         0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
                         0.119047619047619, 0.238095238095238, 0.357142857142857, 0.476190476190476	, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
                         0.119047619047619, 0.238095238095238	, 0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143};


    // spot and interest rate
    market_parameters.S = 1.0;
    market_parameters.r = 0.02;

    // strikes and expiries
    for (int j = 0; j < n_observations; ++j) {
        market_parameters.K[j] = K_over_S[j];
        market_parameters.T[j] = expiries[j];
    }

    //// END INIT MARKET PARAMETERS

    // you may set up your optimal model parameters here:
    HestonParameters optimal_parameters{3.0, 0.10, 0.25, -0.8, 0.08};

    // compute the market_parameters observations with pstar
    std::vector<double> x;
    x.reserve(K_over_S.size());
    for (std::size_t i = 0; i < K_over_S.size(); ++i)
    {
        double current_value = GetHestonEuropeanPrice(optimal_parameters, market_parameters.S, market_parameters.K[i], market_parameters.r, market_parameters.T[i], market_parameters.q);
        x.push_back(current_value);
        std::cout << "K " << market_parameters.K[i] << " T " << market_parameters.T[i] << " value " << current_value << std::endl;
    }

//    // >>> Enter calibrating routine >>>
//    double start_s = clock();
//
//    // algorithm parameters
//    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
//    opts[0]=LM_INIT_MU;
//    // stopping thresholds for
//    opts[1]=1E-10;       // ||J^T e||_inf
//    opts[2]=1E-10;       // ||Dp||_2
//    opts[3]=1E-10;       // ||e||_2
//    opts[4]= LM_DIFF_DELTA; // finite difference if used
//
//    // you may set up your initial point here:
//    double p[5];
//    p[0] = 1.2000;
//    p[1] = 0.20000;
//    p[2] = 0.3000;
//    p[3] = -0.6000;
//    p[4] = 0.2000;
//
//    cout << "\r-------- -------- -------- Heston Model Calibrator -------- -------- --------"<<endl;
//    cout << "Parameters:" << "\t         kappa"<<"\t     vinf"<< "\t       vov"<< "\t      rho" << "\t     v0"<<endl;
//    cout << "\r Initial point:" << "\t"  << scientific << setprecision(8) << p[0]<< "\t" << p[1]<< "\t"<< p[2]<< "\t"<< p[3]<< "\t"<< p[4] << endl;
//    // Calibrate using analytical gradient
//    dlevmar_der(GetHestonPrice, JacHes, p, x, m, n_observations, 100, opts, info, NULL, NULL, (void*) &market_parameters);
//
//    double stop_s = clock();
//
//    cout << "Optimum found:" << scientific << setprecision(8) << "\t"<< p[0]<< "\t" << p[1]<< "\t"<< p[2]<< "\t"<< p[3]<< "\t"<< p[4] << endl;
//    cout << "Real optimum:" << "\t" << pstar[0]<<"\t"<< pstar[1]<< "\t"<< pstar[2]<< "\t"<< pstar[3]<< "\t"<< pstar[4] << endl;
//
//    if (int(info[6]) == 6) {
//        cout << "\r Solved: stopped by small ||e||_2 = "<< info[1] << " < " << opts[3]<< endl;
//    } else if (int(info[6]) == 1) {
//        cout << "\r Solved: stopped by small gradient J^T e = " << info[2] << " < " << opts[1]<< endl;
//    } else if (int(info[6]) == 2) {
//        cout << "\r Solved: stopped by small change Dp = " << info[3] << " < " << opts[2]<< endl;
//    } else if (int(info[6]) == 3) {
//        cout << "\r Unsolved: stopped by itmax " << endl;
//    } else if (int(info[6]) == 4) {
//        cout << "\r Unsolved: singular matrix. Restart from current p with increased mu"<< endl;
//    } else if (int(info[6]) == 5) {
//        cout << "\r Unsolved: no further error reduction is possible. Restart with increased mu"<< endl;
//    } else if (int(info[6]) == 7) {
//        cout << "\r Unsolved: stopped by invalid values, user error"<< endl;
//    }
//
//    cout << "\r-------- -------- -------- Computational cost -------- -------- --------"<<endl;
//    cout << "\r          Time cost: "<< double(stop_s - start_s) / CLOCKS_PER_SEC << " seconds "<<endl;
//    cout << "         Iterations: " << int(info[5]) << endl;
//    cout << "         pv  Evalue: " << int(info[7]) << endl;
//    cout << "         Jac Evalue: "<< int(info[8]) << endl;
//    cout << "# of lin sys solved: " << int(info[9])<< endl; //The attempts to reduce error
//    cout << "\r-------- -------- -------- Residuals -------- -------- --------"<<endl;
//    cout << " \r            ||e0||_2: " << info[0] << endl;
//    cout << "           ||e*||_2: " << info[1]<<endl;
//    cout << "          ||J'e||_inf: " << info[2]<<endl;
//    cout << "           ||Dp||_2: " << info[3]<<endl;
//
//    return 0;
}

/////// RESULTS ////////
// K 0.9371 T 0.119048 value 0.0803314
// K 0.8603 T 0.238095 value 0.154888
// K 0.8112 T 0.357143 value 0.205334
// K 0.776 T 0.47619 value 0.242596
// K 0.747 T 0.595238 value 0.273595
// K 0.7216 T 0.714286 value 0.300806
// K 0.6699 T 1.07143 value 0.358526
// K 0.6137 T 1.42857 value 0.417245
// K 0.9956 T 0.119048 value 0.0429609
// K 0.9868 T 0.238095 value 0.0657904
// K 0.9728 T 0.357143 value 0.0878723
// K 0.9588 T 0.47619 value 0.108093
// K 0.9464 T 0.595238 value 0.126123
// K 0.9358 T 0.714286 value 0.142136
// K 0.9175 T 1.07143 value 0.178126
// K 0.9025 T 1.42857 value 0.208112
// K 1.0427 T 0.119048 value 0.0225191
// K 1.0463 T 0.238095 value 0.0382272
// K 1.0499 T 0.357143 value 0.0507916
// K 1.053 T 0.47619 value 0.0618001
// K 1.0562 T 0.595238 value 0.0715895
// K 1.0593 T 0.714286 value 0.0805466
// K 1.0663 T 1.07143 value 0.104719
// K 1.0766 T 1.42857 value 0.123913
// K 1.2287 T 0.119048 value 0.000331382
// K 1.2399 T 0.238095 value 0.00288361
// K 1.2485 T 0.357143 value 0.00721851
// K 1.2659 T 0.47619 value 0.0113151
// K 1.2646 T 0.595238 value 0.0177755
// K 1.2715 T 0.714286 value 0.0234064
// K 1.2859 T 1.07143 value 0.040997
// K 1.3046 T 1.42857 value 0.0566055
// K 1.3939 T 0.119048 value 4.31064e-07
// K 1.4102 T 0.238095 value 8.3199e-05
// K 1.4291 T 0.357143 value 0.000554006
// K 1.4456 T 0.47619 value 0.00159113
// K 1.4603 T 0.595238 value 0.00317842
// K 1.4736 T 0.714286 value 0.00522975
// K 1.5005 T 1.07143 value 0.0138154
// K 1.5328 T 1.42857 value 0.0232216