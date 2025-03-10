#include "euler_method.hpp"

double euler_method(double t_ini, double i_sol, double t_end, double (*sfn) (double t, double x), size_t N) {

    double roll[2] = {i_sol};
    bool ind = false;

    const double h = (t_end - t_ini) / N;
    for (auto i = 0; i<N; i++) {
        roll[!ind] = roll[ind] + h * sfn(t_ini + i * h, roll[ind]);
        ind=!ind;
    }
    return roll[!ind];
}

double *euler_method_space(double t_ini, double i_sol, double t_end, double (*sfn) (double t, double x), size_t N) {

    double *hist = new double[N];

    hist[0] = i_sol;
    const double h = (t_end - t_ini) / N;
    for (auto i = 0; i<N; i++) {
        hist[i+1] = hist[i] + h * sfn(t_ini + i * h, hist[i]);
    }
    return hist;
}