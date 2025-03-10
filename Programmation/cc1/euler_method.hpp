# pragma once
# include <cstddef>


double euler_method(double t_ini, double i_sol, double t_end, double (*sfn) (double t, double x), size_t N);
double *euler_method_space(double t_ini, double i_sol, double t_end, double (*sfn) (double t, double x), size_t N);