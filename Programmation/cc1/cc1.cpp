#include "newtons.hpp"
#include "euler_method.hpp"
#include <iostream>
#include <cmath>

// f(x) = exp(x) + x**3 - 5
// initial x0 = 0

double f(double x) {
    return std::exp(x) + x*x*x - 5;
}

double df(double x) {
    return std::exp(x) + 3 * x * x;
}

double function(double t, double x) {
    return ((1-std::exp(t))*x)/(1+std::exp(t));
}

int main() {

    std::cout << "newton_1: \n";
    float x0 = 0.0;
    newton_1(x0);
    std::cout << "newton_2: \n";
    newton_2(x0);
    std::cout << "newton_3: \n";
    newton_3(x0);
    std::cout << "euler: \n";

    double *solution;
    solution = euler_method_space(0., 3., 2., &function, 100);
    double err = 0.0;
    double h = 2./100.;
    double t=0.;
    for (int i=0; i<100; i++) {
        t += h;
        err += std::abs(solution[i]- (12*std::exp(t))/ std::pow((1+std::exp(t)),2));
        std::cout << "t="<<t<< ": numerical: " << solution[i] << " real: " <<(12*std::exp(t))/ std::pow((1+std::exp(t)),2) << '\n';
    }
    std::cout << "Error t=[0,2], N=100: " <<  err << std::endl;


    return 0;

}