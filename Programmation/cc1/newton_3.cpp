#include <iostream>
#include "newtons.hpp"

// f(x) = exp(x) + x**3 - 
// initial x0 = 0

double newton_3(double x0) {

    double roll[2] = {x0,0.};
    bool ind = true;
    size_t i = 0;
    const double eps = 1e-7;

    do {
        roll[ind] = roll[!ind] - f(roll[!ind])/df(roll[!ind]);
        std::cout << "i: " << i << ": (xi, f(xi)) = (" << roll[ind] << ", " << f(roll[ind]) << ")\n";
        ind = !ind;
        i++;
    } while (std::abs(roll[0]-roll[1]) > eps);


    return roll[!ind];
}