#include <iostream>
#include "newtons.hpp"

// f(x) = exp(x) + x**3 - 
// initial x0 = 0

double newton_2(double x0) {

    double roll[2] = {x0,0.};
    bool ind = true;

    for (int i = 1; i < 100; i++) {
        roll[ind] = roll[!ind] - f(roll[!ind])/df(roll[!ind]);
        std::cout << "i: " << i << ": (xi, f(xi)) = (" << roll[ind] << ", " << f(roll[ind]) << ")\n";
        ind = !ind;
    }

    return roll[!ind];
}