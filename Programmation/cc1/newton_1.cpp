#include <iostream>
#include "newtons.hpp"

// f(x) = exp(x) + x**3 - 
// initial x0 = 0


double newton_1(double x0) {
    double iterations[100] = {x0};

    for (int i = 1; i < 100; i++) {
        iterations[i] = iterations[i-1] - f(iterations[i-1])/df(iterations[i-1]);
        std::cout << "i: " << i << ": (xi, f(xi)) = (" << iterations[i] << ", " << f(iterations[i]) << ")\n";
    }
    return iterations[99];
}