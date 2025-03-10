#include "integral.hpp"
#include <iostream>
#include <cmath>
#include <cassert>

// Define a test function
double test_function(double x) {
    return x * x;  // f(x) = x^2
}

void test_trapezoidal() {
    Integral integral(0.0, 1.0, test_function);
    double result = integral.trapezoidal(1000);
    double expected = 1.0 / 3.0; // Exact integral of x^2 from 0 to 1
    assert(std::abs(result - expected) < 1e-3);
    std::cout << "Trapezoidal method passed!" << std::endl;
}

void test_rectangular() {
    Integral integral(0.0, 1.0, test_function);
    double result = integral.rectangular(1000);
    double expected = 1.0 / 3.0;
    assert(std::abs(result - expected) < 1e-3);
    std::cout << "Rectangular method passed!" << std::endl;
}

int main() {
    test_trapezoidal();
    test_rectangular();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
