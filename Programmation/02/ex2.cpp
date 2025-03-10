#include "point.hpp"
#include <iostream>

int main(void) {

    std::cout << "Dimension: 3" << '\n';
    Point P{};
    std::cout << P << '\n';
    Point Q{-4};
    std::cout << Q << '\n';
    double arr[3] = {1,2,3};
    Point W{arr};
    std::cout << W << '\n';

    std::cout << "P = " << P;
    P.set(0,10.);
    std::cout << "P.set(0,10), P= " << P << '\n';
    P = Q;
    std::cout << "operator=: P = Q: " << P << " Q: " << Q << "\n";

    P = 1;
    P += Q;
    std::cout << P << "\n";
    P += -Q;
    std::cout << P << "\n";
    P += -Q;
    std::cout << P << "\n";
    return 0;
}