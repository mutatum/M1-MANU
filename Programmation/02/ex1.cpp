#include <iostream>
#include "point2d.hpp"


int main() {

    Point2d p=1;
    std::cout << p << std::endl;

    Point2d q{0,1};
    std::cout << q << std::endl;

    Point2d P,Q,W(1,2);

    std::cout << P << Q << W << std::endl;
    P = Q = W;
    std::cout << P << Q << W << std::endl;

    P = std::move(Q);
    std::cout << P << Q << std::endl;

    P = 1;
    std::cout << "+P and -P operators\n";
    std::cout << P << " | ";
    std::cout << -P << "\n";

    const Point2d T{1,2};

    std::cout << "+= and -= Operators\n";
    std::cout << T << " | P+=T: " << (P+=T) << '\n';
    P -= T;
    P -= T;
    std::cout << T << " | P-=T: " << P << '\n';

    std::cout << "+= and -= for P -= double \n";
    std::cout << "P += 3: " << (P+=3) << "\n\n";

    Q = Point2d{-1,10};
    std::cout << "P + Q operator: " << (P+Q) << "\n\n";

    std::cout << (P+1.) << (1.+P) << "\n\n";
    std::cout << (P-2.) << (2.-P) << "\n\n";

    return 0;
}