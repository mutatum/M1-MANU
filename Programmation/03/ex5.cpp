# include <iostream>
# include "binomial.hpp"

const size_t M=11;
const Binomial<M> triangle;

int main() {
    std::cout << triangle(10,6) << std::endl;
    triangle.print();
    return 0;
}