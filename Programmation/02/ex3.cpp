#include "dynamic_vector.hpp"
#include <iostream>

int main() {

    // Create dynamic_vector from array
    double arr[3] = {1.0, 2.0, 3.0};
    dynamic_vector v1(3, arr);
    
    // Create dynamic_vector with constant value
    dynamic_vector v2(3, 10.0);
    
    std::cout << "Initial vectors:" << std::endl;
    std::cout << "v1 = " << v1 << std::endl;
    std::cout << "v2 = " << v2 << std::endl;
    
    // Test assignment from double (sets all entries)
    v1 = 5.0;
    std::cout << "\nAfter v1 = 5.0:" << std::endl;
    std::cout << "v1 = " << v1 << std::endl;
    
    // Test multiplication
    dynamic_vector v3 = 2 * v1;
    dynamic_vector v4 = v1 * 3;
    std::cout << "\nAfter multiplication:" << std::endl;
    std::cout << "v3 = 2 * v1 = " << v3 << std::endl;
    std::cout << "v4 = v1 * 3 = " << v4 << std::endl;
    
    // Test binary + and -
    dynamic_vector v5 = v1 + v2;
    dynamic_vector v6 = v2 - v1;
    std::cout << "\nAfter binary operations:" << std::endl;
    std::cout << "v5 = v1 + v2 = " << v5 << std::endl;
    std::cout << "v6 = v2 - v1 = " << v6 << std::endl;
    
    // Test compound assignments
    v5 += v1;
    std::cout << "\nAfter v5 += v1:" << std::endl;
    std::cout << "v5 = " << v5 << std::endl;
    
    v6 -= v2;
    std::cout << "\nAfter v6 -= v2:" << std::endl;
    std::cout << "v6 = " << v6 << std::endl;
    
    // Test unary operators
    dynamic_vector v7 = +v1;
    dynamic_vector v8 = -v1;
    std::cout << "\nUnary operations:" << std::endl;
    std::cout << "v7 = +v1 = " << v7 << std::endl;
    std::cout << "v8 = -v1 = " << v8 << std::endl;
    return 0;
}