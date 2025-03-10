#include <iostream>
#include "point.hpp"

int main() {
    std::cout << "=== Testing Point<2> ===\n";
    Point<2> p2a(3.0);
    double arr2[2] = {1.5, 2.5};
    Point<2> p2b(arr2);
    std::cout << "p2a: " << p2a << "\n";  // Expected: (3, 3)
    std::cout << "p2b: " << p2b << "\n";  // Expected: (1.5, 2.5)
    std::cout << "p2a + p2b: " << (p2a + p2b) << "\n";  // Expected: (4.5, 5.5)

    std::cout << "\n=== Testing Point<3> ===\n";
    Point<3> p3a(1.0);
    double arr3[3] = {2.0, 3.0, 4.0};
    Point<3> p3b(arr3);
    std::cout << "p3a: " << p3a << "\n";  // Expected: (1, 1, 1)
    std::cout << "p3b: " << p3b << "\n";  // Expected: (2, 3, 4)
    std::cout << "p3a + p3b: " << (p3a + p3b) << "\n";  // Expected: (3, 4, 5)
    std::cout << "-p3b: " << (-p3b) << "\n";  // Expected: (-2, -3, -4)

    std::cout << "\n=== Testing Point<4> ===\n";
    Point<4> p4a(0.5);
    double arr4[4] = {1.0, 2.0, 3.0, 4.0};
    Point<4> p4b(arr4);
    std::cout << "p4a: " << p4a << "\n";  // Expected: (0.5, 0.5, 0.5, 0.5)
    std::cout << "p4b: " << p4b << "\n";  // Expected: (1, 2, 3, 4)
    std::cout << "p4b - p4a: " << (p4b - p4a) << "\n";  // Expected: (0.5, 1.5, 2.5, 3.5)

    std::cout << "\n=== Testing Point<5> ===\n";
    Point<5> p5a(2.0);
    double arr5[5] = {5.0, 4.0, 3.0, 2.0, 1.0};
    Point<5> p5b(arr5);
    std::cout << "p5a: " << p5a << "\n";  // Expected: (2, 2, 2, 2, 2)
    std::cout << "p5b: " << p5b << "\n";  // Expected: (5, 4, 3, 2, 1)
    std::cout << "p5a += p5b: " << (p5a += p5b) << "\n";  // Expected: (7, 6, 5, 4, 3)

    return 0;
}
