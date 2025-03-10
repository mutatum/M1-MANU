#include <iostream>
#include <cmath>
#include <ranges>
#include <vector>

template<typename T>
T square(T x) { return x*x;}

template<typename T>
T sum(T* array, size_t n) {
    T sum = T(0.0);
    for (size_t i = 0; i<n; i++) {
        sum += array[i];
    }
    return sum;
}

template<typename T>
auto lagrange(T* vx, T*vy, const size_t n){

    return [=] (T x_in) -> T {
        T out = T(0);
        for (size_t i=0; i<n; i++) {
            T prod = T(1);
            for (size_t j=0; j<n; j++) { // Construction Li
                if (i==j) continue;
                prod *= (x_in-vx[j])/(vx[i]-vx[j]);
            }
            out += vy[i] * prod; // yi * Li
        }
        return out;
    };
}


int main() {

    std::cout << square(2.) << "\n";
    std::cout << square(2) << "\n";

    double array[5] = {1,2,3,4,5};
    double chars[5] = {'a','b', 'c', 'd', 'e'};
    std::cout << sum(array, 5) << "\n";
    std::cout << sum(chars, 5) << "\n";

    std::cout << "Lagrange: \n";
    {
        const size_t N = 4;
        double vx[N] = {0.};
        for (size_t i = 0; i<N; i++) vx[i]=(double)i;
        std::for_each(vx, vx+N, [](double &x){x = 1.0+x/4.0;});
        std::cout << vx[2]<<'\n';
        double vy[N];
        std::transform(vx, vx+N, vy, [](double x){return std::exp(x);});
        std::cout << vy[0]<<'\n';
        auto interpolated_exp = lagrange(vx, vy, N);
        std::cout << "exp(1.4)=" << std::exp(1.4) <<'\n';
        std::cout << "Lagrange_exp(1.4)=" << interpolated_exp(1.4) <<'\n';
    }

    return 0;
}