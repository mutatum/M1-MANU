#include "integral.hpp"

Integral::Integral(double a, double b, pfn f) {
    this->_lowbd = a;
    this->_upbd = b;
    this->func = f;
}

double Integral::lowbd() const {
    return this->_lowbd;
}

double Integral::uprbd() const {
    return this->_upbd;
}

void Integral::change_bounds(double a, double b) {
    this->_lowbd = a;
    this->_upbd =  b;
}

double Integral::rectangular(size_t n) const {
    double sum(0.0);
    double h((this->_upbd - this->_lowbd)/n);
    for (size_t i = 0; i<n; i++) {
        sum += h * this->func(this->_lowbd + i*h);
    }
    return sum;
}

double Integral::trapezoidal(size_t n) const {
    double sum(0.0);
    double h((this->_upbd - this->_lowbd)/n);
    for (size_t i = 0; i<n; i++) {
        sum += h * (this->func(this->_lowbd + (i+1.0)*h)+this->func(this->_lowbd + i*h)) / 2.0;
    }
    return sum;
}