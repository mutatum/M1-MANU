# pragma once
# include <cstddef>

typedef double (*pfn)(double);
class Integral {
private:

    double _lowbd;
    double _upbd;
    pfn func;

public:

    Integral(double a, double b, pfn f);
    double lowbd() const;
    double uprbd() const;
    void change_bounds(double a, double b);
    double trapezoidal(size_t n) const;
    double rectangular(size_t n) const;

};
