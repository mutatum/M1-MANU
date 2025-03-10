# pragma once
#include <cstddef>
# include <iostream>


class dynamic_vector {

private:
    size_t ndim=0;
    double* coord;

public:
    dynamic_vector(size_t N, double coord = 0.0);
    dynamic_vector(size_t N, double coord[]);
    dynamic_vector(const dynamic_vector &V);
    ~dynamic_vector();
    void set(size_t i, const double& a);
    double operator[](size_t i) const;
    double& operator[](size_t i);
    dynamic_vector& operator=(const dynamic_vector& V);
    friend dynamic_vector operator*(const double &a, const dynamic_vector &V);
    friend dynamic_vector operator*(const dynamic_vector &V, const double &a);
    const dynamic_vector operator=(const double& i);
    friend dynamic_vector operator+(const dynamic_vector& V, const dynamic_vector& W);
    friend dynamic_vector operator-(const dynamic_vector& V, const dynamic_vector& W);
    dynamic_vector operator+() const;
    dynamic_vector operator-() const;
    dynamic_vector& operator+=(const dynamic_vector &V);
    dynamic_vector& operator-=(const dynamic_vector &V);
    friend std::ostream& operator<<(std::ostream& o, const dynamic_vector &V);
};