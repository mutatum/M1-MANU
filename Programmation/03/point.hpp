#pragma once
#include <iostream>

// Template Class Definition
template<size_t N>
class Point {
protected:
    double coords[N];

public:
    // Constructors
    Point(double aa = 0);
    Point(const double* arr);
    Point(const Point& P);

    // Operators
    double operator[](unsigned long i) const;
    double& operator[](unsigned long i);

    void set(unsigned long i, const double& a);
    Point& operator=(const Point& P);
    
    // Arithmetic Operators
    Point operator-() const;
    Point operator+(const Point& P) const;
    Point operator-(const Point& P) const;
    Point& operator+=(const Point& P);
    Point& operator-=(const Point& P);

    // Friend Functions
    template<size_t M>
    friend std::ostream& operator<<(std::ostream& o, const Point<M>& p);
};

// ===========================================

template<size_t N>
Point<N>::Point(double aa) {
    for (size_t i = 0; i < N; i++) {
        coords[i] = aa;
    }
}

template<size_t N>
Point<N>::Point(const double* arr) {
    for (size_t i = 0; i < N; i++) {
        coords[i] = arr[i];
    }
}

template<size_t N>
Point<N>::Point(const Point<N>& P) {
    for (size_t i = 0; i < N; i++) {
        coords[i] = P.coords[i];
    }
}

template<size_t N>
double Point<N>::operator[](unsigned long i) const {
    return coords[i];
}

template<size_t N>
double& Point<N>::operator[](unsigned long i) {
    return coords[i];
}

template<size_t N>
void Point<N>::set(unsigned long i, const double& a) {
    coords[i] = a;
}

template<size_t N>
Point<N>& Point<N>::operator=(const Point<N>& P) {
    if (this != &P) {
        for (size_t i = 0; i < N; i++) {
            coords[i] = P.coords[i];
        }
    }
    return *this;
}

template<size_t N>
Point<N> Point<N>::operator-() const {
    Point<N> V(*this);
    for (size_t i = 0; i < N; i++) {
        V.coords[i] *= -1.0;
    }
    return V;
}

template<size_t N>
Point<N> Point<N>::operator+(const Point<N>& P) const {
    Point<N> V(*this);
    for (size_t i = 0; i < N; i++) {
        V.coords[i] += P.coords[i];
    }
    return V;
}

template<size_t N>
Point<N> Point<N>::operator-(const Point<N>& P) const {
    Point<N> V(*this);
    for (size_t i = 0; i < N; i++) {
        V.coords[i] -= P.coords[i];
    }
    return V;
}

template<size_t N>
Point<N>& Point<N>::operator+=(const Point<N>& P) {
    for (size_t i = 0; i < N; i++) {
        coords[i] += P.coords[i];
    }
    return *this;
}

template<size_t N>
Point<N>& Point<N>::operator-=(const Point<N>& P) {
    for (size_t i = 0; i < N; i++) {
        coords[i] -= P.coords[i];
    }
    return *this;
}

template<size_t N>
std::ostream& operator<<(std::ostream& o, const Point<N>& p) {
    o << "Point: (";
    for (size_t i = 0; i < N; i++) {
        o << p.coords[i];
        if (i < N - 1) { o << ", "; }
    }
    o << ")";
    return o;
}
