# pragma once
# include <iostream>


class Point {
private:
    static const unsigned long Ndim = 3;
    double coords[Ndim];

public:
    //Constructors
    Point(double aa=0);
    Point(double aa[Point::Ndim]);
    Point(const Point& P);

    // Operators
    double operator[](unsigned long i) const;
    double& operator[](unsigned long i);
    friend std::ostream& operator<<(std::ostream& o, const Point& p);
    void set(unsigned long i, const double& a);
    Point& operator=(const Point& P);
    
    friend Point operator-(const Point& P);
    friend Point operator+(const Point& P);
    friend Point operator+(const Point& Q, const Point& P);
    friend Point operator-(const Point& Q, const Point& P);
    Point operator-=(const Point& P);
    Point operator+=(const Point& P);
};