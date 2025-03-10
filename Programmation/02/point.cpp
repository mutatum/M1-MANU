#include "point.hpp"

// Constructors
Point::Point(double aa) {
  for (unsigned long i = 0; i < Point::Ndim; i++) {
    this->coords[i] = aa;
  }
}

Point::Point(double arr[Point::Ndim]) {
  for (unsigned long i = 0; i < Point::Ndim; i++) {
    this->coords[i] = arr[i];
  }
}

Point::Point(const Point& P) {
  for (unsigned long i = 0; i < Point::Ndim; i++) {
    this->coords[i] = P[i];
  }
}

// Operators
double Point::operator[](unsigned long i) const { return this->coords[i]; }

double &Point::operator[](unsigned long i) { return this->coords[i]; }

std::ostream &operator<<(std::ostream &o, const Point& P) {
  o << "Point: (";
  for (unsigned long i=0; i<Point::Ndim; i++) {
    o << P[i];
    if (i < Point::Ndim-1) {o << ", ";}
  }
  o << ")";
  return o;
}

void Point::set(unsigned long i, const double& a) {
    this->coords[i] = a;
}

Point& Point::operator=(const Point& P) {
  for (unsigned long i=0; i<Point::Ndim; i++) {
    (*this)[i] = P[i];
  }
  return *this;
}

Point operator+(const Point& P) {
  Point V(P);
  return V;
}

Point operator-(const Point& P) {
  Point V(P);

  for (unsigned long i = 0; i < Point::Ndim; i++) {
    V[i] *= -1.;
  }
  return V;
}

Point operator-(const Point& Q, const Point& P) {
  Point V(Q);

  for (unsigned long i = 0; i < Point::Ndim; i++) {
    V[i] -= P[i];
  }
  return V;
}

Point operator+(const Point& Q, const Point& P) {
  Point V(Q);

  for (unsigned long i = 0; i < Point::Ndim; i++) {
    V[i] += P[i];
  }
  return V;
}

Point Point::operator+=(const Point& P) {

  for (unsigned long i = 0; i < Point::Ndim; i++) {
    (*this)[i] -= P[i];
  }
  return *this;
}

Point Point::operator-=(const Point& P) {
  *this = *this -P;
  return *this;
}