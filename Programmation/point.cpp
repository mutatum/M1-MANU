#include "point.hpp"
#include <cmath>
#include <iostream>

Point::Point(double x, double y) : M_x(x), M_y(y) {};
Point::Point(Point &p) : M_x(p.M_x), M_y(p.M_y) {};
Point::~Point() {};
double Point::X() const { return this->M_x; }
double Point::Y() const { return this->M_y; }
double &Point::X() { return this->M_x; }
double &Point::Y() { return this->M_y; }
void Point::zero() { this->M_x = this->M_y = 0; }
double Point::norm() const { return std::sqrt(M_x + M_y); }
Point &Point::operator=(const Point &other) {
  M_x = other.M_x;
  M_y = other.M_y;
  std::cout << "copy assignment" << std::endl;

  return *this;
}

Point &Point::operator=(Point &&other) noexcept {
  if (&other == this)
    return *this;
  M_x = other.M_x;
  M_y = other.M_y;
  other.M_x = 0.;
  other.M_y = 0.;
  std::cout << "move assignment" << std::endl;

  return *this;
}

Point &Point::operator=(double xx) {
  M_x = M_y = xx;
  return *this;
}

Point::operator double() const { return this->norm(); }

const Point operator-(const Point &p) {
    const Point ret(-p.M_x, -p.M_y);
    return ret;
}

std::ostream &operator<<(std::ostream &os, const Point &p) {
  os << "(x, y) = (" << p.M_x << ", " << p.M_y << ")";
  return os;
}