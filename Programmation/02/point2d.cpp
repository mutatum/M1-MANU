#include "point2d.hpp"
#include <cmath>
#include <iostream>

Point2d::Point2d(double x, double y) : M_x(x), M_y(y) {};
Point2d::Point2d(Point2d &p) : M_x(p.M_x), M_y(p.M_y) {};
Point2d::~Point2d() {};
double Point2d::X() const { return this->M_x; }
double Point2d::Y() const { return this->M_y; }
double &Point2d::X() { return this->M_x; }
double &Point2d::Y() { return this->M_y; }
void Point2d::zero() { this->M_x = this->M_y = 0; }
double Point2d::norm() const { return std::sqrt(M_x + M_y); }
Point2d &Point2d::operator=(const Point2d &other) {
  M_x = other.M_x;
  M_y = other.M_y;
  std::cout << "copy assignment" << std::endl;

  return *this;
}

Point2d &Point2d::operator=(Point2d &&other) noexcept {
  if (&other == this)
    return *this;
  M_x = other.M_x;
  M_y = other.M_y;
  other.M_x = 0.;
  other.M_y = 0.;
  std::cout << "move assignment" << std::endl;

  return *this;
}

Point2d &Point2d::operator=(double xx) {
  M_x = M_y = xx;
  return *this;
}

Point2d::operator double() const { return this->norm(); }

const Point2d Point2d::operator-(const Point2d &p) {
    const Point2d ret(-p.M_x, -p.M_y);
    return ret;
}

const Point2d Point2d::operator+(const Point2d &p) {
    const Point2d ret(p.M_x, p.M_y);
    return ret;
}

const Point2d& Point2d::operator+=(const Point2d& p) {
  this->M_x += p.M_x, this->M_y += p.M_y;
  return *this;
}

const Point2d& Point2d::operator-=(const Point2d& p) {
  this->M_x -= p.M_x, this->M_y -= p.M_y;
  return *this;
}

const Point2d& Point2d::operator+=(const double aa) {
  this->M_x += aa; this->M_y += aa;
  return *this;
}

const Point2d& Point2d::operator-=(const double aa) {
  this->M_x += aa; this->M_y += aa;
  return *this;
}

const Point2d operator+(const Point2d& p, const Point2d& q) {
  return Point2d{p.M_x+q.M_x, p.M_y + q.M_y};
}

const Point2d operator-(const Point2d& p, const Point2d& q) {
  return Point2d{p.M_x-q.M_x, p.M_y - q.M_y};
}

const Point2d operator+(const Point2d &p, double aa) {
  return Point2d(p.M_x + aa, p.M_y + aa);
}
const Point2d operator+(double aa, const Point2d &p) {
  return Point2d(p.M_x + aa, p.M_y + aa);
}
const Point2d operator-(const Point2d &p, double aa) {
  return Point2d(p.M_x - aa, p.M_y - aa);
}
const Point2d operator-(double aa, const Point2d &p) {
  return Point2d(p.M_x - aa, p.M_y - aa);
}

std::ostream &operator<<(std::ostream &os, const Point2d &p) {
  os << "(x, y) = (" << p.M_x << ", " << p.M_y << ")";
  return os;
}