#pragma once
#include <iostream>

class Point {
public:
  Point(double x = .0, double y = .0);
  Point(Point &p);
  ~Point();
  double X() const;
  double &X();
  double Y() const;
  double &Y();
  void zero();
  double norm() const;
  friend std::ostream &operator<<(std::ostream &os, const Point &p);
  Point &operator=(const Point &other);
  Point &operator=(Point &&other) noexcept;
  Point &operator=(double xx);
  operator double() const;
  friend const Point operator-(const Point& p);

private:
  double M_x;
  double M_y;
};