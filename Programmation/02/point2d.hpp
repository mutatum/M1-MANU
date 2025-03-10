#pragma once
#include <iostream>

class Point2d {
public:
  Point2d(double x = .0, double y = .0);
  Point2d(Point2d &p);
  ~Point2d();
  double X() const;
  double &X();
  double Y() const;
  double &Y();
  void zero();
  double norm() const;
  friend std::ostream &operator<<(std::ostream &os, const Point2d &p);
  Point2d &operator=(const Point2d &other);
  Point2d &operator=(Point2d &&other) noexcept;
  Point2d &operator=(double xx);
  operator double() const;
  const Point2d operator-(const Point2d& p);
  const Point2d operator+(const Point2d& p);
  const Point2d& operator+=(const Point2d& p);
  const Point2d& operator-=(const Point2d& p);
  const Point2d& operator+=(const double aa);
  const Point2d& operator-=(const double aa);
  friend const Point2d operator+(const Point2d& p, const Point2d& q);
  friend const Point2d operator-(const Point2d& p, const Point2d& q);
  friend const Point2d operator-(const Point2d& p, const double aa);
  friend const Point2d operator+(const Point2d& p, const double aa);
  friend const Point2d operator-(const double aa, const Point2d& p);
  friend const Point2d operator+(const double aa, const Point2d& p);

private:
  double M_x;
  double M_y;
};