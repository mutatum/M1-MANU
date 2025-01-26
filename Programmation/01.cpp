#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;
typedef unsigned long long ull;

// Exercise 1
template <typename T> T power(T x, size_t a) {
  T res(1);
  while (a > 0) {
    res *= x;
    --a;
  }
  return res;
}

// Exercise 2
ull factorial(size_t n) {
  static vector<ull> factorials = {1, 1};

  if (n < factorials.size()) {
    return factorials[n];
  }
  for (auto i = factorials.size(); i <= n; i++) {
    factorials.push_back(static_cast<ull>(i) * factorials[i - 1]);
  }
  return factorials[n];
}

ull binomial(size_t n, size_t k) {
  return factorial(n) / (factorial(k) * factorial(n - k));
}

template <typename T> class Matrix {
  size_t rows, cols;
  vector<T> data;

public:
  Matrix(size_t r, size_t c, T initialValue = T())
      : rows(r), cols(c), data(r * c, initialValue) {}

  T &operator()(size_t i, size_t j) { return data[i * rows + j]; }
  T operator()(size_t i, size_t j) const { return data[i * rows + j]; }

  int numRows() const { return rows; }
  int numCols() const { return cols; }

  void print() const {
    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < cols; j++) {
        cout << (*this)(i, j) << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
};

// Exercise 3
void print_pascal_triangle_fact(size_t n) {

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j <= i; j++) {
      cout << factorial(i) / (factorial(i - j) * factorial(j)) << " ";
    }
    cout << endl;
  }
}

void print_pascal_triangle_sum(size_t n) {
  Matrix<ull> triangle(n, n, 0);

  for (size_t i = 0; i <= n; i++) {
    for (size_t j = 0; j < i; j++) {
      if (j == 0 || j == i - 1) {
        triangle(n - i + j, j) = 1;
      } else {
        triangle(n - i + j, j) =
            triangle(n - i + j, j - 1) + triangle(n - i + j + 1, j);
      }
    }
  }
  triangle.print();
}

// Exercise 4
size_t calculate_m(double x) {
  size_t m = 0;
  while (abs(x) >= .5) {
    x /= 2;
    m++;
  }

  return m;
}

double expTaylor(double x) {
  double N(10);

  size_t m = calculate_m(x);
  x /= power(2.0, m);
  double e = 1 + x / N;
  for (size_t t = 0; t < 10; t++) {
    e *= x / (N - static_cast<double>(t));
    e++;
  }
  return power(e, power(2UL, m));
}

// Exercise 5
double Pade(double x, size_t N) {

  double P = 0.0;
  for (long n = static_cast<long>(N); n >= 0; --n) {
    P = static_cast<double>(binomial(N, n)) /
            (static_cast<double>(factorial(n)) *
             static_cast<double>(binomial(2 * N, n))) +
        x * P;
  }

  return P;
}

double expPade(double x) {

  size_t m = calculate_m(x);
  x /= power(2.0, m);

  return power((Pade(x, 5) / Pade(-x, 5)), power(2UL, m));
}

// Exercise 6
void CalculateRealAndImaginary(double r, double theta, double &pReal,
                               double &pImaginary) {

  pReal = r * cos(theta);
  pImaginary = r * sin(theta);
}

// Exercise 7
std::pair<double, double> quadraticRoot(const double &a, const double &b,
                                        const double &c) {

  double delta{b * b - 4 * a * c};

  if (delta < 0) {
    cerr << "Delta is negative: " << delta << endl;
    return std::make_pair(0.0, 0.0);
  } else {
    return std::make_pair((-b - sqrt(delta)) / 2, (-b + sqrt(delta)) / 2);
  }
}

// Exercise 8
double scalarProduct(const double *u, const double *v, size_t n) {
  double res(0.0);
  for (size_t ii = 0; ii <= n; ii++) {
    res += u[ii] * v[ii];
  }
  return res;
}

// Exercise 9
double myFunction(double x) { return x * x; }

double myOtherFunction(double x) { return x * x * x; }

// Exercise 10
double rectangular(double a, double b, double (*f)(double), size_t N) {
  double S(0.0);
  double step((b - a) / N);
  for (auto i = 0; i < N; i++) {
    S += f(a + i * step) * step;
  }
  return S;
}

double trapezoidal(double a, double b, double (*f)(double), size_t N) {
  double S(0.0);
  double step((b - a) / N);
  for (auto i = 0; i < N; i++) {
    S += step * (f(a + (i + 1) * step) + f(a + i * step)) / 2.0;
  }
  return S;
}

int main() {

  cout << "Exercice 1: " << endl;
  cout << power(2, 3) << endl;

  cout << "\nExercice 2: " << endl;
  cout << factorial(0) << endl;
  cout << factorial(1) << endl;
  cout << factorial(2) << endl;
  cout << factorial(10) << endl;
  cout << factorial(20) << endl;

  Matrix<int> test(2, 2, 0);
  test.print();
  cout << endl;
  test(1, 1) = 2;
  test(1, 0) = -1;
  test.print();
  cout << endl;

  cout << "\nExercice 3: " << endl;
  cout << "Pascal's triangle" << endl;
  print_pascal_triangle_fact(10);
  cout << endl;
  print_pascal_triangle_sum(10);
  cout << "\nExercice 4: ";
  cout << endl;
  cout << "Exp Taylor:" << endl;
  cout << endl << std::scientific << std::setprecision(14);
  cout << expTaylor(0);
  cout << endl;
  cout << expTaylor(1);
  cout << endl;
  cout << expTaylor(10);
  cout << endl;
  cout << expTaylor(30);
  cout << endl << endl;
  cout << "\nExercice 5: " << endl;
  cout << expPade(0);
  cout << endl;
  cout << expPade(1);
  cout << endl;
  cout << expPade(10);
  cout << endl;
  cout << expPade(30);
  cout << endl;

  cout << "\nComparing precision. Exp(14): " << std::exp(14) << endl;
  cout << " delta expTaylor10: " << std::abs(expTaylor(14)-std::exp(14)) << endl;
  cout << " delta expPade5 " << std::abs(expPade(14) - std::exp(14)) << endl;

  cout << "\nExercice 6: " << endl;
  double r{1.0};
  double theta{3.14 / 4};
  double pReal, pImaginary;

  CalculateRealAndImaginary(r, theta, pReal, pImaginary);
  cout << std::fixed << pReal << " Im: " << pImaginary;
  cout << endl;

  cout << "\nExercice 7: " << endl;
  const auto [x1, x2] = quadraticRoot(1, 4, 3);

  cout << x1 << " " << x2 << endl;

  cout << "\nExercice 8: " << endl;
  double u[] = {1, 2, 3, 4, 5};
  double v[] = {2, 2, 2, 2, 2};

  cout << scalarProduct(u, v, 5) << endl;

  cout << "\nExercice 9: " << endl;
  double (*p_function)(double x);

  p_function = myFunction;
  cout << p_function(3.0) << endl;
  p_function = myOtherFunction;
  cout << p_function(3.0) << endl;

  auto Id = [](double x) { return x; };

  cout << "\nExercice 10: " << endl;
  cout << "Rectangular" << endl;
  cout << rectangular(0, 1, Id, 10) << endl;
  cout << rectangular(0, 1, expPade, 10) << endl;
  cout << "Trapezoidal" << endl;
  cout << trapezoidal(0, 1, Id, 10) << endl;
  cout << trapezoidal(0, 1, expPade, 10) << endl;

  return 0;
}