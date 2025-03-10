#include <iostream>
#include <cmath>
#include <type_traits>

size_t calculate_m(double x)
{
  size_t m = 0;
  while (abs(x) >= .5)
  {
    x /= 2;
    m++;
  }

  return m;
}

template<typename T>
T binomial(size_t n, size_t k)
{
  if (k>n-k) return binomial(n, n-k);
  T result = 1.0;
  for (size_t i=0; i<k; i++) {
    result *= static_cast<T>(n-i) / static_cast<T>(i+1);
  }
  return result;
}

template<typename T>
T factorial(size_t n) {
  T result = T(1);
  for (size_t i = 2; i<=n; i++) {
    result *= (T)i;
  }
  return result;
}

template<typename T, size_t N>
T Pade<N>(T x)
{
  T P = T(0);
  for (size_t n = N; n >= 0; --n)
  {
    P = binomial(N, n) /
            (factorial(n)) *
            (binomial(2 * N, n)) +
        x * P;
  }

  return P;
}

double expPade(double x)
{

  size_t m = calculate_m(x);
  x /= std::pow(2.0, m);

  return std::pow((Pade<5>(x) / Pade<5>(-x)), std::pow(2UL, m));
}

int main()
{
  cout << "Exp Taylor:" << endl;
  cout << endl
       << std::scientific << std::setprecision(14);
  cout << expTaylor(0);
  cout << endl;
  cout << expTaylor(1);
  cout << endl;
  cout << expTaylor(5);
  cout << endl;
  cout << expTaylor(30);
  cout << endl
       << endl;
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
  cout << " delta expTaylor10: " << std::abs(expTaylor(14) - std::exp(14)) << endl;
  cout << " delta expPade5 " << std::abs(expPade(14) - std::exp(14)) << endl;

  return 0;
}