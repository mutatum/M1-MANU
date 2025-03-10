# pragma once
# include <cstddef>
# include <iostream>


template<size_t M>
class Binomial{
    private:
        double coeff[M][M]={};

    public:
        Binomial();
        double operator()(size_t n, size_t k) const;
        void print() const;
};

// ==============================

template<size_t M>
Binomial<M>::Binomial() {
  for (size_t i = 1; i <= M; i++) {
    for (size_t j = 0; j < i; j++) {
      if (j == 0 || j == i - 1) {
        this->coeff[M - i + j][j] = 1;
      } else {
        this->coeff[M - i + j][j] =
            this->coeff[M - i + j][j - 1] + this->coeff[M - i + j + 1][j];
      }
    }
  }
}

template<size_t M>
double Binomial<M>::operator()(size_t n, size_t k) const {
    return this->coeff[M-n+k-1][k];
}

template<size_t M>
void Binomial<M>::print() const {
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < M; j++) {
            std::cout << this->coeff[i][j] << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}