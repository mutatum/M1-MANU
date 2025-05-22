
#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using Eigen::Vector3d;

/* U = { rho,
         rho * v,
         rho * E} */

/* F(U) = { rho * v,
            rho * v tensor v + p * I,
            rho * E * v + p * v } */

template <size_t n>
double legendre(double x)
{
    switch (n)
    {
    case 0:
        return 1;
    case 1:
        return x;
    default:
        return (1 / (n + 1)) * ((2 * n + 1) * legendre<n - 1>(x) - n * legendre<n - 2>(x));
    }
}

// Quarteroni Sacco Saleri Méthodes Numériques p.348
// Noeuds et poids des formules de Gauss-Legendre
// returns [x,w] où w est un vector ligne de poids
template <std::size_t n>
Eigen::Matrix<double, 2, n> zplege()
{
    Eigen::Matrix<double, 2, n> M;
    Eigen::Matrix<double, n, n> JacM;
    M.setZero();
    JacM.setZero();

    // if (n<=1) return null; // error
    double b;
    for (std::size_t i = 0; i < n - 1; i++)
    {
        b = 1. / (4 - 1. / std::pow(i + 1, 2));
        JacM(i, i + 1) = JacM(i + 1, i) = std::sqrt(b); // sur et sous diag
    }
    std::cout << JacM << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, n, n>> eigensolver(JacM);
    if (eigensolver.info() != Eigen::Success)
        abort();
    auto values = eigensolver.eigenvalues();
    auto vectors = eigensolver.eigenvectors();
    std::cout << "\nThe eigenvalues of JacM are:\n"
              << values << std::endl;
    std::cout << "\nThe eigenvectors of JacM are:\n"
              << vectors << std::endl;
    for (std::size_t i = 0; i < n; i++)
    {
        M(0, i) = values(i);
        M(1, i) = vectors(0, i) * vectors(0, i) * 2;
    }
    return M;
}

template <size_t n>
Eigen::Matrix<double, n + 1, n + 1> build_mass_matrix()
{
    Eigen::Matrix<double, n + 1, n + 1> M;
    M.setZero();
    for (std::size_t i = 0; i < n + 1; i++)
    {
        M(i, i) = 2. / (2. * i + 1.); // <Lk, Lm> = delta{k,m} / (k+1/2)
    }
    return M;
}

inline double eq_of_state(double gamma, double rho, double e)
{
    return (gamma - 1) * rho * e;
}

Vector3d physical(Vector3d U)
{
    Vector3d F;
    double u = U[1] / U[0];
    double E = U[2] / U[0];
    double e = E - .5 * (u * u);
    double p = eq_of_state(1.4, U[0], e);
    F[0] = U[1]; // rho * v
    F[1] = U[1] * u + p;
    F[2] = (E + p) * u;
    return F;
}

/*  A note about expression templates

This is an advanced topic that we explain on this page, but it is useful to just mention it now.
In Eigen, arithmetic operators such as operator+ don't perform any computation by themselves,
they just return an "expression object" describing the computation to be performed.
The actual computation happens later, when the whole expression is evaluated, typically in operator=.
While this might sound heavy, any modern optimizing compiler is able to optimize away that
abstraction and the result is perfectly optimized code.
*/
int main()
{
    constexpr int ndof = 3;

    Eigen::Array<double, 1, 10> cells;  // 10 cellules
    auto M = build_mass_matrix<ndof>(); // M is n+1 square
    std::cout << "Mass Matrix:\n"
              << M << std::endl;
    std::cout << M.inverse() << std::endl;
    // std::cout << M*M.inverse() << std::endl; // this is indeed the inverse
    return 0;
}