#pragma once
#include <cstddef>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>

template <typename T>
class mesh_1d
{

protected:
    std::vector<T> vertices_coordinate_;
    std::vector<T> cells_center_;
    std::size_t n_vertices_{0};
    std::size_t n_cells_{1};
    T dx_;

public:
    mesh_1d(T x_left, T x_right, size_t n_cells);
    ~mesh_1d();

    inline const size_t n_vertices() const { return this->n_vertices_; }
    inline const size_t n_cells() const { return this->n_cells_; }
    inline const T dx() const { return this->dx_; }

    inline const T xc(size_t ii) const { return this->cells_center_.at(ii); }
    inline const std::vector<T> &centers() const { return this->cells_center_; }; // Not necessarily const..
    void print();
};

template <typename T>
mesh_1d<T>::mesh_1d(T x_left, T x_right, size_t n_cells) : cells_center_(std::vector<T>(n_cells, T(0))),
                                                           vertices_coordinate_(std::vector<T>(n_cells + 1, T(0))),
                                                           n_vertices_(n_cells + 1),
                                                           n_cells_(n_cells),
                                                           dx_((x_right-x_left)/n_cells)
{
    assert(x_left < x_right);
    assert(n_cells > 0);

    for (std::size_t ii = 0; ii < n_cells; ii++)
    {
        this->vertices_coordinate_.at(ii) = x_left + ii * this->dx_;
        this->cells_center_.at(ii) = x_left + (ii+.5) * this->dx_;
    }
    this->vertices_coordinate_.at(n_cells) = x_right;
}

template <typename T>
mesh_1d<T>::~mesh_1d()
{
}

template <typename T>
void mesh_1d<T>::print()
{
    std::cout.precision(std::ceil(std::max(4, int(-std::log10(this->dx_)))));
    std::cout << "1D Mesh representing ["
              << *this->vertices_coordinate_.begin()
              << ", " << this->vertices_coordinate_.back() << "] with "
              << n_cells() << " cells and dx of: " << dx() << "\n";
    std::cout << "\nVertices: \n{";
    for (const auto& v : this->vertices_coordinate_)
    {
        std::cout << v;
        if (&v != &this->vertices_coordinate_.back())
            std::cout << ", ";
    }
    std::cout << "}\n\nCells center: \n{";
    for (const auto& c : this->cells_center_)
    {
        std::cout << c;
        if (&c != &(this->cells_center_.back()))
            std::cout << ", ";
    }
    std::cout << "}\n";
}