#pragma once
#include <cstddef>
#include <vector>
#include <cmath>
#include <algorithm>
#include "mesh_1d.hpp"

enum function
{
    GAUSS_PULSE,
};

template <typename T>
class residual;

template <typename T>
class field
{
public:
    field(function label, mesh_1d<T> const &msh);

    inline const size_t n_cells() const { return this->n_cells_; };
    inline const std::vector<T> &operator()() const { return this->values_; };
    inline const T operator()(size_t ii) const { return this->values_.at(ii); };
    inline const std::string label() const { return this->label_; }

    field<T> &operator+=(residual<T> const &res);

protected:
    std::string label_{"IDENTITY"};
    std::vector<T> values_;
    size_t n_cells_{1};
};

template <typename T>
field<T>::field(function fname, mesh_1d<T> const &msh) : n_cells_(msh.n_cells()),
                                                         values_(std::vector<T>(msh.n_cells()))
{
    T (*func)(T) = [](T xx) -> T
    { return xx; };
    switch (fname)
    {
    case GAUSS_PULSE:
        this->label_ = "GAUSS_PULSE";
        func = [](T xx) -> T
        { return std::exp(-std::pow((xx - 5.), 2)); };
        break;
    }

    auto centers = msh.centers();
    std::transform(centers.begin(), centers.end(), this->values_.begin(), func);
}
template <typename T>
field<T> &field<T>::operator+=(residual<T> const &res)
{
    for (size_t ii = 0; ii < this->n_cells_; ii++)
    {
        this->values_.at(ii) += res(ii);
    }
    return *this;
}