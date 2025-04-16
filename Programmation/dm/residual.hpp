#pragma once
#include "mesh_1d.hpp"
#include "field.hpp"

template <typename T>
class residual
{
public:
    residual(mesh_1d<T> const &mesh) : mesh_(mesh), values_(mesh.n_cells()) {}
    void assemble_from_field(field<T> const &uh);

    const mesh_1d<T> &mesh() const { return mesh_; }
    inline const size_t n_cells() const { return mesh_.n_cells(); }

    T operator()(size_t ii) const { return values_.at(ii); }

    template <typename S>
    friend residual<S> operator*(S lambda, residual<S> const &res);

protected:
    std::vector<T> values_;
    mesh_1d<T> const &mesh_;
};

// Burgers and LF numerical flux specific code
template <typename T>
T flux_Burgers(T const &u)
{
    return 0.5 * u * u;
}

template <typename T, typename Func>
T LF_numerical_flux(T const &ul, T const &ur, T max_vel, Func flux)
{
    return 0.5 * (flux(ul) + flux(ur)) - 0.5 * max_vel * (ur - ul);
}

template <typename T>
void residual<T>::assemble_from_field(field<T> const &uh)
{
    assert(uh.n_cells() > 1);
    T max_vel{0.};
    for (auto value : uh())
        max_vel = max_vel < std::abs(value) ? std::abs(value) : max_vel;
    assert(max_vel >= 0);
    for (size_t ii = 1; ii < uh().size() - 1; ii++)
    {
        this->values_.at(ii) = -(LF_numerical_flux(uh(ii), uh(ii + 1), max_vel, flux_Burgers<T>) - LF_numerical_flux(uh(ii - 1), uh(ii), max_vel, flux_Burgers<T>))/this->mesh_.dx();
    }
    this->values_.at(0) = this->values_.at(1);
    size_t last = uh.n_cells() - 1;
    this->values_.at(last) = this->values_.at(last - 1);
}

template <typename S>
residual<S> operator*(S lambda, residual<S> const &res)
{
    residual<S> mres(res);
    for (auto &val : mres.values_)
    {
        val *= lambda;
    }
    return mres;
}