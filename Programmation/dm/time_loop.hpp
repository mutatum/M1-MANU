#pragma once
#include "mesh_1d.hpp"
#include "field.hpp"
#include "residual.hpp"

template <typename T, typename X> // Should be some sort of `template <typename T, typename X>`, differenciating time from space
class time_loop
{
protected:
    size_t iteration_counter_{0};
    T dt_{0.};
    T physical_time_{0.};
    T final_time_{0.};
    // bool last_iteration_{false};
    double cfl_number_ = 0.5;

    residual<X> residual_;   // Should residual be able to change outside time_loop?
    mesh_1d<X> const &mesh_; // Mesh can change outside

public:
    time_loop(mesh_1d<X> const &msh, T final_time) : mesh_(msh), // I somewhat dislike this representation in light of ex6
                                                     final_time_(final_time),
                                                     residual_(msh) {};
    ~time_loop() {};

    inline const T final_time() const { return final_time_; };
    void run(field<X> &uh);
    inline void compute_dt(field<X> &uh);
};

template <typename T, typename X>
inline void time_loop<T, X>::compute_dt(field<X> &uh)
{
    T m{0.};
    for (auto value : uh())
        m = m < std::abs(value) ? std::abs(value) : m;
    assert(m >= 0.);

    this->dt_ = std::min(this->final_time_ - this->physical_time_, // can be 0, shouldn't be negative
                         this->cfl_number_ * (this->mesh_.dx() / m));
    assert(this->dt_ >= 0);
    assert(this->dt_ < final_time_);
}

template <typename T, typename X>
void time_loop<T, X>::run(field<X> &uh)
{
    this->compute_dt(uh);
    while (this->final_time_ - this->physical_time_ > 1e-12)
    {
        this->residual_.assemble_from_field(uh);
        uh += this->dt_ * this->residual_;
        this->physical_time_ += this->dt_;
        this->compute_dt(uh);
    }
}