#pragma once
#include <string>
#include <fstream>
#include "mesh_1d.hpp"
#include "field.hpp"

template <typename T>
class output_writer
{
public:
    output_writer(mesh_1d<T> const &msh, std::string filename) : mesh_(msh), radical_name_(filename) {};
    void write_solution(field<T> solution, std::string loop_counter);
    inline const mesh_1d<T> &mesh() const { return mesh_; };

protected:
    mesh_1d<T> const &mesh_;
    std::string radical_name_;
};

template <typename T>
void output_writer<T>::write_solution(field<T> solution, std::string loop_counter)
{
    std::ofstream ofs{solution.label() + '_' + loop_counter + ".txt"}; // implicit ofstream::trunc
    if (!ofs)
        std::cerr << "couldn't open " << solution.label() << " for writing\n";
    ofs << loop_counter << '\n';
    ofs << "X\tY" << '\n';
    ofs.precision(12);
    for (size_t ii{0}; ii < mesh_.n_cells(); ii++)
    {
        ofs << mesh_.xc(ii) << ' ' << solution(ii) << '\n';
    }
    ofs << std::endl;
    ofs.close();
}