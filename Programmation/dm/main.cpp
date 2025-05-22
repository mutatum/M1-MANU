#include "mesh_1d.hpp"
#include "field.hpp"
#include "output_writer.hpp"
#include "time_loop.hpp"

/*template <typename T>
T flux_Burgers(T const &u) // Useless
{
    return 0.5 * u * u;
}

template <typename T, typename Func>
T LF_numerical_flux(T const &ul, T const &ur, T max_vel, Func flux)
{
    return 0.5 * (flux(ul) + flux(ur)) - 0.5 / *(ur - ul);
}*/

/* 
   If function parameters are to be arbitrary, we need an instance of an object
 that has generic (virtual I guess) fluxes and numerical_fluxes, capable of
 calling numerical_fluxes with it's corresponding parameters.
 To rephrase, if received params of num_flux can change, then calling must account for it.
 It seems to me that rather than variadic parameters or bundled together params into an object,
 it is more readable and better suited overall (readability, correspondence to the equation formulation ..)
 to have a class that represents both problem and perhaps the numerical flux. 

 ##Later: Gemini 2.5 Pro suggested something akin to this)
*/

template <typename X>
class flux_problem
{
public:
    //    X flux()
protected:
}; // Starting to think about ex6)

/*
    Author: Choquert Ronan
    Completion: exercices 0-5 done on 0% internet (new appartment and no cellular data left), two C++ books (no complete reference book :( )
    Although 5 is partial since I couldn't resist thinking about ex6.
    I just now prompted Gemini that gave me a very good code structure proposal. I'll start implementing it,
    keeping as a rule 0 copy pasting. I'm conflicted as I very heavily dislike having code that doesn't
    originate from my mind, but I'm afraid I can't erase from my memory how seducing the proposed layout seemed.
    I'll post this version with it's mistakes and inefficiencies along side the cleaner one.
*/

int main()
{
    double x_left = 0.; // needed initialization parameters
    double x_right = 15.;
    size_t Ne = 1000;
    double final_time = 6.;


    mesh_1d<double> msh(x_left, x_right, Ne); // mesh construction
    //msh.print();
    field<double> field{GAUSS_PULSE, msh};                // initial field defined and initialized on mesh
    output_writer<double> out_stream(msh, "GAUSS_PULSE"); // initialize writer

    out_stream.write_solution(field, "initial"); // write initial field

    time_loop<double, double> loop_in_time(msh, final_time); // initialize time−loop
    loop_in_time.run(field);                                 // launch time−loop

    out_stream.write_solution(field,"final");// write final field
    
    return 0;
}