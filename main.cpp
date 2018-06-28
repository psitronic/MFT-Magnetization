#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "fn.h"
#include "fn.C"

#define TINIT 1.0
#define TMAX 300.0
#define TSTEP 1.0
#define hINIT 0.01
#define hMAX 0.15
#define hSTEP 0.02

int main()
{
    const double J = 1.83; // The angular momentum
    const double Tc = 90.0; // The Curie temperature
    double T;
    double h; // The normalized field h=H/He, He - internal filed
    double m; // The reduced magnetization

    h = hINIT;

    do
    {
	T = TINIT;      

	do
	{
// Determination of the reduced magnetization m from the mean-field theory
// by solving of the self-consistency equation m(h,T)=B(x), where B(x)
// is the Brillouin function (see also function FCN_Mag_FM() in fn.C)

	    magnetization_FM(h, T, Tc, J, &m);

	    std::cout << T << "\t" << m << std::endl;
	    
	} while((T += TSTEP) <= TMAX);

    } while((h += hSTEP) <= hMAX);

    return 0;
}
