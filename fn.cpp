// The Brilloin function
inline double brillouin(double x, double J)
{
    double T1, T2;

    T1 = (tanh(((2.0 * J + 1.0) * x) / (2.0 * J)));
    T2 = (tanh(x / (2.0 * J)));

    return ((2.0 * J + 1.0)/(2.0 * J * T1)) - (1.0 / (2.0 * J * T2));
}

inline double FCN_Mag_FM(double m, void *params)
{
    struct params *p = (struct params *) params;

    double h = p->h;
    double Temp = p->Temp;
    double Tc = p->Tc;
    double J = p->J;
    double x;

//In the FM case, according to the mean-field theory, x=((m+h)/t * 3J/(J+1)), t=T/Tc

    x = ((3.0 * J)/(J + 1.0)) * ((Tc * (m + h)) / Temp);

    return (brillouin(x, J) - m);
}

// Solving of the self-consistency equation m(h,T)=B(x), where B(x)
// is the Brillouin function (see also function FCN_Mag_FM())
void magnetization_FM(double h, double Temp, double Tc, double J, double *M)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;
    double x_lo = 0.000001, x_hi = 1.0;
    gsl_function F;
    struct params params = {h, Temp, Tc, J};

    F.function = &FCN_Mag_FM;
    F.params = &params;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    do
    {
	iter++;
	status = gsl_root_fsolver_iterate (s);
	r = gsl_root_fsolver_root (s);
	x_lo = gsl_root_fsolver_x_lower (s);
	x_hi = gsl_root_fsolver_x_upper (s);
	status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

	if (status == GSL_SUCCESS)
	    *M = r;
    }
    while(status == GSL_CONTINUE && iter < max_iter);
}

/*void magnetization(double h, double T, double Tc, double J, double *M)
{
    double x = 0.0;
    double F_left, F_right;
   
    F_left = FCN_Mag(h, T, Tc, J, x);

    do 
    {
	F_right = FCN_Mag(h, T, Tc, J, x);

	if ((F_left * F_right) < 0)
	{
	    *M = x;
	    break;
	}
	else
	{
	    F_left = F_right;
	    *M = 0.0;
	}

    } while((x += 0.01) <= 1.1);

}
*/
