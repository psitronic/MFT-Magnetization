struct params
{
    double h, Temp, Tc, J;
};
double FCN_Mag_FM(double m, void *params);
double brillouin(double x, double J);
void magnetization_FM(double h, double T, double Tc, double J, double *m);
