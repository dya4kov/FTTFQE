#include "../header/ODESystemList.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TF::Potential::system(double t, const double y[], double f[], void *params) {
    f[0] = y[1];
    if (t > 0) f[1] = *parameter*t*FermiDirac::Half::function(y[0]/t);
    else f[1] = 1e+20;
    return GSL_SUCCESS;
}
int TF::Potential::jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) { return 0; }
double* TF::Potential::parameter = new double;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TF::Energy::system(double t, const double y[], double f[], void *params) {
    f[0] = parameter[0]*y[1];
    if (t > 0) {
        f[1] = parameter[1]*t*FermiDirac::Half::function(y[0]/t);
        f[2] = parameter[2]*FermiDirac::ThreeHalf::function(y[0]/t)*t*t;
    }
    else { 
        f[1] = f[2] = 1e+20;
    }
    return GSL_SUCCESS;
}
int TF::Energy::jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) { return 0; }
double* TF::Energy::parameter = new double[dimension];
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TF::Entropy::system(double t, const double y[], double f[], void *params) {
    f[0] = parameter[0]*y[1];
    if (t > 0) {
        f[1] = parameter[1]*t*FermiDirac::Half::function(y[0]/t);
        f[2] = parameter[2]*FermiDirac::ThreeHalf::function(y[0]/t)*t*t;
    }
    else { 
        f[1] = f[2] = 1e+20;
    }
   return GSL_SUCCESS;
}
int TF::Entropy::jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) { return 0; }
double* TF::Entropy::parameter = new double[dimension];
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TFCorrection::Potential::system(double t, const double y[], double f[], void *params) {
    f[0] = y[1];
    f[2] = y[3];
    if (t > 0) {
        f[1] = *parameter*t*FermiDirac::Half::function(y[0]/t);
        f[3] = *parameter*(FermiDirac::MHalf::function(y[0]/t)/2.0*y[2] + t*Y::derivative(y[0]/t));
    }
    else  {
        f[1] = f[3] = 1e+20;
    }
    return GSL_SUCCESS;
}
int TFCorrection::Potential::jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) { return 0; }
double* TFCorrection::Potential::parameter = new double;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TFCorrection::Energy::system(double t, const double y[], double f[], void *params)
{
    static int sign;
    f[0] = y[1];
    f[2] = y[3];
    if (t > 0) {
        f[1] = parameter[0]*t*FermiDirac::Half::function(y[0]/t);
        f[3] = parameter[0]*(FermiDirac::MHalf::function(y[0]/t)/2.0*y[2] + t*Y::derivative(y[0]/t));
        f[4] = -parameter[1]*t*(0.5*y[2]*FermiDirac::Half::function(y[0]/t) + t*Y::function(y[0]/t));
        f[4] > 0 ? sign = 1 : sign = -1;
    }
    else {
        f[1] = f[3] = 1e+20;
        f[4] = sign*1e+20;
    }
    return GSL_SUCCESS;
}

int TFCorrection::Energy::jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) { return 0; }
double* TFCorrection::Energy::parameter = new double[2];
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TFCorrection::Entropy::system(double t, const double y[], double f[], void *params)
{
    static int sign;
    f[0] = y[1];
    f[2] = y[3];
    if (t > 0) {
        f[1] = parameter[0]*t*FermiDirac::Half::function(y[0]/t);
        f[3] = parameter[0]*(FermiDirac::MHalf::function(y[0]/t)/2.0*y[2] + t*Y::derivative(y[0]/t));
        f[4] = -parameter[1]*t*(0.5*y[2]*FermiDirac::Half::function(y[0]/t) + 2*t*Y::function(y[0]/t));
        f[4] > 0 ? sign = 1 : sign = -1;
    }
    else {
        f[1] = f[3] = 1e+20;
        f[4] = sign*1e+20;
    }
    return GSL_SUCCESS;
}

int TFCorrection::Entropy::jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) { return 0; }
double* TFCorrection::Entropy::parameter = new double[2];
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
