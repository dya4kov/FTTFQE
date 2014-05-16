
#include "../header/SpecialFunctions.h"

double FermiDirac::DMHalf::function(double x) {
    static double p1[5] = {-1.253314128820,
                  -1.723663557701,
                  -6.559045729258E-01,
                  -6.342283197616E-02,
                  -1.488383106116E-05};
   static double q1[5] = {1.0,
                  2.191780925980,
                  1.605812955406,
                  4.443669527481E-01,
                  3.624232288112E-02};
   static double p2[5] = {1.0738127694,
                   5.6003303660,
                   3.6882211270,
                   1.1743392816,
                   2.3641935527E-01};
   static double q2[5] = {1.0,
                   4.6031840667,
                   4.3075910674E-01,
                   4.2151132145E-01,
                   1.1832601601E-02};
   static double p3[5] = {-8.222559330E-01,
                   -3.620369345E+01,
                   -3.015385410E+03,
                   -7.049871579E+04,
                   -5.698145924E+04}; 
   static double q3[5] = {1.0,
                    3.935689841E+01,
                    3.568756266E+03,
                    4.181893625E+04,
                    3.385138907E+05};
   static double* p;
   static double* q;
   double upSum = 0, downSum = 0;
   double DupSum = 0, DdownSum = 0;
   int i;
   if (x < 1.0)
   {
      p = p1; q = q1;
      for (i = 0; i < 5; i++)
      {
         upSum += p[i]*exp(i*x);
         downSum += q[i]*exp(i*x);
         DupSum += p[i]*i*exp(i*x);
         DdownSum += q[i]*i*exp(i*x);

      }
      return exp(x)*gsl_sf_gamma(1.0/2.0) + 2*exp(2*x)*upSum/downSum + exp(2*x)*(DupSum*downSum - upSum*DdownSum)/downSum/downSum;
   }
   if (x < 4.0)
   {
      p = p2; q = q2;
      for (i = 0; i < 5; i++)
      {
         upSum += p[i]*pow(x, i);
         downSum += q[i]*pow(x, i);
         DupSum += p[i]*i*pow(x, i - 1);
         DdownSum += q[i]*i*pow(x, i - 1);
      }
      return (DupSum*downSum - upSum*DdownSum)/downSum/downSum;
   }
   else
   {
      p = p3; q = q3;
      for (i = 0; i < 5; i++)
      {
         upSum += p[i]*pow(x, -2*i);
         downSum += q[i]*pow(x, -2*i);
         DupSum += -p[i]*2*i*pow(x, -2*i - 1);
         DdownSum += -q[i]*2*i*pow(x, -2*i - 1);

      }
      return 1/sqrt(x) - 3.0/2.0*pow(x, -5.0/2.0)*upSum/downSum + pow(x, -3.0/2.0)*(DupSum*downSum - upSum*DdownSum)/downSum/downSum;
   }
}

double FermiDirac::MHalf::function(double x) {
    if (x < -100.0) return 0;
    else return gsl_sf_gamma(1.0/2.0)*gsl_sf_fermi_dirac_mhalf(x);
}

double FermiDirac::Half::function(double x) {
    if (x < -100.0) return 0;
    else return gsl_sf_gamma(3.0/2.0)*gsl_sf_fermi_dirac_half(x);
}

double FermiDirac::ThreeHalf::function(double x) {
    if (x < -100.0) return 0;
    else return gsl_sf_gamma(5.0/2.0)*gsl_sf_fermi_dirac_3half(x);
}

double FermiDirac::MHalfSquared::function(double x) {
    double f = FermiDirac::MHalf::function(x); 
    return f*f;
}


gsl_function Y::F = {integrationFunction, NULL};
gsl_integration_workspace* Y::w = gsl_integration_workspace_alloc(5000);

double Y::integrationFunction(double x, void* params) {
    return FermiDirac::MHalfSquared::function(x);
}

double Y::derivative(double x) {
    return 7.0/4.0*FermiDirac::MHalfSquared::function(x) + 
            FermiDirac::Half::function(x)*FermiDirac::DMHalf::function(x)/2.0;
}

double Y::function(double x) {
    static double integral = 0;
    static double abserr = 0;

    gsl_integration_qags(&F, -200, x, 0, 1E-07, 5000, w, &integral, &abserr);

    return FermiDirac::Half::function(x)*FermiDirac::MHalf::function(x)/2.0 + 3.0/2.0*integral;
}
