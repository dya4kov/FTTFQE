#pragma once

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <math.h>
/**
* @brief List of the Fermi-Dirac functions.
* @details
* Fermi-Dirac functions have the following representation:
* @f[
*	I_k(x) = \int_0^\infty\frac{y^kdy}{1 + exp(y - x)}.
* @f]
* Derivative of the Fermi-Dirac function at k > 0:
* @f[
*	I_k'(x) = kI_{k - 1}(x).
* @f]
*/
namespace FermiDirac {
	/**
	* @brief Derivative of the function @f$ I_{-1/2}(x) @f$.
	* @details We need to use good approximation of the function
	* @f$ I_{-1/2}(x) @f$ to calculate its derivative, 
	* because in this case @f$ k < 0 @f$.
	* @return Value of the function at the point x.
	*/
    struct DMHalf {
        static double function(double x);
    };
	/**
	* @brief Function @f$ I_{-1/2}(x) @f$.
	* @return Value of the function at the point x.
	*/
    struct MHalf {
        static double function(double x);
    };
	/**
	* @brief Function @f$ I_{1/2}(x) @f$.
	* @return Value of the function at the point x.
	*/
    struct Half {
        static double function(double x);
    };
	/**
	* @brief Function @f$ I_{3/2}(x) @f$.
	* @return Value of the function at the point x.
	*/
    struct ThreeHalf {
        static double function(double x);
    };
	/**
	* @brief Function @f$ I_{-1/2}^2(x) @f$.
	* @return Value of the function at the point x.
	*/
    struct MHalfSquared {
        static double function(double x);
    };
}

/**
* @brief Representation of the function @f$ Y(x) @f$.
* @details Function @f$ Y(x) = \frac12 I_{1/2}(x) I_{-1/2}(x) 
* + \frac32 \int_{-\infty}^x I_{-1/2}^2(t)dt @f$
* is integrated numerically using GSL qags integration.
*/
struct Y {
	/**
	* @return Value of the function @f$ Y(x) @f$.
	*/
    static double function(double x);
	/**
	* @return Value of the derivative @f$ Y'(x) @f$.
	*/
    static double derivative(double x);
	/**
	* @brief Integration function for GSL integrator.
	*/
    static double integrationFunction(double x, void* params);
	/**
	* @brief Workspace for GSL integrator.
	*/
    static gsl_integration_workspace* w;
	/**
	* @brief GSL function for integrator.
	*/
    static gsl_function F;
};