#pragma once

#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include "ODESystem.h"
#include "Solution.h"

#define rk2 gsl_odeiv2_step_rk2
#define rkf45 gsl_odeiv2_step_rkf45
#define rk4 gsl_odeiv2_step_rk4
#define rkck gsl_odeiv2_step_rkck
#define rk8pd gsl_odeiv2_step_rk8pd
#define rk2imp gsl_odeiv2_step_rk2imp
#define rk4imp gsl_odeiv2_step_rk4imp
#define bsimp gsl_odeiv2_step_bsimp
#define rk1imp gsl_odeiv2_step_rk1imp
#define msadams gsl_odeiv2_step_msadams
#define msbdf gsl_odeiv2_step_msbdf

typedef gsl_odeiv2_step_type stepType;

/**
* @author Sergey Dyachkov
* @version 1.0
* @date 22.12.2013
* @bug currently unknown
* @brief A class solving an ODE to evaluate a specific quantity.
* @details An object-oriented superstructure of GSL ODE solver.
* A solver for a specific equation for some quantity. 
* All parameters about system are taken directly from struct
* with quantity description.
*/

template <class Quantity>
class ODESolver {
public:
   /**
	* @brief An ODESolver constructor.
	* @param [in] sType The type of GSL solver. Default is rkf45.
	* @param [in] hstart Start step size. Default is 1e-6.
	* @param [in] epsAbsolute The acceptable absolute error level. Default is 1e-6.
	* @param [in] epsRelative The acceptable relative error level. By default is unsetted.
	* @details Constructs a GSL solver with desired parameters of accuracy and step size
	* for a specific quantity. It makes both the driver with adaptive step control and 
	* single step evolver.
	*/
    ODESolver(const stepType* sType = rkf45,
                             double hstart = 1e-6,
                             double epsAbsolute = 1e-6, 
                             double epsRelative = 0.0)
    {
        system = new ODESystem<Quantity>;
        solution = new Solution(system->getDimension());
        epsAbs = epsAbsolute;
        epsRel = epsRelative;
        h = hstart;
        this->sType = sType;
        step = gsl_odeiv2_step_alloc(sType, system->getDimension());
        control = gsl_odeiv2_control_y_new(epsAbs, epsRel);
        evolve = gsl_odeiv2_evolve_alloc(system->getDimension());
        driver = gsl_odeiv2_driver_alloc_y_new(system->toGSL(), sType, h, epsAbs, epsRel);
        y0 = new double[system->getDimension()];
        for (int i = 0; i < system->getDimension(); ++i) {
            y0[i] = 0;
        }	
    }
   /**
	* @brief An ODESolver destructor.
	* @details Frees all the memory of GSL solver functions and system.
	*/
    ~ODESolver(void) {
        gsl_odeiv2_driver_free(driver);
        gsl_odeiv2_evolve_free(evolve);
        gsl_odeiv2_control_free(control);
        gsl_odeiv2_step_free(step);
        delete system;
        delete solution;
    }
   /**
	* @brief Sets absolute and relative accuracy of the solver.
	* @param [in] epsAbsolute Absolute acceptable error.
	* @param [in] epsRelative Relative acceptable error.
	*/
    void setAccuracy(double epsAbsolute, double epsRelative) {
        epsAbs = epsAbsolute;
        epsRel = epsRelative;
    }
   /**
	* @brief Sets step size.
	* @param [in] h Step size.
	* @details If left limit of integration is bigger than right than step is set to negative.
	*/
    void setStep(double h) {
        gsl_odeiv2_driver_free(driver);
        if (this->h > 0) {
            this->h = fabs(h);
            driver = gsl_odeiv2_driver_alloc_y_new(system->toGSL(), sType, this->h, epsAbs, epsRel);
        }
        else {
            this->h = -fabs(h);
            driver = gsl_odeiv2_driver_alloc_y_new(system->toGSL(), sType, this->h, epsAbs, epsRel);
        }
    }
   /**
    * @brief Sets limits of integration.
	* @param [in] xFrom Start point of integration.
	* @param [in] xTo End point of integration.
	* @details If start point is bigger than end point the step size is turned to negative. 
    */
    void setLimits(double xFrom, double xTo) {
        this->xFrom = xFrom;
        this->xTo = xTo;
        solution->x = xFrom;
        if (xFrom > xTo) {
            gsl_odeiv2_driver_free(driver);
            h = -fabs(h);
            driver = gsl_odeiv2_driver_alloc_y_new(system->toGSL(), sType, h, epsAbs, epsRel);
        }
    }
   /**
	* @brief Sets the initial values of the solution.
	* @param [in] Pointer to array of initials values of type double.
	*/
    void setInitials(double* y0) {
        for (int i = 0; i < system->getDimension(); ++i) {
            this->y0[i] = y0[i];
            solution->y[i] = y0[i];
        }
    }
   /**
    * @brief Apply single step with adaptive step size control.
	* @return The struct with current solution y and point x.
	*/
    Solution evolveApply() {
        status = gsl_odeiv2_evolve_apply(evolve, 
                                         control, 
                                         step, 
                                         system->toGSL(), 
                                         &solution->x, 
                                         xTo, &h, 
                                         solution->y);
        return *solution;
    }
   /**
	* @brief Apply single step with fixed step.
	* @param [in] fixedStep Value of fixed step.
	* @return The struct with current solution y and point x.
	*/
    Solution evolveApplyFixedStep(const double fixedStep) {
            status = gsl_odeiv2_evolve_apply_fixed_step(evolve, 
                                                        control, 
                                                        step, 
                                                        system->toGSL(), 
                                                        &solution->x, 
                                                        fixedStep, 
                                                        solution->y);
            return *solution;
    }
   /**
	* @brief Apply driver with fixed step.
	* @param [in] fixedStep Value of fixed step.
	* @return The struct solution y.
	*/
    Solution driverApplyFixedStep(const double fixedStep) {
        unsigned long n = static_cast<unsigned long>(floor((xTo - xFrom)/fixedStep));
        status = gsl_odeiv2_driver_apply_fixed_step(driver, 
                                                    &solution->x, 
                                                    fixedStep, n, 
                                                    solution->y);

        return *solution;
    }
   /**
	* @brief Apply driver with adaptive step size conrol.
	* @param [in] hMin Value of minimum step step. Default is provided by GSL.
	* @param [in] nMax Maximum number of steps. Default is unset.
	* @return The struct solution y.
	*/
    Solution driverApply(double hMin = 0, unsigned long nMax = 0) {
        gsl_odeiv2_driver_set_hmin(driver, hMin);
        gsl_odeiv2_driver_set_nmax(driver, nMax);
        status = gsl_odeiv2_driver_apply(driver, &solution->x, xTo, solution->y);
        return *solution;
    }
   /**
	* @brief Reset solver to its start conditions.
	* @param [in] h Value of new start step size. Default is 1e-6.
	*/
    void reset(double h = 1e-6) {
        for (int i = 0; i < system->getDimension(); ++i) {
                solution->y[i] = this->y0[i];
        }
        solution->x = xFrom;
        gsl_odeiv2_driver_reset(driver);
        this->h = h;
    }

private:

    ODESystem<Quantity>* system;

    double xFrom;
    double xTo;
    double h;
    double epsAbs;
    double epsRel;
    double* y0;
    int status;

    const stepType* sType;

    gsl_odeiv2_step* step;
    gsl_odeiv2_control* control;
    gsl_odeiv2_evolve* evolve;
    gsl_odeiv2_driver* driver;

    Solution* solution;

};