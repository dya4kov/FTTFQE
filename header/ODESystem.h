#pragma once
#include <gsl/gsl_odeiv2.h>

/**
 * @author Sergey Dyachkov
 * @version 1.0
 * @date 22.12.2013
 * @bug currently unknown
 * @brief A class with parameters for ODE system.
 */

template <class Quantity>
class ODESystem {	
public: 
   /**
	* @brief An ODE system constructor.
	* @details Takes a quantity type as a template parameter
	* and makes a GSL ODE system from list with corresponding
	* quantity.
	*/
    ODESystem() {
        system = new gsl_odeiv2_system;
        system->dimension = dimension = Quantity::dimension;
        system->function = &Quantity::system;
        if (Quantity::isJacobian) system->jacobian = &Quantity::jacobian;
        else system->jacobian = NULL;
        system->params = NULL;
        *Quantity::parameter = 0;
    }
   /**
	* @brief An ODE system destructor.
	* @details Delete a GSL ODE system.
	*/
    ~ODESystem() {
        delete system;
    }
   /**
	* @brief A member returning a pointer to GSL system.
	*/
    gsl_odeiv2_system* toGSL() {
        return system;
    }
   /**
	* @brief A member returning a dimension of the system.
	*/
    int getDimension() {
        return dimension;
    }
private:
    gsl_odeiv2_system* system;
    int dimension;
};