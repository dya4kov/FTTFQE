#pragma once

#include <gsl/gsl_matrix.h>
#include "SpecialFunctions.h"
#include <gsl/gsl_errno.h>
/**
 * @brief List of the ODE systems for Thomas-Fermi model
 */
namespace TF {
   /**
	* @brief ODE system for Thomas-Fermi potential.
	* @details In this system @f$ y_0(t) = \phi(t) @f$,
	* @f$ y_1(t) = \phi'(t) @f$, @f$ p @f$ - parameter:
	* @f{eqnarray*}{
	*  f_0 &=& y_1, \\
	*  f_1 &=& ptI_{1/2}
	*	  \left(
	*		\frac{y_0}{t}
	*	  \right).
	* @f}
	*/
    struct Potential {
	   /**
		* @brief Representation of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] f[] Current derivatives (left-hand side of the ODE).
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int system(double t, const double y[], double f[], void *params);
		/**
		* @brief Representation of jacobian of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
	   /**
		* @brief Is jacobian needed for this system?
		* @details For this system jacobian set as false. 
		*/
        static const bool isJacobian = false;
	   /**
		* @brief Dimension of the system.
		* @details Dimension is set to 2.
		*/
        static const int dimension = 2;
	   /**
	    * @brief Parameter @f$ p @f$ for the ODE.
		*/
        static double* parameter;
    };
	/**
	* @brief ODE system for Thomas-Fermi energy.
	* @details In this system 
	* @f$ y_0(t) = \phi(t) @f$,
	* @f$ y_1(t) = \phi'(t) @f$,
	* @f$ y_2(t) = E(t) @f$:
	* @f{eqnarray*}{
	*  f_0 &=& p_0 y_1, \\
	*  f_1 &=& p_1 t I_{1/2}
	*	  \left(
	*		\frac{y_0}{t}
	*	  \right), \\
	*  f_2 &=& p_2 t^2 I_{3/2}
	*				\left(
	*				  \frac{y_0}{t}
	*				\right).
	* @f}
	*/
    struct Energy {
		/**
		* @brief Representation of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] f[] Current derivatives (left-hand side of the ODE).
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int system(double t, const double y[], double f[], void *params);
	   /**
		* @brief Representation of jacobian of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
	   /**
		* @brief Is jacobian needed for this system?
		* @details For this system jacobian set as false.
		*/
        static const bool isJacobian = false; 
	   /**
		* @brief Dimension of the system.
		* @details Dimension is set to 3.
		*/
        static const int dimension = 3;
	   /**
		* @brief Array of parameters for the ODE.
		*/
        static double* parameter;
    };
	/**
	* @brief ODE system for Thomas-Fermi entropy.
	* @details In this system
	* @f$ y_0(t) = \phi(t) @f$,
	* @f$ y_1(t) = \phi'(t) @f$,
	* @f$ y_2(t) = S(t) @f$:
	* @f{eqnarray*}{
	*  f_0 &=& p_0 y_1, \\
	*  f_1 &=& p_1 t I_{1/2}
	*	  \left(
	*		\frac{y_0}{t}
	*	  \right), \\
	*  f_2 &=& p_2 t^2 I_{3/2}
	*				\left(
	*				  \frac{y_0}{t}
	*				\right)
	* @f}
	*/
    struct Entropy {
		/**
		* @brief Representation of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] f[] Current derivatives (left-hand side of the ODE).
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int system(double t, const double y[], double f[], void *params);
		/**
		* @brief Representation of jacobian of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		/**
		* @brief Is jacobian needed for this system?
		* @details For this system jacobian set as false.
		*/
        static const bool isJacobian = false;
		/**
		* @brief Dimension of the system.
		* @details Dimension is set to 3.
		*/
        static const int dimension = 3;
		/**
		* @brief Array of parameters for the ODE.
		*/
        static double* parameter;
    };
}

/**
* @brief List of the ODE systems for quantum and exchange corrections to the Thomas-Fermi model
*/
namespace TFCorrection {
	/**
	* @brief ODE system for quantum and exchange corrections to the Thomas-Fermi potential.
	* @details In this system 
	* @f$ y_0(t) = \phi(t) @f$,
	* @f$ y_1(t) = \phi'(t) @f$,
	* @f$ y_2(t) = \psi(t) @f$,
	* @f$ y_3(t) = \psi'(t) @f$,
	* @f$ p @f$ - parameter:
	* @f{eqnarray*}{
	*  f_0 &=& y_1, \\
	*  f_1 &=& ptI_{1/2}
	*	  \left(
	*		\frac{y_0}{t}
	*	  \right), \\
	*  f_2 &=& y_2, \\
	*  f_3 &=& p\left[
	*		       \frac12 y_2 I_{-1/2}\left(\frac{y_0}{t}\right) + tY'(t)
	*           \right].
	* @f}
	*/
    struct Potential {
		/**
		* @brief Representation of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] f[] Current derivatives (left-hand side of the ODE).
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int system(double t, const double y[], double f[], void *params);
		/**
		* @brief Representation of jacobian of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		/**
		* @brief Is jacobian needed for this system?
		* @details For this system jacobian set as false.
		*/
        static const bool isJacobian = false;
		/**
		* @brief Dimension of the system.
		* @details Dimension is set to 4.
		*/
        static const int dimension = 4;
		/**
		* @brief Parameter @f$ p @f$ for the ODE.
		*/
        static double* parameter;
    };
	/**
	* @brief ODE system for quantum and exchange corrections to the Thomas-Fermi energy.
	* @details In this system
	* @f$ y_0(t) = \phi(t) @f$,
	* @f$ y_1(t) = \phi'(t) @f$,
	* @f$ y_2(t) = \psi(t) @f$,
	* @f$ y_3(t) = \psi'(t) @f$,
	* @f$ y_4(t) = \Delta E(t) @f$,
	* @f$ p_0, p_1 @f$ - parameters:
	* @f{eqnarray*}{
	*  f_0 &=& y_1, \\
	*  f_1 &=& p_0 tI_{1/2}
	*	  \left(
	*		\frac{y_0}{t}
	*	  \right), \\
	*  f_2 &=& y_2, \\
	*  f_3 &=& p_0\left[
	*		       \frac12 y_2 I_{-1/2}\left(\frac{y_0}{t}\right) + tY'(t)
	*           \right], \\
	*  f_4 &=& -p_1 t \left[
	*				\frac12 y_2 I_{1/2}\left(\frac{y_0}{t}\right) + tY(t)
	*           \right].
	* @f}
	*/
    struct Energy {
		/**
		* @brief Representation of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] f[] Current derivatives (left-hand side of the ODE).
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int system(double t, const double y[], double f[], void *params);
		/**
		* @brief Representation of jacobian of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		/**
		* @brief Is jacobian needed for this system?
		* @details For this system jacobian set as false.
		*/
        static const bool isJacobian = false;
		/**
		* @brief Dimension of the system.
		* @details Dimension is set to 5.
		*/
        static const int dimension = 5;
		/**
		* @brief Array of parameters @f$ p_0 @f$ and @f$ p_1 @f$ for the ODE.
		*/
        static double* parameter;
    };
	/**
	* @brief ODE system for quantum and exchange corrections to the Thomas-Fermi entropy.
	* @details In this system
	* @f$ y_0(t) = \phi(t) @f$,
	* @f$ y_1(t) = \phi'(t) @f$,
	* @f$ y_2(t) = \psi(t) @f$,
	* @f$ y_3(t) = \psi'(t) @f$,
	* @f$ y_4(t) = \Delta S(t) @f$,
	* @f$ p_0, p_1 @f$ - parameters:
	* @f{eqnarray*}{
	*  f_0 &=& y_1, \\
	*  f_1 &=& p_0 tI_{1/2}
	*	  \left(
	*		\frac{y_0}{t}
	*	  \right), \\
	*  f_2 &=& y_2, \\
	*  f_3 &=& p_0\left[
	*		       \frac12 y_2 I_{-1/2}\left(\frac{y_0}{t}\right) + tY'(t)
	*           \right], \\
	*  f_4 &=& -p_1 t \left[
	*				\frac12 y_2 I_{1/2}\left(\frac{y_0}{t}\right) + 2tY(t)
	*           \right].
	* @f}
	*/
    struct Entropy {
		/**
		* @brief Representation of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] f[] Current derivatives (left-hand side of the ODE).
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int system(double t, const double y[], double f[], void *params);
		/**
		* @brief Representation of jacobian of the ODE system.
		* @param [in] t Current point.
		* @param [in] y[] Current solution.
		* @param [in] *params Void pointer to some parameters of the ODE.
		*/
        static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		/**
		* @brief Is jacobian needed for this system?
		* @details For this system jacobian set as false.
		*/
        static const bool isJacobian = false;
		/**
		* @brief Dimension of the system.
		* @details Dimension is set to 5.
		*/
        static const int dimension = 5;
		/**
		* @brief Array of parameters @f$ p_0 @f$ and @f$ p_1 @f$ for the ODE.
		*/
        static double* parameter;
    };
}

