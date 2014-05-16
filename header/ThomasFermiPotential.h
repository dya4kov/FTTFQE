#pragma once
#include "ODESolver.h"
#include "ODESystemList.h"
#include <fstream>
#include "Solution.h"
/**
 * @brief This class implements interface for calculating Thomas-Fermi potential.
 * @details Thomas-Fermi potential is calculated from the following boundary problem:
 * @f{eqnarray*}{
 *  \frac{d^2\phi}{dx^2} &=& a x I_{1/2}
 *	  \left(
 *		\frac{\phi(x)}{x}
 *	  \right), \\
 *  \phi(0) &=& \frac{1}{T}\sqrt{\frac{4\pi}{3V}}, \\
 *  \phi(1) &=& \phi'(1).
 * @f}
 * Here @f$ a = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$.
 */
class TFPotential {
public:
   /**
	* @brief Constructor of Thomas-Fermi potential.
	* @details Allocates memory for potential solver
	* and loads table of precalculated values of potential at
	* @f$ x = 1 @f$.
	*/
    TFPotential(void);
   /**
	* @brief Destructor of Thomas-Fermi potential.
	*/
    ~TFPotential(void);
   /**
	* @brief Value of potential @f$ \phi @f$ at @f$ x = 1 @f$.
	* @details Solving boundary problem by shooting method with 
	* good start values of @f$ \phi(1) @f$ from table for fast 
	* convergency.
	* @return @f$ \phi(1) @f$
	*/
    double valueAt_1();
	/**
	* @brief Value of potential @f$ \phi @f$ at some point @f$ x @f$.
	* @return @f$ \phi(x) @f$
	*/
    double valueAt_x(double x);
	/**
	* @brief Value of derivative of potential @f$ \phi'(0) @f$.
	* @return @f$ \phi'(0) @f$
	*/
    double derivativeAt_0();
	/**
	* @brief Set volume and temperature values for calculating potential.
	*/
    void setParameters(double V, double T);
	/**
	* @brief Set precision for calculating potential.
	*/
    void setPrecision(double eps);

private:
    void setPhiTable();
    void calculate();
    void setInitialShotParameters(double& p1, double& p2);

    const int vTableSize;
    const int tTableSize;
    const double lgV0;
    const double lgT0;

    const double lgVStep;
    const double lgTStep;

    double phi_0;
    double phi_1;
    double V;
    double T;
    double eps;

    double** phiTable;

    ODESolver<TF::Potential>* solver;

    bool calculated;

};

