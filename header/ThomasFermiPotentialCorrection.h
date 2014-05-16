#pragma once
#include "ODESolver.h"
#include "ODESystemList.h"
#include "ThomasFermiPotential.h"
#include <fstream>
#include "Solution.h"
/**
* @brief This class implements interface for calculating correction to the Thomas-Fermi potential.
* @details correction to the Thomas-Fermi potential is calculated from the following boundary problem:
* @f{eqnarray*}{
*  \frac{d^2\phi}{dx^2} &=& a x I_{1/2}
*	  \left(
*		\frac{\phi(x)}{x}
*	  \right), \\
*  \frac{d^2\psi}{dx^2} &=& a \left[
*		       \frac12 \psi(x) I_{-1/2}\left(\frac{\phi(x)}{x}\right) + xY'(x)
*           \right], \\
*  \phi(1) &=& \phi'(1) = \phi_1, \\
*  \psi(1) &=& \psi'(1), \\
*  \psi(0) &=& 0.
* @f}
* Here @f$ a = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$, @f$ \phi @f$ - Thomas-Fermi potential, 
* @f$ \psi @f$ - correction to Thomas-Fermi potential.
*/
class TFPotentialCorrection
{
public:
	/**
	* @brief Constructor of correction to the Thomas-Fermi potential.
	* @details Allocates memory for potential correction solver
	* and loads table of precalculated values of potential correction at
	* @f$ x = 1 @f$.
	*/
    TFPotentialCorrection(void);
	/**
	* @brief Destructor of correction to the Thomas-Fermi potential.
	*/
    ~TFPotentialCorrection(void);
	/**
	* @brief Value of correction to potential @f$ \psi @f$ at @f$ x = 1 @f$.
	* @details Solving boundary problem by shooting method with
	* good start values of @f$ \psi(1) @f$ from table for fast
	* convergency.
	* @return @f$ \psi(1) @f$
	*/
    double valueAt_1();
	/**
	* @brief Value of correction to potential @f$ \psi @f$ at some point @f$ x @f$.
	* @return @f$ \psi(x) @f$
	*/
    double valueAt_x(double x);
	/**
	* @brief Value of derivative of correction to potential @f$ \psi'(0) @f$.
	* @return @f$ \psi'(0) @f$
	*/
    double derivativeAt_0();
	/**
	* @brief Set volume and temperature values for calculating correction to potential.
	*/
    void setParameters(double V, double T);
	/**
	* @brief Set precision for calculating correction to potential.
	*/
    void setPrecision(double eps);

private:

    void setPsiTable();
    void calculate();
    void setInitialShotParameters(double& p1, double& p2);

    const int vTableSize;
    const int tTableSize;
    const double lgV0;
    const double lgT0;
    const double lgVStep;
    const double lgTStep;

    double psi_1;
    double V;
    double T;
    double eps;

    double** psiTable;

    ODESolver<TFCorrection::Potential>* solver;
    TFPotential* potential;

    bool calculated;
};

