#pragma once
#include "Model.h"
#include "ThomasFermiPotential.h"
#include "ThomasFermiPotentialCorrection.h"

/**
* @brief This class implements interface for
* quantum and exchange corrections to the Thomas-Fermi model.
* @details It implements virtual methods for calculating thermodynamic values.
*/
class ThomasFermiCorrection : public Model {
public:
	/**
	* @brief A constructor of corrections to the Thomas-Fermi model.
	* @details It allocates memory for ODE solvers of 
	* corrections to the Thomas-Fermi energy and entropy
	* and for Thomas-Fermi potential and correction to potential.
	*/
    ThomasFermiCorrection(Charge Z = 1.0);
	/**
	* @brief A destructor of corrections to the Thomas-Fermi model.
	* @details It frees all the allocated memory in the constructor.
	*/
    ~ThomasFermiCorrection(void);
	/**
	* @brief Get name for corrections to the Thomas-Fermi model as "TFC".
	* @return "TFC".
	*/
    std::string getName();

protected:
	/**
	* @brief Implemented method for calculating correction to pressure.
	* @details The formula for calculation:
	* @f[
	*	\Delta P = \frac{3T^2}{\pi^3}I_{1/2}(\phi(1))\psi(1) + Y(\phi(1)).
	* @f]
	*/
    void calculatePressure();
	/**
	* @brief Implemented method for calculating correction to energy.
	* @details The calculation goes through solving the ODE with parameters
	* @f$ p_0 = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ p_1 = \pi^{-3}VT^{2} @f$,
	* @f$ \Delta E_0 = -0.269900170 @f$. Correction to Thomas-Fermi energy 
	* is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& p_0 x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  \frac{d^2\psi}{dx^2} &=& p_0 \left[
	*		       \frac12 \psi(x) I_{-1/2}\left(\frac{\phi(x)}{x}\right) + xY'(x)
	*           \right], \\
	*  \Delta E'(x) &=& -p_1 x \left[
	*						\frac12 \psi(x) I_{1/2}\left(\frac{\phi(x)}{x}\right) + xY(x)
	*					\right], \\
	* \phi(1) &=& \phi'(1), \\
	* \psi(1) &=& \psi'(1), \\
	* \Delta E(1) &=& \frac{1}{3\pi}\sqrt{\frac{T}{2}}\psi'(0) - \Delta E_0,
	* @f}
	* Finally, @f$ \Delta E = \Delta E(0) @f$.
	*/
    void calculateEnergy();
	/**
	* @brief Implemented method for calculating correction to entropy.
	* @details The calculation goes through solving the ODE with parameters
	* @f$ p_0 = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ p_1 = \pi^{-3}VT @f$.
	* Correction to Thomas-Fermi entropy is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& p_0 x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  \frac{d^2\psi}{dx^2} &=& p_0 \left[
	*		       \frac12 \psi(x) I_{-1/2}\left(\frac{\phi(x)}{x}\right) + xY'(x)
	*           \right], \\
	*  \Delta S'(x) &=& -p_1 x \left[
	*						\frac12 \psi(x) I_{1/2}\left(\frac{\phi(x)}{x}\right) + 2xY(x)
	*					\right], \\
	* \phi(1) &=& \phi'(1), \\
	* \psi(1) &=& \psi'(1), \\
	* \Delta S(1) &=& \frac{1}{3\pi\sqrt{2T}}\psi'(0),
	* @f}
	* Finally, @f$ \Delta S = \Delta S(0) @f$.
	*/
    void calculateEntropy();
	/**
	* @brief Implemented method for calculating correction to chemical potential.
	* @details The formula for calculation:
	* @f[
	*	\Delta \mu = \frac{1}{3\pi}\sqrt{\frac{T}{2}}
	*    \left[\frac12 I_{-1/2}(\phi(1)) + \psi(1)\right].
	* @f]
	*/
    void calculateChemPotential();
    
	/**
	* @brief Implemented method for calculating thermal pressure in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*	\Delta P_T = \Delta P - \left. \Delta P \right|_{T = 0}.
	* @f]
	*/
    void calculateThermalP();
	/**
	* @brief Implemented method for calculating thermal pressure in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	* \Delta E_T = \Delta E - \left.\Delta E\right|_{T = 0}.
	* @f]
	*/
    void calculateThermalE();
	/**
	* @brief Implemented method for calculating thermal pressure in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*	\Delta S_T = \Delta S - \left.\Delta S\right|_{T = 0}.
	* @f]
	*/
    void calculateThermalS();
	/**
	* @brief Implemented method for calculating thermal pressure in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*   \Delta \mu_T = \Delta \mu - \left.\Delta \mu\right|_{T = 0}.
	* @f]
	*/
    void calculateThermalCP();
	/**
	* @brief Implemented method which sets all necessary parameters for Thomas-Fermi potential and corrections to it.
	*/
    void getReady();

private:

    TFPotential* potential;
    TFPotential* coldPotential;
    TFPotentialCorrection* potentialCorrection;
    TFPotentialCorrection* coldPotentialCorrection;
    ODESolver<TFCorrection::Energy>* energySolver;
    ODESolver<TFCorrection::Entropy>* entropySolver;

    Temperature& makeChargeless(Temperature& T);
    Volume& makeChargeless(Volume& V);
    Pressure& makeChargeless(Pressure& P);
    Energy& makeChargeless(Energy& E);
    Entropy& makeChargeless(Entropy& E);
    ChemicalPotential& makeChargeless(ChemicalPotential& M);

    Temperature& makeChargeful(Temperature& T);
    Volume& makeChargeful(Volume& V);
    Pressure& makeChargeful(Pressure& P);
    Energy& makeChargeful(Energy& E);
    Entropy& makeChargeful(Entropy& E);
    ChemicalPotential& makeChargeful(ChemicalPotential& M);

    double coldP();
    double coldE();
    double coldS();
    double coldM();
};

