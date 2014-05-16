#pragma once

#include "ThermodynamicFunctions.h"
#include <string>
/**
 * @brief A class provides interface for models.
 * @details The idea is to set parameters V and T in atomic units, precision and get 
 * values of the thermodynamic functions: pressure, energy, entropy and chemical potential and
 * correspoinding thermal parts.
 */
class Model {

public:
   /**
	* @brief A constructor of a model.
	* @param Z An atomic number of element. Default is 1.
	*/
    Model(Charge Z = 1.0);
   /**
	* @brief A destructor of a model.
	*/
    ~Model();
   /**
	* @brief Sets precision.
	* @param precision Acceptable relative error level.
	*/
    void setPrecision(const double precision);
   /**
	* @brief Sets parameters.
	* @param V Volume of the 1 atom.
	* @param T Temperature.
	*/
    void setParameters(Volume V, Temperature T);
   /**
	* @brief Sets atomic number.
	* @param Z Atomic number.
	*/
    void setAtomicNumber(Charge Z);
   /**
	* @brief Get name of model.
	* @return String name of the model.
	*/
    virtual std::string getName();
   /**
	* @brief Get current temperature.
	*/
    Temperature getTemperature(); 
   /**
	* @brief Get current volume.
	*/
    Volume getVolume();
	/**
	* @brief Get current pressure.
	*/
    Pressure getPressure();
	/**
	* @brief Get current energy.
	*/
    Energy getEnergy();
	/**
	* @brief Get current entropy.
	*/
    Entropy getEntropy();
	/**
	* @brief Get current chemical potential.
	*/
    ChemicalPotential getChemicalPotential();
	/**
	* @brief Get current thermal pressure.
	*/
    Pressure getThermalPressure();
	/**
	* @brief Get current thermal energy.
	*/
    Energy getThermalEnergy();
	/**
	* @brief Get current thermal entropy.
	*/
    Entropy getThermalEntropy();
	/**
	* @brief Get current thermal chemical potential.
	*/
    ChemicalPotential getThermalChemicalPotential();

protected:

	/**
	* @brief Current temperature in model.
	*/
    Temperature T;
	/**
	* @brief Current volume in model.
	*/
    Volume V;
	/**
	* @brief Current pressure in model.
	*/
    Pressure P;
	/**
	* @brief Current energy in model.
	*/
    Energy E;
	/**
	* @brief Current entropy in model.
	*/
    Entropy S;
	/**
	* @brief Current chemical potential in model.
	*/
    ChemicalPotential M;
	/**
	* @brief Current thermal pressure in model.
	*/
    Pressure Pthermal;
	/**
	* @brief Current thermal energy in model.
	*/
    Energy Ethermal;
	/**
	* @brief Current thermal entropy in model.
	*/
    Entropy Sthermal;
	/**
	* @brief Current thermal chemical potential in model.
	*/
    ChemicalPotential Mthermal;
	/**
	* @brief Atomic number of element in the model.
	*/
    Charge Z;
	/**
	* @brief Full virtual method for calculating pressure.
	*/
    virtual void calculatePressure();
	/**
	* @brief Full virtual method for calculating energy.
	*/
    virtual void calculateEnergy();
	/**
	* @brief Full virtual method for calculating entropy.
	*/
    virtual void calculateEntropy();
	/**
	* @brief Full virtual method for calculating chemical potential.
	*/
    virtual void calculateChemPotential();
	/**
	* @brief Full virtual method for calculating thermal pressure.
	*/
    virtual void calculateThermalP();
	/**
	* @brief Full virtual method for calculating thermal energy.
	*/
    virtual void calculateThermalE();
	/**
	* @brief Full virtual method for calculating thermal entropy.
	*/
    virtual void calculateThermalS();
	/**
	* @brief Full virtual method for calculating thermal chemical potential.
	*/
    virtual void calculateThermalCP();
	/**
	* @brief Full virtual method which sets all the parameters before calculating.
	*/
    virtual void getReady();
	/**
	* @brief Bool variable of state for calculating pressure.
	* @details Is set to true if pressure has already calculated for the current parameters.
	*/
    bool calculatedPressure;
	/**
	* @brief Bool variable of state for calculating energy.
	* @details Is set to true if has already calculated for the current parameters.
	*/
    bool calculatedEnergy;
	/**
	* @brief Bool variable of state for calculating entropy.
	* @details Is set to true if entropy has already calculated for the current parameters.
	*/
    bool calculatedEntropy;
	/**
	* @brief Bool variable of state for calculating chemical potential.
	* @details Is set to true if chemical potential has already calculated for the current parameters.
	*/
    bool calculatedChemPot;
	/**
	* @brief Bool variable of state for calculating thermal pressure.
	* @details Is set to true if thermal pressure has already calculated for the current parameters.
	*/
    bool calculatedThermalP;
	/**
	* @brief Bool variable of state for calculating thermal energy.
	* @details Is set to true if thermal energy has already calculated for current parameters.
	*/
    bool calculatedThermalE;
	/**
    * @brief Bool variable of state for calculating thermal entropy.
	* @details Is set to true if entropy has already calculated for current parameters.
	*/
    bool calculatedThermalS;
	/**
	* @brief Bool variable of state for calculating thermal chemical potential.
	* @details Is set to true if thermal chemical potential has already calculated for current parameters.
	*/
    bool calculatedThermalCP;
	/**
	* @brief Bool variable of state for all the parameters for ODE and Thomas-Fermi potential are set.
	*/
    bool readyToCalculate;
	/**
	* @brief Bool variable which indicates that parameters V and T have been set.
	*/
    bool parameterSet;
	/**
	* @brief Bool variable which indicates that new value of precision have been set.
	*/
    bool precisionSet;
	/**
	* @brief Value of acceptable precision.
	*/
    double eps;
	/**
	* @brief Value of temperature which corresponds to zero.
	*/
    const double coldT;
};