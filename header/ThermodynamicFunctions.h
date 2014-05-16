#pragma once
#include <math.h>
#include "Units.h"
/**
 * @brief This class provides interface for thermodynamic functions.
 * @details Contains value, unit and scale (linear or logarithmic) of the thermodynamic function. 
 */
class ThermodynamicFunction {
public:
   /**
	* @brief Constructor of the thermodynamic function.
	*/
    ThermodynamicFunction();
   /**
	* @brief Destructor of the thermodynamic function.
	*/
    ~ThermodynamicFunction();
   /**
	* @brief Get current value of the thermodynamic function.
	*/
    double getValue();
   /**
	* @brief Get current scale of the thermodynamic function.
	*/
    Scaling getCurrentScale();
   /**
	* @brief Get current unit of the thermodynamic function.
	*/
    Unit getCurrentUnit();
protected:
    double defaultValue;
    double currentValue;
    Unit unit;
    Scaling scale;
};
/**
* @brief This class implements thermodynamic function interface for the temperature.
*/
class Temperature : public ThermodynamicFunction {
public:
	/**
	* @brief A default constructor of thermodynamic function.
	*/
    Temperature();
	/**
	* @brief An initialization by other temperature function.
	*/
    Temperature(const Temperature& T);
	/**
	* @brief Transform current units to new units.
	* @param unit New unit: Atomic, eV or K.
	* @return Reference to this object.
	*/
    Temperature& transformToUnits(Unit unit);
	/**
	* @brief Transform current scale to new scale.
	* @param scale New scale: lin or log.
	* @return Reference to this object.
	*/
    Temperature& transformToScale(Scaling scale);
	/**
	* @brief Operator '=' for temperature.
	* @return Reference to this object.
	*/
    Temperature& operator=(const Temperature& T);
	/**
	* @brief Sets a new value.
	* @param value A new value.
	* @param unit Unit type of a new value. Atomic by default.
	* @param scale Scaling of a new value. Linear by default.
	*/
    void setValue(double value, Unit unit = Atomic, Scaling scale = lin);
 };
/**
* @brief This class implements thermodynamic function interface for the volume.
*/
class Volume : public ThermodynamicFunction {
    public:
		/**
		* @brief A default constructor of thermodynamic function.
		*/
        Volume();
		/**
		* @brief An initialization by other volume function.
		*/
        Volume(const Volume& V);
		/**
		* @brief Transform current units to new units.
		* @param unit New unit: Atomic, cmc or mc.
		* @return Reference to this object.
		*/
        Volume& transformToUnits(Unit unit);
		/**
		* @brief Transform current scale to new scale.
		* @param scale New scale: lin or log.
		* @return Reference to this object.
		*/
        Volume& transformToScale(Scaling scale);
		/**
		* @brief Operator '=' for volume.
		* @return Reference to this object.
		*/
        Volume& operator=(const Volume& V);
		/**
		* @brief Sets a new value.
		* @param value A new value.
		* @param unit Unit type of a new value. Atomic by default.
		* @param scale Scaling of a new value. Linear by default.
		*/
        void setValue(double value, Unit unit = Atomic, Scaling scale = lin);
};
/**
* @brief This class implements thermodynamic function interface for the density.
*/
class Density : public ThermodynamicFunction {
public:
	/**
	* @brief A default constructor of thermodynamic function.
	*/
    Density();
	/**
	* @brief An initialization by other density function.
	*/
    Density(const Density& D);
	/**
	* @brief Transform current units to new units.
	* @param unit New unit: Atomic, gOverCmc or kgOverMc.
	* @return Reference to this object.
	*/
    Density& transformToUnits(Unit unit);
	/**
	* @brief Transform current scale to new scale.
	* @param scale New scale: lin or log.
	* @return Reference to this object.
	*/
    Density& transformToScale(Scaling scale);
	/**
	* @brief Operator '=' for density.
	* @return Reference to this object.
	*/
    Density& operator=(const Density& D);
	/**
	* @brief Sets a new value.
	* @param value A new value.
	* @param unit Unit type of a new value. gOverCmc by default.
	* @param scale Scaling of a new value. Linear by default.
	*/
    void setValue(double value, Unit unit = gOverCmc, Scaling scale = lin);
};
/**
* @brief This class implements thermodynamic function interface for the concentration.
*/
class Concentration : public ThermodynamicFunction {
public:
	/**
	* @brief A default constructor of thermodynamic function.
	*/
    Concentration();
	/**
	* @brief An initialization by other concentration function.
	*/
    Concentration(const Concentration& C);
	/**
	* @brief Transform current units to new units.
	* @param unit New unit: Atomic, oneOverCmc or oneOverMc.
	* @return Reference to this object.
	*/
    Concentration& transformToUnits(Unit unit);
	/**
	* @brief Transform current scale to new scale.
	* @param scale New scale: lin or log.
	* @return Reference to this object.
	*/
    Concentration& transformToScale(Scaling scale);
	/**
	* @brief Operator '=' for concentration.
	* @return Reference to this object.
	*/
    Concentration& operator=(const Concentration& C);
	/**
	* @brief Sets a new value.
	* @param value A new value.
	* @param unit Unit type of a new value. Atomic by default.
	* @param scale Scaling of a new value. Linear by default.
	*/
    void setValue(double value, Unit unit = Atomic, Scaling scale = lin);
};
/**
* @brief This class implements thermodynamic function interface for the pressure.
*/
class Pressure : public ThermodynamicFunction {
public:
	/**
	* @brief A default constructor of thermodynamic function.
	*/
    Pressure();
	/**
	* @brief An initialization by other pressure function.
	*/
    Pressure(const Pressure& P);
	/**
	* @brief Transform current units to new units.
	* @param unit New unit: Atomic, Pa, GPa or MBar.
	* @return Reference to this object.
	*/
    Pressure& transformToUnits(Unit unit);
	/**
	* @brief Transform current scale to new scale.
	* @param scale New scale: lin or log.
	* @return Reference to this object.
	*/
    Pressure& transformToScale(Scaling scale);
	/**
	* @brief Operator '=' for pressure.
	* @return Reference to this object.
	*/
    Pressure& operator=(const Pressure& P);
	/**
	* @brief Sets a new value.
	* @param value A new value.
	* @param unit Unit type of a new value. Atomic by default.
	* @param scale Scaling of a new value. Linear by default.
	*/
    void setValue(double value, Unit unit = Atomic, Scaling scale = lin);
};
/**
* @brief This class implements thermodynamic function interface for the energy.
*/
class Energy : public ThermodynamicFunction {
public:
	/**
	* @brief A default constructor of thermodynamic function.
	*/
    Energy();
	/**
	* @brief An initialization by other energy function.
	*/
    Energy(const Energy& E);
	/**
	* @brief Transform current units to new units.
	* @param unit New unit: Atomic or eV.
	* @return Reference to this object.
	*/
    Energy& transformToUnits(Unit unit);
	/**
	* @brief Transform current scale to new scale.
	* @param scale New scale: lin or log.
	* @return Reference to this object.
	*/
    Energy& transformToScale(Scaling scale);
	/**
	* @brief Operator '=' for energy.
	* @return Reference to this object.
	*/
    Energy& operator=(const Energy& E);
	/**
	* @brief Sets a new value.
	* @param value A new value.
	* @param unit Unit type of a new value. Atomic by default.
	* @param scale Scaling of a new value. Linear by default.
	*/
    void setValue(double value, Unit unit = Atomic, Scaling scale = lin);
};
/**
* @brief This class implements thermodynamic function interface for the entropy.
*/
class Entropy : public ThermodynamicFunction {
public:
	/**
	* @brief A default constructor of thermodynamic function.
	*/
    Entropy();
	/**
	* @brief An initialization by other energy function.
	*/
    Entropy(const Entropy& S);
	/**
	* @brief Transform current units to new units.
	* @param unit New unit: only atomic.
	* @return Reference to this object.
	*/
    Entropy& transformToUnits(Unit unit);
	/**
	* @brief Transform current scale to new scale.
	* @param scale New scale: lin or log.
	* @return Reference to this object.
	*/
    Entropy& transformToScale(Scaling scale);
	/**
	* @brief Operator '=' for entropy.
	* @return Reference to this object.
	*/
    Entropy& operator=(const Entropy& S);
	/**
	* @brief Sets a new value.
	* @param value A new value.
	* @param unit Unit type of a new value. Atomic by default.
	* @param scale Scaling of a new value. Linear by default.
	*/
    void setValue(double value, Unit unit = Atomic, Scaling scale = lin);
};
/**
* @brief This class implements thermodynamic function interface for the chemical potential.
*/
class ChemicalPotential : public ThermodynamicFunction {
public:
	/**
	* @brief A default constructor of thermodynamic function.
	*/
    ChemicalPotential();
	/**
	* @brief An initialization by other chemical potential function.
	*/
    ChemicalPotential(const ChemicalPotential& M);
	/**
	* @brief Transform current units to new units.
	* @param unit New unit: only atomic.
	* @return Reference to this object.
	*/
    ChemicalPotential& transformToUnits(Unit unit);
	/**
	* @brief Transform current scale to new scale.
	* @param scale New scale: lin or log.
	* @return Reference to this object.
	*/
    ChemicalPotential& transformToScale(Scaling scale);
	/**
	* @brief Operator '=' for chemical potential.
	* @return Reference to this object.
	*/
    ChemicalPotential& operator=(const ChemicalPotential& M);
	/**
	* @brief Sets a new value.
	* @param value A new value.
	* @param unit Unit type of a new value. Atomic by default.
	* @param scale Scaling of a new value. Linear by default.
	*/
    void setValue(double value, Unit unit = Atomic, Scaling scale = lin);
};
