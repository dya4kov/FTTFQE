#pragma once
#include "ThomasFermi.h"
#include "ThomasFermiCorrection.h"
#include "Units.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cctype>

typedef std::string string;
typedef double Data;
/**
* @brief Class for reading input file and generating output.
* @details Format of the input file:
*    - first line should indicate the names of models for calculating:
*	   -# "model: TF;"\n
*           Calculates only Thomas-Fermi functions.
*	   -# "model: TF, TFC;"\n
*	       Calculates both Thomas-Fermi functions and corrections to them.
*	 - next lines should indicate:
*	   -# "start temperature: [Unit] (Scaling) value;"
*	   -# "start volume/density/concentration: [Unit] (Scaling) value;"
*	   -# "end temperature: [Unit] (Scaling) value;"
*	   -# "end volume/density/concentration: [Unit] (Scaling) value;"
*      -# "number of points in T: value;"
*      -# "number of points in V/D/C: value;"
*      -# "atomic number: value;"
*      -# "atomic mass: value;"\n
*         Is used only with density.
*	   -# "calculate P/PT/E/ET/S/ST/M/MT: [Unit] (Scaling);"
*      -# "presicion: value;"\n
*		  By default is 1e-5.
*      -# Other options:\n
*		  "calculate PVT: [-] lin;" - calculates thermal equation of state;\n
*		  "calculate caloric: [-] lin;" - calculates caloric equation of state;\n
*		  "calculate meanCharge: [-] lin;" - calculates effective electron charge 
*		   at the bound of an atom;\n
*		  "calculate region of validity: [-] lin;" - calculates region of validity of the
*         Thomas-Fermi model.\n
*		  "calculate region of validity thermal: [-] lin;" - calculates region of validity 
*		  of the thermal contribution to the Thomas-Fermi model.
*
*    The order of lines doesn't matter. ';' is necessary.
*/
class Interface
{
public:
	/**
	* @brief Constructor of this class.
	*/
    Interface(void);
	/**
	* @brief Destructor of this class.
	*/
    ~Interface(void);
	/**
	* @brief Read and parse file.
	* @param filename Name of the input file.
	*/
    void Read(string filename);
	/**
	* @brief Generating the output file.
	*/
    void Print();

private:

    enum { 
        TF, 
        TFC, 
        Size, 
        Error 
    };

    enum ParameterType {
        undefinedType,
        startT,
        endT,
        startV,
        endV,
        startD,
        endD,
        startC,
        endC,
        charge,
        mass,
        pointsT,
        pointsV,
        pointsD,
        pointsC,
        pressure,
        action,
        energy,
        entropy,
        Pthermal,
        Ethermal,
        Sthermal,
        CPthermal,
        chemicalPotential,
        precision,
        meanCharge,
        PVT,
		caloric,
		validityBoundary,
		thermalValidityBoundary,
        final
    };

    struct ParameterData {
        Scaling scale;
        Unit unit;
        Data data;
    };

    string inputFilename;

    Model* models[Size];
    std::map<ParameterType, ParameterData> parameters;
    std::vector<Volume> V;
    std::vector<Temperature> T;
    std::vector<Density> D;
    std::vector<Concentration> C;
    std::ofstream out;
    int digitsToPrint;
    
    bool volumeAssigned;
    bool densityAssigned;
    bool concentrationAssigned;
    bool temperatureAssigned;
	bool precisionAssigned;

    int getID(string modelName) const;
    Model* initializeModel(string modelName) const;

    ParameterType defineType(string parameterName) const;
    Scaling defineScale(string scale) const;
    Unit defineUnit(string unit) const;

    bool assignVolume();
    bool assignTemperature();
    bool assignDensity();
    bool assignConcentration();
    bool assignPrecision();

    void printValue(double value);
    void printFirstLine();
    void printPressure(Model* model, Unit unit, Scaling scale);
    void printEnergy(Model* model, Unit unit, Scaling scale);
    void printEntropy(Model* model, Unit unit, Scaling scale);
    void printChemicalPotential(Model* model, Unit unit, Scaling scale);
    void printThermalPressure(Model* model, Unit unit, Scaling scale);
    void printThermalEnergy(Model* model, Unit unit, Scaling scale);
    void printThermalEntropy(Model* model, Unit unit, Scaling scale);
    void printThermalChemicalPotential(Model* model, Unit unit, Scaling scale);
    void printMeanCharge();
    void printPVT();
	void printCaloric();
	void printValidityBoundary();
	void printThermalValidityBoundary();
};

