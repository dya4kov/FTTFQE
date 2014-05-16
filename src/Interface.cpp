#include "../header/Interface.h"

Interface::Interface(void) {
    for (int i = 0; i < Size; ++i) {
        models[i] = NULL;
    }
}

Interface::~Interface(void) {
    for (int i = 0; i < Size; ++i) {
        if (models[i] != NULL) delete models[i];
    }
}

int Interface::getID(string modelName) const {
    if (!modelName.compare("TF")) return TF;
    if (!modelName.compare("TFC")) return TFC;
    return Error;
}

Model* Interface::initializeModel(string modelName) const {
    if (!modelName.compare("TF")) return new ThomasFermi;
    if (!modelName.compare("TFC")) return new ThomasFermiCorrection;
    return NULL;
}

Interface::ParameterType Interface::defineType(string parameterName) const {
    if (!parameterName.compare("start temperature")) return startT;
    if (!parameterName.compare("start volume")) return startV;
    if (!parameterName.compare("start density")) return startD;
    if (!parameterName.compare("start concentration")) return startC;
    if (!parameterName.compare("end temperature")) return endT;
    if (!parameterName.compare("end volume")) return endV;
    if (!parameterName.compare("end density")) return endD;
    if (!parameterName.compare("end concentration")) return endC;
    if (!parameterName.compare("atomic number")) return charge;
    if (!parameterName.compare("atomic mass")) return mass;
    if (!parameterName.compare("number of points in T")) return pointsT;
    if (!parameterName.compare("number of points in V")) return pointsV;
    if (!parameterName.compare("number of points in D")) return pointsD;
    if (!parameterName.compare("number of points in C")) return pointsC;
    if (!parameterName.compare("calculate P")) return pressure;
    if (!parameterName.compare("calculate E")) return energy;
    if (!parameterName.compare("calculate S")) return entropy;
    if (!parameterName.compare("calculate M")) return chemicalPotential;
    if (!parameterName.compare("calculate PT")) return Pthermal;
    if (!parameterName.compare("calculate ET")) return Ethermal;
    if (!parameterName.compare("calculate ST")) return Sthermal;
    if (!parameterName.compare("calculate MT")) return CPthermal;
    if (!parameterName.compare("calculate meanCharge")) return meanCharge;
	if (!parameterName.compare("calculate caloric")) return caloric;
    if (!parameterName.compare("precision")) return precision;
    if (!parameterName.compare("calculate PVT")) return PVT;
	if (!parameterName.compare("region of validity")) return validityBoundary;
	if (!parameterName.compare("region of validity thermal")) return thermalValidityBoundary;
    return undefinedType;
}

Unit Interface::defineUnit(string unit) const {
    if (!unit.compare("-")) return dimensionless;
    if (!unit.compare("eV")) return eV;
    if (!unit.compare("hartree")) return Atomic;
    if (!unit.compare("gcc")) return gcc;
    if (!unit.compare("K")) return Kelvin;
    if (!unit.compare("cmc")) return cmc;
    if (!unit.compare("gOverCmc")) return gOverCmc;
    if (!unit.compare("kgOverMc")) return kgOverMc;
    if (!unit.compare("oneOverCmc")) return oneOverCmc;
    if (!unit.compare("oneOverMc")) return oneOverMc;
    if (!unit.compare("GPa")) return GPa;
    if (!unit.compare("MBar")) return MBar;
    if (!unit.compare("mc")) return mc;
    if (!unit.compare("Pa")) return Pa;

    return undefinedUnit;
}

Scaling Interface::defineScale(string scale) const {
    if (!scale.compare("lin")) return lin;
    if (!scale.compare("lg")) return lg;

    return lin;
}

void Interface::Read(string filename) {
    std::cout << "reading of input file...   ";
    const int bufsize = 200;
    char buffer[bufsize];

    const char endLine = ';';
    const char startUnit = '[';
    const char endUnit = ']';
    const char startScale = '(';
    const char endScale = ')';
    const char endType = ':';
    const char modelDelimeter = ',';

    string type;
    string unit;
    string value;
    string scale;
    string model;
    string currentLine;

    string::iterator begin;
    string::iterator end;

    ParameterData currentData;

    inputFilename = filename;
    std::ifstream inputData(filename.c_str(), std::ios::in);

    if (!inputData.is_open()) {
        std::cout << "cannot find file: " << filename << std::endl;
        exit(0);
    }

    //reading first line to define models
    inputData.getline(buffer, bufsize, endType);
    inputData.getline(buffer, bufsize, endLine); 
    currentLine.assign(buffer);
    currentLine += endLine;
    begin = currentLine.begin();
    end = currentLine.begin();

    while (begin != currentLine.end()) {

        while (!isalpha(*begin) && *begin != endLine) ++begin;
        if (begin == currentLine.end()) break;
        while (*end != endLine && *end != modelDelimeter) ++end;

        model.assign(begin, end);

        if (models[getID(model)] == NULL) {
                models[getID(model)] = initializeModel(model);
        }
        begin = end + 1;
        end = end + 1;
    }

    //define other parameters
    inputData.getline(buffer, bufsize, endType);
    while (!inputData.eof()) {

        //assign parameter type
        currentLine.assign(buffer);
        begin = currentLine.begin();
        end = currentLine.end() - 1;
        while (!isalpha(*begin) && begin != currentLine.end()) ++begin;
        while (!isalpha(*end) && end != currentLine.begin()) --end;
        type.assign(begin, end + 1);

        inputData.getline(buffer, bufsize, endLine);
        currentLine.assign(buffer);
        //assign unit type
        begin = currentLine.begin();
        end = currentLine.end() - 1;

        while(*begin != startUnit && begin != currentLine.end() - 1) ++begin;
        while(*end != endUnit && end != currentLine.begin()) --end;
        if (begin == currentLine.end() - 1) unit = "-";
        else unit.assign(begin + 1, end);
        //assign scaling
        begin = currentLine.begin();
        end = currentLine.end() - 1;

        while(*begin != startScale && begin != currentLine.end() - 1) ++begin;
        while(*end != endScale && end != currentLine.begin()) --end;
        if (begin == currentLine.end() - 1) scale = "lin";
        else scale.assign(begin + 1, end);
        //assign value
        end = currentLine.end() - 1;

        while(!isdigit(*end) && end != currentLine.begin()) --end;
        if (end != currentLine.begin()) {
            begin = end;
            while(!isspace(*begin) && begin != currentLine.begin()) --begin;
            begin += 1;
            value.assign(begin, end + 1);
        }
        else {
            value = "0";
        }

        currentData.scale = defineScale(scale);
        currentData.unit = defineUnit(unit);
        currentData.data = atof(value.c_str());  

        parameters[defineType(type)] = currentData;
        inputData.getline(buffer, bufsize, endType);
    }

    inputData.close();
    std::cout << "completed" << std::endl;
}

bool Interface::assignVolume() {
    std::vector<Volume>::iterator vIter;
    std::map<ParameterType, ParameterData>::iterator paramIter;
    double vStep;
    int numberOfPoints;
    Volume firstV;
    Volume lastV;
    paramIter = parameters.find(startV);
    if (paramIter != parameters.end()) {
        firstV.setValue(paramIter->second.data, 
        paramIter->second.unit,
        paramIter->second.scale);
        V.push_back(firstV);
        paramIter = parameters.find(endV);
        if (paramIter != parameters.end()) {
            lastV.setValue(paramIter->second.data, 
                           paramIter->second.unit,
                           paramIter->second.scale);
            paramIter = parameters.find(pointsV);
            if (paramIter != parameters.end()) {
                numberOfPoints = static_cast<int>(floor(paramIter->second.data));
                    vStep = (lastV.getValue() - firstV.getValue())/
                            ((numberOfPoints >= 2 ? numberOfPoints : 2) - 1);
                for (int i = 1; i < numberOfPoints; ++i) {
                    lastV.setValue(firstV.getValue() + i*vStep,
                    firstV.getCurrentUnit(),
                    firstV.getCurrentScale());
                    V.push_back(lastV);
                }
            }
        }
    }
    else {
        return volumeAssigned = false;
    }
    return volumeAssigned = true;
}

bool Interface::assignTemperature() {
    std::vector<Temperature>::iterator tIter;
    std::map<ParameterType, ParameterData>::iterator paramIter;
    double tStep;
    int numberOfPoints;
    Temperature firstT;
    Temperature lastT;
    paramIter = parameters.find(startT);
    if (paramIter != parameters.end()) {
        firstT.setValue(paramIter->second.data, 
                        paramIter->second.unit,
                        paramIter->second.scale);
        T.push_back(firstT);
        paramIter = parameters.find(endT);
        if (paramIter != parameters.end()) {
            lastT.setValue(paramIter->second.data, 
                           paramIter->second.unit,
                           paramIter->second.scale);
            paramIter = parameters.find(pointsT);
            if (paramIter != parameters.end()) {
                numberOfPoints = static_cast<int>(floor(paramIter->second.data));
                tStep = (lastT.getValue() - firstT.getValue())/
                        ((numberOfPoints >= 2 ? numberOfPoints : 2) - 1);
                for (int i = 1; i < numberOfPoints; ++i) {
                    lastT.setValue(firstT.getValue() + i*tStep,
                                   firstT.getCurrentUnit(),
                                   firstT.getCurrentScale());
                    T.push_back(lastT);
                }
            }
        }
    }
    else {
        return temperatureAssigned = false;
    }
    return temperatureAssigned = true;
}

bool Interface::assignDensity() {
    std::vector<Density>::iterator dIter;
    std::vector<Volume>::iterator vIter;
    std::map<ParameterType, ParameterData>::iterator paramIter;
    double dStep;
    int numberOfPoints;
    double A;
    Density firstD;
    Density lastD;
    Volume currentV;
    Density currentD;
    paramIter = parameters.find(startD);
    if (paramIter != parameters.end()) {
        firstD.setValue(paramIter->second.data, 
        paramIter->second.unit,
        paramIter->second.scale);
        //assign mass
        paramIter = parameters.find(mass);
        if (paramIter != parameters.end()) {
            A = paramIter->second.data;
        }
        else {
            A = 1.0;
        }
        currentD = firstD;
        currentD.transformToScale(lin).transformToUnits(gOverCmc);
        currentV.setValue(A/(Avogadro*1.4818e-25)/currentD.getValue());
        D.push_back(firstD);
        V.push_back(currentV);
        paramIter = parameters.find(endD);
        if (paramIter != parameters.end()) {
            lastD.setValue(paramIter->second.data, 
                           paramIter->second.unit,
                           paramIter->second.scale);
            paramIter = parameters.find(pointsD);
            if (paramIter != parameters.end()) {
                numberOfPoints = static_cast<int>(floor(paramIter->second.data));
                    dStep = (lastD.getValue() - firstD.getValue())/
                            ((numberOfPoints >= 2 ? numberOfPoints : 2) - 1);
                for (int i = 1; i < numberOfPoints; ++i) {
                    lastD.setValue(firstD.getValue() + i*dStep,
                    firstD.getCurrentUnit(),
                    firstD.getCurrentScale());
                    currentD = lastD;
                    currentD.transformToScale(lin).transformToUnits(gOverCmc);
                    currentV.setValue(A/(Avogadro*1.4818e-25)/currentD.getValue());
                    D.push_back(lastD);
                    V.push_back(currentV);
                }
            }
        }
    }
    else {
        return densityAssigned = false;
    }
    return densityAssigned = true;
}

bool Interface::assignConcentration() {
    std::vector<Concentration>::iterator cIter;
    std::vector<Volume>::iterator vIter;
    std::map<ParameterType, ParameterData>::iterator paramIter;
    double cStep;
    int numberOfPoints;
    double mass;
    Concentration firstC;
    Concentration lastC;
    Volume currentV;
    Concentration currentC;
    paramIter = parameters.find(startC);
    if (paramIter != parameters.end()) {
        firstC.setValue(paramIter->second.data, 
        paramIter->second.unit,
        paramIter->second.scale);
        currentC = firstC;
        currentC.transformToScale(lin).transformToUnits(Atomic);
        currentV.setValue(1.0/currentC.getValue());
        C.push_back(firstC);
        V.push_back(currentV);
        paramIter = parameters.find(endC);
        if (paramIter != parameters.end()) {
            lastC.setValue(paramIter->second.data, 
                           paramIter->second.unit,
                           paramIter->second.scale);
            paramIter = parameters.find(pointsC);
            if (paramIter != parameters.end()) {
                numberOfPoints = static_cast<int>(floor(paramIter->second.data));
                    cStep = (lastC.getValue() - firstC.getValue())/
                            ((numberOfPoints >= 2 ? numberOfPoints : 2) - 1);
                for (int i = 1; i < numberOfPoints; ++i) {
                    lastC.setValue(firstC.getValue() + i*cStep,
                                   firstC.getCurrentUnit(),
                                   firstC.getCurrentScale());
                    currentC = lastC;
                    currentC.transformToScale(lin).transformToUnits(Atomic);
                    currentV.setValue(1.0/currentC.getValue());
                    C.push_back(lastC);
                    V.push_back(currentV);
                }
            }
        }
    }
    else {
        return concentrationAssigned = false;
    }
    return concentrationAssigned = true;
}

bool Interface::assignPrecision() {
    std::map<ParameterType, ParameterData>::iterator paramIter;
    double epsilon;
    digitsToPrint = 0;
    paramIter = parameters.find(precision);
    if (paramIter != parameters.end()) {
    epsilon = paramIter->second.data;
        while (epsilon < 1.0) {
                epsilon *= 10;
                digitsToPrint++;
        }
		return precisionAssigned = true;
    }
    else {
        digitsToPrint = 5;
		return precisionAssigned = false;
    }
}

void Interface::Print() {

    if (!assignVolume()) 
        if (!assignDensity())
            if (!assignConcentration()) {
                std::cerr << "cannot find volume, "
                          << "density or concentration points"
                          << std::endl;
                exit(0);
            }
                
    if (!assignTemperature()) {
        std::cerr << "cannot find temperature points"
                  << std::endl;
        exit(0);
    }
	if (!assignPrecision()) {
		std::cout << "precision set to default value 1e-5"
			      << std::endl;
	}

    std::map<ParameterType, ParameterData>::iterator paramIter;
    string outFileName = inputFilename + "_out.dat";
    out.open(outFileName.c_str(), std::ios::out);

    for (int i = 0; i < Size; ++i) {

        if (models[i] != NULL) {
            out << "results for model " << models[i]->getName() << ":" << std::endl << std::endl;
            std::cout << "starting calculation for model " << models[i]->getName() << ":" << std::endl;

            paramIter = parameters.find(charge);
            if (paramIter != parameters.end()) {
                models[i]->setAtomicNumber(paramIter->second.data);
            }
            else {
                models[i]->setAtomicNumber(1.0);
            }

            paramIter = parameters.find(precision);
            if (paramIter != parameters.end()) {
                models[i]->setPrecision(paramIter->second.data);
            }
            else {
                models[i]->setPrecision(1e-5);
            }

            paramIter = parameters.find(pressure);
            if (paramIter != parameters.end()) {
                    printPressure(models[i], 
                                  paramIter->second.unit,
                                  paramIter->second.scale);
            }
            paramIter = parameters.find(energy);
            if (paramIter != parameters.end()) {
                    printEnergy(models[i],
                                paramIter->second.unit,
                                paramIter->second.scale);
            }
            paramIter = parameters.find(entropy);
            if (paramIter != parameters.end()) {
                    printEntropy(models[i],
                                 paramIter->second.unit,
                                 paramIter->second.scale);
            }
            paramIter = parameters.find(chemicalPotential);
            if (paramIter != parameters.end()) {
                    printChemicalPotential(models[i],
                                           paramIter->second.unit,
                                           paramIter->second.scale);
            }
                    
            paramIter = parameters.find(Pthermal);
            if (paramIter != parameters.end()) {
                    printThermalPressure(models[i], 
                                  paramIter->second.unit,
                                  paramIter->second.scale);
            }
            paramIter = parameters.find(Ethermal);
            if (paramIter != parameters.end()) {
                    printThermalEnergy(models[i],
                                paramIter->second.unit,
                                paramIter->second.scale);
            }
            paramIter = parameters.find(Sthermal);
            if (paramIter != parameters.end()) {
                    printThermalEntropy(models[i],
                                 paramIter->second.unit,
                                 paramIter->second.scale);
            }
            paramIter = parameters.find(CPthermal);
            if (paramIter != parameters.end()) {
                    printThermalChemicalPotential(models[i],
                                           paramIter->second.unit,
                                           paramIter->second.scale);
            }
            
            std::cout << "calculation for model " << models[i]->getName() << " finished" << std::endl;
        }
        
    }
    
    paramIter = parameters.find(meanCharge);
    if (paramIter != parameters.end()) {
                printMeanCharge();
    }
    paramIter = parameters.find(PVT);
    if (paramIter != parameters.end()) {
                printPVT();
    }
	paramIter = parameters.find(validityBoundary);
    if (paramIter != parameters.end()) {
                printValidityBoundary();
    }
	paramIter = parameters.find(thermalValidityBoundary);
    if (paramIter != parameters.end()) {
                printThermalValidityBoundary();
    }
	paramIter = parameters.find(caloric);
    if (paramIter != parameters.end()) {
                printCaloric();
    }
	 
    out.close();
}

void Interface::printValue(double value) {
    out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
    out.width(20);
    out.precision(digitsToPrint);
    out.fill(' ');
    out << value;
}

void Interface::printFirstLine() {
    if (volumeAssigned) {
        std::vector<Volume>::iterator vIter;
        out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
        out.width(20);
        out.precision(digitsToPrint);
        out.fill(' ');
        out << "T\\V"; 
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            printValue(vIter->getValue());
        }
        out << std::endl;
    }
    else if (densityAssigned) {
        std::vector<Density>::iterator dIter;
        out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
        out.width(20);
        out.precision(digitsToPrint);
        out.fill(' ');
        out << "T\\Density"; 
        for (dIter = D.begin(); dIter != D.end(); ++dIter) {
            printValue(dIter->getValue());
        }
        out << std::endl;
    }
    else if (concentrationAssigned) {
        std::vector<Concentration>::iterator cIter;
        out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
        out.width(20);
        out.precision(digitsToPrint);
        out.fill(' ');
        out << "T\\Concentration"; 
        for (cIter = C.begin(); cIter != C.end(); ++cIter) {
            printValue(cIter->getValue());
        }
        out << std::endl;
    }
}

void Interface::printPressure(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating pressure: ";
    out << "pressure:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            model->setParameters(*vIter, *tIter);
            printValue(model->getPressure().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printThermalPressure(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating thermal part of pressure: ";
    out << "thermal part of pressure:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            model->setParameters(*vIter, *tIter);
            printValue(model->getThermalPressure().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printEnergy(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating energy: ";
    out << "energy:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            model->setParameters(*vIter, *tIter);
            printValue(model->getEnergy().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printThermalEnergy(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating thermal part of energy: ";
    out << "thermal part of energy:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            model->setParameters(*vIter, *tIter);
            printValue(model->getThermalEnergy().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printEntropy(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating entropy: ";
    out << "entropy:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
                model->setParameters(*vIter, *tIter);
                printValue(model->getEntropy().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printThermalEntropy(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating thermal part of entropy: ";
    out << "thermal part of entropy:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
                model->setParameters(*vIter, *tIter);
                printValue(model->getThermalEntropy().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printChemicalPotential(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating chemical potential: ";
    out << "chemical potential:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            model->setParameters(*vIter, *tIter);
            printValue(model->getChemicalPotential().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printThermalChemicalPotential(Model* model, Unit unit, Scaling scale) {
    std::cout << "calculating thermal part of chemical potential: ";
    out << "thermal part of chemical potential:" << std::endl << std::endl;

    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;

    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            model->setParameters(*vIter, *tIter);
            printValue(model->getThermalChemicalPotential().transformToUnits(unit).transformToScale(scale).getValue());
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printMeanCharge() {
    std::cout << "calculating MeanCharge: ";
    out << "meanCharge:" << std::endl << std::endl;
    ThomasFermi* tf = (ThomasFermi*) models[TF];
    ThomasFermiCorrection* tfc = (ThomasFermiCorrection*) models[TFC];
    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;
    double meanCharge;
    double M1 = 0;
    double M2 = 0;
    double Temp;
    
    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            if (tf != NULL) { 
                tf->setParameters(*vIter, *tIter);
                M1 = tf->getChemicalPotential().getValue();
            }
            if (tfc != NULL) {
                tfc->setParameters(*vIter, *tIter);
                M2 = tfc->getChemicalPotential().getValue();
            }
            Temp = tIter->transformToUnits(Atomic).transformToScale(lin).getValue();
            meanCharge = pow(Temp, 1.5)*sqrt(2.0)/M_PI/M_PI*vIter->getValue()
            * FermiDirac::Half::function((M1 + M2)/Temp);
            printValue(meanCharge);
            
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printPVT() {
    std::cout << "calculating PVT: ";
    out << "PVT:" << std::endl << std::endl;
    ThomasFermi* tf = (ThomasFermi*) models[TF];
    ThomasFermiCorrection* tfc = (ThomasFermiCorrection*) models[TFC];
    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;
    double P1 = 0;
    double P2 = 0;
    double Temp;
    double Vol;
    double PVTvalue;
    
    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            if (tf != NULL) { 
                tf->setParameters(*vIter, *tIter);
                P1 = tf->getPressure().getValue();
            }
            if (tfc != NULL) {
                tfc->setParameters(*vIter, *tIter);
                P2 = tfc->getPressure().getValue();
            }
            Temp = tIter->transformToUnits(Atomic).transformToScale(lin).getValue();
            Vol = vIter->transformToUnits(Atomic).transformToScale(lin).getValue();
            PVTvalue = (P1 + P2 )*Vol/Temp;
            printValue(PVTvalue);
				std::cout << "PV/T = " << PVTvalue << std::endl;        
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printCaloric() {
    std::cout << "calculating caloric: ";
    out << "caloric:" << std::endl << std::endl;
    ThomasFermi* tf = (ThomasFermi*) models[TF];
    ThomasFermiCorrection* tfc = (ThomasFermiCorrection*) models[TFC];
    std::vector<Volume>::iterator vIter;
    std::vector<Temperature>::iterator tIter;
    double P1 = 0;
    double P2 = 0;
	double E1 = 0;
	double E2 = 0;
    double Temp;
    double Vol;
    double caloricValue;
    
    printFirstLine();

    for (tIter = T.begin(); tIter != T.end(); ++tIter) {
        printValue(tIter->getValue());
        for (vIter = V.begin(); vIter != V.end(); ++vIter) {
            if (tf != NULL) { 
                tf->setParameters(*vIter, *tIter);
                P1 = tf->getPressure().getValue();
					 E1 = tf->getEnergy().getValue();
            }
            if (tfc != NULL) {
                tfc->setParameters(*vIter, *tIter);
                P2 = tfc->getPressure().getValue();
					 E2 = tfc->getEnergy().getValue();
            }
            Temp = tIter->transformToUnits(Atomic).transformToScale(lin).getValue();
            Vol = vIter->transformToUnits(Atomic).transformToScale(lin).getValue();
            caloricValue = E1 + E2 - 1.5*(P1 + P2)*Vol;
            printValue(caloricValue);
				std::cout << "caloric value = " << caloricValue << std::endl;
            
        }
        out << std::endl;
    }

    out << std::endl;
    std::cout << "finished" << std::endl;
}

void Interface::printValidityBoundary() {
	ThomasFermi* tf = (ThomasFermi*) models[TF];
	ThomasFermiCorrection* tfc = (ThomasFermiCorrection*) models[TFC];
	std::cout << "calculating validity boundary" << std::endl;
	
	double lgTright = 0.0;
	double lgTleft = -4.0;
	double lgTcurrent;
	double lgNcurrent = 18.5;
	
	double E;
	double dE;

	double M;
	double dM;
	double Ne;

	double eps = 1.0;

	Temperature T;
	Volume V;

	//calculate start point

	V.setValue(pow(10, 25 - lgNcurrent - log10(1.4818)));
	while (fabs(lgTright - lgTleft) > 1e-5) {
		lgTcurrent = 0.5*(lgTright + lgTleft);
		T.setValue(pow(10, lgTcurrent));

		tf->setParameters(V, T);
		tfc->setParameters(V, T);

		E = tf->getPressure().getValue();
		dE = tfc->getPressure().getValue();

		if (fabs(dE/E) < eps) {
			lgTright = lgTcurrent;
		}
		else {
			lgTleft = lgTcurrent;
		}
	}

	lgTcurrent = 0.5*(lgTleft + lgTright);
	
	double distance = 0.1; // set step from point to point on lgN, lgT plane
	double phi1 = 0.01; // min angle to find next point
	double phi2 = M_PI - 0.01; //max angle to find next point
	double phiCurrent; // angle corresponds the exact point of the boundary curve
	double dLgT;
	double dLgN;

	while (lgTcurrent > -4.0) {

		M = tf->getChemicalPotential().getValue();
		dM = tfc->getChemicalPotential().getValue();
		Ne = sqrt(2 * pow(10, 1.5*lgTcurrent)) / M_PI / M_PI*FermiDirac::Half::function((M) / pow(10, lgTcurrent));

		std::cout << "lgN = " << lgNcurrent;
		std::cout << "; lgT = " << lgTcurrent;
		std::cout << "; Ne = " << log10(Ne / 1.4818*1e+25);
		std::cout << "; dE/E = " << fabs(dE / E) << std::endl;

		printValue(lgNcurrent);
		printValue(lgTcurrent + log10(27.229));
		printValue(log10(Ne / 1.4818) + 25);

		out << std::endl;

		while (fabs(phi1 - phi2) > 1e-5) {
			phiCurrent = 0.5*(phi1 + phi2);
			dLgT = distance*cos(phiCurrent);
			dLgN = distance*sin(phiCurrent);
			T.setValue(pow(10, lgTcurrent + dLgT));
			V.setValue(pow(10, 25 - lgNcurrent - dLgN - log10(1.4818)));

			tf->setParameters(V, T);
			tfc->setParameters(V, T);

			E = tf->getPressure().getValue();
			dE = tfc->getPressure().getValue();

			if (fabs(dE/E) < eps) {
				phi1 = phiCurrent;
			}
			else {
				phi2 = phiCurrent;
			}
		}
		phi1 = 0.01;
		phi2 = M_PI + 0.01;
		lgNcurrent += dLgN;
		lgTcurrent += dLgT;
	}
}

void Interface::printThermalValidityBoundary() {
	ThomasFermi* tf = (ThomasFermi*) models[TF];
	ThomasFermiCorrection* tfc = (ThomasFermiCorrection*) models[TFC];
	std::cout << "calculating thermal validity boundary" << std::endl;
	
	double lgTright = 0.0;
	double lgTleft = -4.0;
	double lgTcurrent;
	double lgNcurrent = 18.5;
	
	double E;
	double dE;

	double M;
	double dM;
	double Ne;

	double eps = 1.0;

	Temperature T;
	Volume V;

	//calculate start point
	V.setValue(pow(10, 25 - lgNcurrent - log10(1.4818)));
	while (fabs(lgTright - lgTleft) > 1e-5) {
		lgTcurrent = 0.5*(lgTright + lgTleft);
		T.setValue(pow(10, lgTcurrent));

		tf->setParameters(V, T);
		tfc->setParameters(V, T);

		E = tf->getThermalEnergy().getValue();
		dE = tfc->getThermalEnergy().getValue();

		if (fabs(dE/E) < eps) {
			lgTright = lgTcurrent;
		}
		else {
			lgTleft = lgTcurrent;
		}
	}

	lgTcurrent = 0.5*(lgTleft + lgTright);

	double distance = 0.1; // set step from point to point on lgN, lgT plane
	double phi1 = 0.05; // min angle to find next point
	double phi2 = M_PI + 0.05; //max angle to find next point
	double phiCurrent; // angle corresponds the exact point of the boundary curve
	double dLgT;
	double dLgN;

	while (lgTcurrent > -3.5) {

		M = tf->getChemicalPotential().getValue();
		dM = tfc->getChemicalPotential().getValue();
		Ne = sqrt(2 * pow(10, 1.5*lgTcurrent)) / M_PI / M_PI*FermiDirac::Half::function((M) / pow(10, lgTcurrent));

		std::cout << "lgN = " << lgNcurrent;
		std::cout << "; lgT = " << lgTcurrent + log10(27.229);
		std::cout << "; Ne = " << log10(Ne / 1.4818*1e+25);
		std::cout << "; dE/E = " << fabs(dE / E);
		std::cout << "; P + dP = " << fabs(E + dE) << std::endl;

		printValue(lgNcurrent);
		printValue(lgTcurrent + log10(27.229));
		printValue(log10(Ne / 1.4818) + 25);
		printValue(dE / E);
		out << std::endl;

		while (fabs(phi1 - phi2) > 1e-5) {
			phiCurrent = 0.5*(phi1 + phi2);
			dLgT = distance*cos(phiCurrent);
			dLgN = distance*sin(phiCurrent);
			T.setValue(pow(10, lgTcurrent + dLgT));
			V.setValue(pow(10, 25 - lgNcurrent - dLgN - log10(1.4818)));

			tf->setParameters(V, T);
			tfc->setParameters(V, T);

			E = tf->getThermalEnergy().getValue();
			dE = tfc->getThermalEnergy().getValue();

			if (fabs(dE/E) < eps) {
				phi1 = phiCurrent;
			}
			else {
				phi2 = phiCurrent;
			}
		}
		phi1 = 0.01;//acos(dLgT / distance) - 1.0 / 3.0*acos(dLgT / distance);
		phi2 = M_PI - 0.01;//acos(dLgT / distance) + 1.0 / 3.0*acos(dLgT / distance);
		lgNcurrent += dLgN;
		lgTcurrent += dLgT;

	}
}
