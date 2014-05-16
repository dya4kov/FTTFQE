#include "../header/Model.h"
#include <iostream>

Model::Model(Charge Z) : coldT(1e-7) {
    this->Z = Z;
    calculatedPressure = false;
    calculatedEnergy = false;
    calculatedEntropy = false;
    calculatedChemPot = false;
    calculatedThermalP = false;
    calculatedThermalE = false;
    calculatedThermalS = false;
    calculatedThermalCP = false;
    parameterSet = false;
    precisionSet = false;
    readyToCalculate = false;
}

Model::~Model() {

}

void Model::calculatePressure() {}
void Model::calculateEnergy() {}
void Model::calculateEntropy() {}
void Model::calculateChemPotential() {}

void Model::calculateThermalP() {}
void Model::calculateThermalE() {}
void Model::calculateThermalS() {}
void Model::calculateThermalCP() {}

void Model::getReady() {}

void Model::setParameters(Volume V, Temperature T) {
    this->V = V;
    this->T = T;
    parameterSet = true;
    readyToCalculate = false;
    calculatedChemPot = false;
    calculatedEnergy = false;
    calculatedEntropy = false;
    calculatedPressure = false;
    calculatedThermalP = false;
    calculatedThermalE = false;
    calculatedThermalS = false;
    calculatedThermalCP = false;

}

void Model::setPrecision(const double eps) {
    this->eps = eps;
    precisionSet = true;
    readyToCalculate = false;
    calculatedChemPot = false;
    calculatedEnergy = false;
    calculatedEntropy = false;
    calculatedPressure = false;
    calculatedThermalP = false;
    calculatedThermalE = false;
    calculatedThermalS = false;
    calculatedThermalCP = false;
}

void Model::setAtomicNumber(Charge Z) {
    this->Z = Z;
}

std::string Model::getName() { return " "; }

Temperature Model::getTemperature() {
    return T;
}

Volume Model::getVolume() {
    return V;
}

Pressure Model::getPressure() {
    if (!calculatedPressure) calculatePressure();
    return P;
}

Energy Model::getEnergy() {
    if (!calculatedEnergy) calculateEnergy();
    return E;
}

Entropy Model::getEntropy() {
    if (!calculatedEntropy) calculateEntropy();
    return S;
}

ChemicalPotential Model::getChemicalPotential() {
    if (!calculatedChemPot) calculateChemPotential();
    return M;
}

Pressure Model::getThermalPressure() {
    if (!calculatedThermalP) calculateThermalP();
    return Pthermal;
}

Energy Model::getThermalEnergy() {
    if (!calculatedThermalE) calculateThermalE();
    return Ethermal;
}

Entropy Model::getThermalEntropy() {
    if (!calculatedThermalS) calculateThermalS();
    return Sthermal;
}

ChemicalPotential Model::getThermalChemicalPotential() {
    if (!calculatedThermalCP) calculateThermalCP();
    return Mthermal;
}
