#include "../header/ThermodynamicFunctions.h"

ThermodynamicFunction::ThermodynamicFunction() {
    defaultValue = 0;
    currentValue = 0;
    unit = Atomic;
    scale = lin;
}

ThermodynamicFunction::~ThermodynamicFunction() {

}

Scaling ThermodynamicFunction::getCurrentScale() {
    return scale;
}

Unit ThermodynamicFunction::getCurrentUnit() {
    return unit;
}

double ThermodynamicFunction::getValue() {
    return currentValue;
}

Temperature::Temperature() : ThermodynamicFunction() {}
Volume::Volume() : ThermodynamicFunction() {}
Pressure::Pressure() : ThermodynamicFunction() {}
Energy::Energy() : ThermodynamicFunction() {}
Entropy::Entropy() : ThermodynamicFunction() {}
ChemicalPotential::ChemicalPotential() : ThermodynamicFunction() {}
Density::Density() : ThermodynamicFunction() {}
Concentration::Concentration() : ThermodynamicFunction() {}

Temperature::Temperature(const Temperature& T) {
    defaultValue = T.defaultValue;
    currentValue = T.currentValue;
    unit = T.unit;
    scale = T.scale;
}

Volume::Volume(const Volume& V) {
    defaultValue = V.defaultValue;
    currentValue = V.currentValue;
    unit = V.unit;
    scale = V.scale;
}
Energy::Energy(const Energy& E) {
    defaultValue = E.defaultValue;
    currentValue = E.currentValue;
    unit = E.unit;
    scale = E.scale;
}
Pressure::Pressure(const Pressure& P) {
    defaultValue = P.defaultValue;
    currentValue = P.currentValue;
    unit = P.unit;
    scale = P.scale;
}
Entropy::Entropy(const Entropy& S) {
    defaultValue = S.defaultValue;
    currentValue = S.currentValue;
    unit = S.unit;
    scale = S.scale;
}
ChemicalPotential::ChemicalPotential(const ChemicalPotential& M) {
    defaultValue = M.defaultValue;
    currentValue = M.currentValue;
    unit = M.unit;
    scale = M.scale;
}

Density::Density(const Density& D) {
    defaultValue = D.defaultValue;
    currentValue = D.currentValue;
    unit = D.unit;
    scale = D.scale;
}

Concentration::Concentration(const Concentration& C) {
    defaultValue = C.defaultValue;
    currentValue = C.currentValue;
    unit = C.unit;
    scale = C.scale;
}

Temperature& Temperature::operator=(const Temperature& T) {
    defaultValue = T.defaultValue;
    currentValue = T.currentValue;
    unit = T.unit;
    scale = T.scale;
    return *this;
}

Volume& Volume::operator=(const Volume& V) {
    defaultValue = V.defaultValue;
    currentValue = V.currentValue;
    unit = V.unit;
    scale = V.scale;
    return *this;
}
Energy& Energy::operator=(const Energy& E) {
    defaultValue = E.defaultValue;
    currentValue = E.currentValue;
    unit = E.unit;
    scale = E.scale;
    return *this;
}
Pressure& Pressure::operator=(const Pressure& P) {
    defaultValue = P.defaultValue;
    currentValue = P.currentValue;
    unit = P.unit;
    scale = P.scale;
    return *this;
}
Entropy& Entropy::operator=(const Entropy& S) {
    defaultValue = S.defaultValue;
    currentValue = S.currentValue;
    unit = S.unit;
    scale = S.scale;
    return *this;
}
ChemicalPotential& ChemicalPotential::operator=(const ChemicalPotential& M) {
    defaultValue = M.defaultValue;
    currentValue = M.currentValue;
    unit = M.unit;
    scale = M.scale;
    return *this;
}

Density& Density::operator=(const Density& D) {
    defaultValue = D.defaultValue;
    currentValue = D.currentValue;
    unit = D.unit;
    scale = D.scale;
    return *this;
}

Concentration& Concentration::operator=(const Concentration& C) {
    defaultValue = C.defaultValue;
    currentValue = C.currentValue;
    unit = C.unit;
    scale = C.scale;
    return *this;
}

void Temperature::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    switch (unit) {
        case Atomic : 
            defaultValue = newValue;
            break;
        case eV : 
            defaultValue = newValue/27.229;
            break;
        case Kelvin : 
            defaultValue = newValue*8.6173324e-5/27.229;
            break;
    }

    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

void Volume::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    switch (unit) {
        case Atomic : 
            defaultValue = newValue;
            break;
        case cmc : 
            defaultValue = newValue/1.4818e-25;
            break;
        case mc : 
            defaultValue = newValue/1.4818e-31;
            break;
        }

    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

void Pressure::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    switch (unit) {
        case Atomic : 
            defaultValue = newValue;
            break;
        case MBar : 
            defaultValue = newValue/294.18;
            break;
        case GPa : 
            defaultValue = newValue/29418.0;
            break;
		  case Pa:
			   defaultValue = newValue/29418e+9;
				break;
    }

    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

void Energy::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    switch (unit) {
        case Atomic : 
            defaultValue = newValue;
            break;
        case eV : 
            defaultValue = newValue/27.229;
            break;
    }

    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

void Entropy::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    defaultValue = newValue;
    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

void ChemicalPotential::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    switch (unit) {
        case Atomic : 
            defaultValue = newValue;
            break;
        case eV : 
            defaultValue = newValue/27.229;
            break;
    }

    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

void Density::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    switch (unit) {
        case gOverCmc : 
            defaultValue = newValue;
            break;
        case kgOverMc : 
            defaultValue = newValue*1000;
            break;
        }

    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

void Concentration::setValue(double value, Unit unit, Scaling scale) {
    double newValue;
    if (scale == lg) newValue = pow(10, value);
    else newValue = value;

    switch (unit) {
        case Atomic : 
            defaultValue = newValue;
            break;
        case oneOverCmc : 
            defaultValue = newValue*1.4818e-25;
            break;
        case oneOverMc : 
            defaultValue = newValue*1.4818e-31;
            break;
        }

    currentValue = value;
    this->unit = unit;
    this->scale = scale;
}

Temperature& Temperature::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case Atomic : 
            currentValue = defaultValue;
            break;
        case eV : 
            currentValue = defaultValue*27.229;
            break;
        case Kelvin : 
            currentValue = defaultValue*27.229/8.6173324e-5;
            break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

Volume& Volume::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case Atomic : 
            currentValue = defaultValue;
            break;
        case cmc : 
            currentValue = defaultValue*1.4818e-25;
            break;
        case mc :
            currentValue = defaultValue*1.4818e-31;
            break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

Pressure& Pressure::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case Atomic : 
            currentValue = defaultValue;
            break;
        case GPa : 
            currentValue = defaultValue*29418;
            break;
        case MBar :
            currentValue = defaultValue*294.18;
            break;
		  case Pa:
			   currentValue = defaultValue*29418e+9;
				break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

Energy& Energy::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case Atomic : 
            currentValue = defaultValue;
            break;
        case eV : 
            currentValue = defaultValue*27.229;
            break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

Entropy& Entropy::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case Atomic : 
            currentValue = defaultValue;
            break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

ChemicalPotential& ChemicalPotential::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case Atomic : 
            currentValue = defaultValue;
            break;
        case eV : 
            currentValue = defaultValue*27.229;
            break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

Density& Density::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case gOverCmc :
            currentValue = defaultValue;
            break;
        case kgOverMc :
            currentValue = defaultValue*1000;
            break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

Concentration& Concentration::transformToUnits(Unit newUnit) {
    switch (newUnit) {
        case Atomic : 
            currentValue = defaultValue;
            break;
        case oneOverCmc : 
            currentValue = defaultValue/1.4818e-25;
            break;
        case oneOverMc : 
            currentValue = defaultValue/1.4818e-31;
            break;
    }
    if (scale == lg) currentValue = log10(currentValue);
    return *this;
}

Temperature& Temperature::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg) {
            currentValue = log10(fabs(currentValue));
        }
        else currentValue = pow(10, currentValue);
    scale = newScale;
    return *this;
}

Volume& Volume::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg)
            currentValue = log10(fabs(currentValue));
        else currentValue = pow(10, currentValue);
    scale = newScale;
    return *this;
}

Pressure& Pressure::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg)
            currentValue = log10(fabs(currentValue));
        else currentValue = pow(10, currentValue);
    scale = newScale;
    return *this;
}

Energy& Energy::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg)
            currentValue = log10(fabs(currentValue));
        else currentValue = pow(10, currentValue);
        scale = newScale;
    return *this;
}

Entropy& Entropy::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg)
            currentValue = log10(fabs(currentValue));
        else currentValue = pow(10, currentValue);
    scale = newScale;
    return *this;
}

ChemicalPotential& ChemicalPotential::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg)
            currentValue = log10(fabs(currentValue));
        else currentValue = pow(10, currentValue);
    scale = newScale;
    return *this;
}

Density& Density::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg)
            currentValue = log10(fabs(currentValue));
        else currentValue = pow(10, currentValue);
    scale = newScale;
    return *this;
}

Concentration& Concentration::transformToScale(Scaling newScale) {
    if (newScale != scale)
        if (newScale == lg)
            currentValue = log10(fabs(currentValue));
        else currentValue = pow(10, currentValue);
    scale = newScale;
    return *this;
}