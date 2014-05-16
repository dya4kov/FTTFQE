#include "../header/ThomasFermi.h"


ThomasFermi::ThomasFermi(Charge Z) : Model(Z) {
    potential = new TFPotential;
    coldPotential = new TFPotential;
    energySolver = new ODESolver<TF::Energy>;
    entropySolver = new ODESolver<TF::Entropy>;
}

ThomasFermi::~ThomasFermi(void) {
    delete potential;
    delete coldPotential;
    delete energySolver;
    delete entropySolver;
}

std::string ThomasFermi::getName() {
    return "TF";
}

Temperature& ThomasFermi::makeChargeless(Temperature& T) {
    T.setValue(T.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -4.0/3.0));
    return T;
}

Volume& ThomasFermi::makeChargeless(Volume& V) {
    V.setValue(V.transformToScale(lin).transformToUnits(Atomic).getValue()*Z);
    return V;
}

Pressure& ThomasFermi::makeChargeless(Pressure& P) {
    P.setValue(P.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -10.0/3.0));
    return P;
}

Energy& ThomasFermi::makeChargeless(Energy& E) {
    E.setValue(E.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -7.0/3.0));
    return E;
}

Entropy& ThomasFermi::makeChargeless(Entropy& S) {
    S.setValue(S.transformToScale(lin).transformToUnits(Atomic).getValue()/Z);
    return S;
}

ChemicalPotential& ThomasFermi::makeChargeless(ChemicalPotential& M) {
    M.setValue(M.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -4.0/3.0));
    return M;
}

Temperature& ThomasFermi::makeChargeful(Temperature& T) {
    T.setValue(T.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 4.0/3.0));
    return T;
}

Volume& ThomasFermi::makeChargeful(Volume& V) {
    V.setValue(V.transformToScale(lin).transformToUnits(Atomic).getValue()/Z);
    return V;
}

Pressure& ThomasFermi::makeChargeful(Pressure& P) {
    P.setValue(P.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 10.0/3.0));
    return P;
}

Energy& ThomasFermi::makeChargeful(Energy& E) {
    E.setValue(E.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 7.0/3.0));
    return E;
}

Entropy& ThomasFermi::makeChargeful(Entropy& S) {
    S.setValue(S.transformToScale(lin).transformToUnits(Atomic).getValue()*Z);
    return S;
}

ChemicalPotential& ThomasFermi::makeChargeful(ChemicalPotential& M) {
    M.setValue(M.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 4.0/3.0));
    return M;
}

void ThomasFermi::getReady() {
    potential->setParameters(makeChargeless(V).getValue(), makeChargeless(T).getValue());
    potential->setPrecision(eps);
    makeChargeful(V);
    makeChargeful(T);
    readyToCalculate = true;
}

void ThomasFermi::calculatePressure() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    makeChargeless(V);
    double value = 
            pow(2*T.transformToScale(lin).getValue(), 5.0/2.0)/6.0/M_PI/M_PI
            *FermiDirac::ThreeHalf::function(potential->valueAt_1());
    P.setValue(value, Atomic);
    calculatedPressure = true;
    makeChargeful(P);
    makeChargeful(T);
    makeChargeful(V);
}

void ThomasFermi::calculateThermalP() {
    if (!calculatedPressure) calculatePressure();
    makeChargeless(T);
    makeChargeless(V);
    makeChargeless(P);
    double value = P.getValue();
    Pthermal.setValue(value - coldP(), Atomic);
    calculatedThermalP= true;
    makeChargeful(P);
    makeChargeful(T);
    makeChargeful(V);
    makeChargeful(Pthermal);
}

void ThomasFermi::calculateEnergy() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    makeChargeless(V);
    double value;
    double* initials = new double[TF::Energy::dimension];

    initials[0] = potential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = 2.0*sqrt(2.0)*V.getValue()*pow(T.getValue(), 5.0/2.0)/M_PI/M_PI
                * FermiDirac::ThreeHalf::function(potential->valueAt_1())
                + 0.76874512422;

    energySolver->setLimits(1.0, 0.0);
    energySolver->setInitials(initials);
    energySolver->setAccuracy(0.0, eps/10);

    TF::Energy::parameter[0] = 1.0;
    TF::Energy::parameter[1] = pow(2.0, 7.0/6.0)
                             * pow(3.0, 2.0/3.0)
                             * pow(M_PI, -5.0/3.0)
                             * sqrt(T.getValue())*pow(V.getValue(), 2.0/3.0);
    TF::Energy::parameter[2] = 3.0*sqrt(2.0)*V.getValue()*pow(T.getValue(), 5.0/2.0)/M_PI/M_PI;

    value = energySolver->driverApply().y[2];

    E.setValue(value);

    calculatedEnergy = true;
    energySolver->reset();
    makeChargeful(E);
    makeChargeful(T);
    makeChargeful(V);
}

void ThomasFermi::calculateThermalE() {
    if (!calculatedEnergy) calculateEnergy();
    makeChargeless(T);
    makeChargeless(V);
    makeChargeless(E);
    double value = E.getValue();
    Ethermal.setValue(value - coldE());
    calculatedThermalE = true;
    makeChargeful(T);
    makeChargeful(V);
    makeChargeful(E);
    makeChargeful(Ethermal);
}

void ThomasFermi::calculateEntropy() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    makeChargeless(V);
    double value;
    double* initials = new double[TF::Entropy::dimension];

    initials[0] = potential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = 4.0*sqrt(2.0*T.getValue())*V.getValue()*T.getValue()/M_PI/M_PI
                * FermiDirac::ThreeHalf::function(potential->valueAt_1())
                - potential->derivativeAt_0();

    entropySolver->setLimits(1.0, 0.0);
    entropySolver->setInitials(initials);
    entropySolver->setAccuracy(0.0, eps/10);

    TF::Entropy::parameter[0] = 1.0;
    TF::Entropy::parameter[1] = pow(2.0, 7.0/6.0)
                              * pow(3.0, 2.0/3.0)
                              * pow(M_PI, -5.0/3.0)
                              * sqrt(T.getValue())*pow(V.getValue(), 2.0/3.0);
    TF::Entropy::parameter[2] = 7.0*sqrt(2.0*T.getValue())*V.getValue()*T.getValue()/M_PI/M_PI;

    value = entropySolver->driverApply().y[2];

    S.setValue(value);

    calculatedEntropy = true;
    entropySolver->reset();
    makeChargeful(S);
    makeChargeful(T);
    makeChargeful(V);
}

void ThomasFermi::calculateThermalS() {
    if (!calculatedEntropy) calculateEntropy();
    makeChargeless(T);
    makeChargeless(V);
    makeChargeless(S);
    double value = S.getValue();
    Sthermal.setValue(value - coldS());
    calculatedThermalS = true;
    makeChargeful(T);
    makeChargeful(V);
    makeChargeful(S);
    makeChargeful(Sthermal);
}

void ThomasFermi::calculateChemPotential() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    double value = T.getValue()*potential->valueAt_1();
    M.setValue(value, Atomic);
    calculatedChemPot = true;
    makeChargeful(M);
    makeChargeful(T);
}

void ThomasFermi::calculateThermalCP() {
    if (!calculatedChemPot) calculateChemPotential();
    makeChargeless(T);
    makeChargeless(V);
    makeChargeless(M);
    double value = M.getValue();
    Mthermal.setValue(value - coldM());
    calculatedThermalCP = true;
    makeChargeful(T);
    makeChargeful(V);
    makeChargeful(M);
    makeChargeful(Mthermal);
}

double ThomasFermi::coldP() {
    coldPotential->setParameters(V.getValue(), coldT);
    coldPotential->setPrecision(eps);
    double value = 
        pow(2*coldT, 5.0/2.0)/6.0/M_PI/M_PI
        *FermiDirac::ThreeHalf::function(coldPotential->valueAt_1());
    return value;
}

double ThomasFermi::coldE() {
    double value;
    double* initials = new double[TF::Energy::dimension];

    coldPotential->setParameters(V.getValue(), coldT);
    coldPotential->setPrecision(eps);
    
    initials[0] = coldPotential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = 2.0*sqrt(2.0)*V.getValue()*pow(coldT, 5.0/2.0)/M_PI/M_PI
                * FermiDirac::ThreeHalf::function(coldPotential->valueAt_1())
                + 0.76874512422;

    energySolver->setLimits(1.0, 0.0);
    energySolver->setInitials(initials);
    energySolver->setAccuracy(0.0, eps/10);

    TF::Energy::parameter[0] = 1.0;
    TF::Energy::parameter[1] = pow(2.0, 7.0/6.0)
                             * pow(3.0, 2.0/3.0)
                             * pow(M_PI, -5.0/3.0)
                             * sqrt(coldT)*pow(V.getValue(), 2.0/3.0);
    TF::Energy::parameter[2] = 3.0*sqrt(2.0)*V.getValue()*pow(coldT, 5.0/2.0)/M_PI/M_PI;

    value = energySolver->driverApply().y[2];

    energySolver->reset();
    return value;
}

double ThomasFermi::coldS() {
    double value;
    double* initials = new double[TF::Entropy::dimension];
    
    coldPotential->setParameters(V.getValue(), coldT);
    coldPotential->setPrecision(eps);

    initials[0] = coldPotential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = 4.0*sqrt(2.0*coldT)*V.getValue()*coldT/M_PI/M_PI
                * FermiDirac::ThreeHalf::function(coldPotential->valueAt_1())
                - coldPotential->derivativeAt_0();

    entropySolver->setLimits(1.0, 0.0);
    entropySolver->setInitials(initials);
    entropySolver->setAccuracy(0.0, eps/10);

    TF::Entropy::parameter[0] = 1.0;
    TF::Entropy::parameter[1] = pow(2.0, 7.0/6.0)
                              * pow(3.0, 2.0/3.0)
                              * pow(M_PI, -5.0/3.0)
                              * sqrt(coldT)*pow(V.getValue(), 2.0/3.0);
    TF::Entropy::parameter[2] = 7.0*sqrt(2.0*coldT)*V.getValue()*coldT/M_PI/M_PI;

    value = entropySolver->driverApply().y[2];

    entropySolver->reset();
    return value;         
}

double ThomasFermi::coldM() {
    coldPotential->setParameters(V.getValue(), coldT);
    coldPotential->setPrecision(eps);
    
    return coldT*coldPotential->valueAt_1();
}