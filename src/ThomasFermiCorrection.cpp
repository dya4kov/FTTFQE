#include "../header/ThomasFermiCorrection.h"

ThomasFermiCorrection::ThomasFermiCorrection(Charge Z) : Model(Z) {
    potential = new TFPotential;
    potentialCorrection = new TFPotentialCorrection;
    coldPotential = new TFPotential;
    coldPotentialCorrection = new TFPotentialCorrection;
    energySolver = new ODESolver<TFCorrection::Energy>;
    entropySolver = new ODESolver<TFCorrection::Entropy>;
}

ThomasFermiCorrection::~ThomasFermiCorrection(void) {
    delete potential;
    delete potentialCorrection;
    delete coldPotential;
    delete coldPotentialCorrection;
    delete energySolver;
    delete entropySolver;
}

std::string ThomasFermiCorrection::getName() {
    return "TFC";
}

Temperature& ThomasFermiCorrection::makeChargeless(Temperature& T) {
    T.setValue(T.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -4.0/3.0));
    return T;
}

Volume& ThomasFermiCorrection::makeChargeless(Volume& V) {
    V.setValue(V.transformToScale(lin).transformToUnits(Atomic).getValue()*Z);
    return V;
}

Pressure& ThomasFermiCorrection::makeChargeless(Pressure& P) {
    P.setValue(P.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -8.0/3.0));
    return P;
}

Energy& ThomasFermiCorrection::makeChargeless(Energy& E) {
    E.setValue(E.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -5.0/3.0));
    return E;
}

Entropy& ThomasFermiCorrection::makeChargeless(Entropy& S) {
    S.setValue(S.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -1.0/3.0));
    return S;
}

ChemicalPotential& ThomasFermiCorrection::makeChargeless(ChemicalPotential& M) {
    M.setValue(M.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, -2.0/3.0));
    return M;
}

Temperature& ThomasFermiCorrection::makeChargeful(Temperature& T) {
    T.setValue(T.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 4.0/3.0));
    return T;
}

Volume& ThomasFermiCorrection::makeChargeful(Volume& V) {
    V.setValue(V.transformToScale(lin).transformToUnits(Atomic).getValue()/Z);
    return V;
}

Pressure& ThomasFermiCorrection::makeChargeful(Pressure& P) {
    P.setValue(P.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 8.0/3.0));
    return P;
}

Energy& ThomasFermiCorrection::makeChargeful(Energy& E) {
    E.setValue(E.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 5.0/3.0));
    return E;
}

Entropy& ThomasFermiCorrection::makeChargeful(Entropy& S) {
    S.setValue(S.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 1.0/3.0));
    return S;
}

ChemicalPotential& ThomasFermiCorrection::makeChargeful(ChemicalPotential& M) {
    M.setValue(M.transformToScale(lin).transformToUnits(Atomic).getValue()*pow(Z, 2.0/3.0));
    return M;
}

void ThomasFermiCorrection::getReady() {
    potential->setParameters(makeChargeless(V).getValue(), makeChargeless(T).getValue());
    potential->setPrecision(eps);
    potentialCorrection->setParameters(V.getValue(), T.getValue());
    potentialCorrection->setPrecision(eps);
    makeChargeful(V);
    makeChargeful(T);
    readyToCalculate = true;
}

void ThomasFermiCorrection::calculatePressure() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    makeChargeless(V);
    double value = pow(T.getValue(), 2.0)/3.0/M_PI/M_PI/M_PI
                 * (FermiDirac::Half::function
                        (
                                potential->valueAt_1()
                        )
                    * potentialCorrection->valueAt_1() 
                    + Y::function(potential->valueAt_1())
                 );
    P.setValue(value, Atomic);
    calculatedPressure = true;
    makeChargeful(P);
    makeChargeful(T);
    makeChargeful(V);
}

void ThomasFermiCorrection::calculateThermalP() {
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

void ThomasFermiCorrection::calculateEnergy() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    makeChargeless(V);
    double value;
    double* initials = new double[TFCorrection::Energy::dimension];

    initials[0] = potential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = potentialCorrection->valueAt_1();
    initials[3] = initials[2];
    initials[4] = sqrt(T.getValue()*0.5)/3.0/M_PI*potentialCorrection->derivativeAt_0() + 0.269900170;

    energySolver->setLimits(1.0, 0.0);
    energySolver->setInitials(initials);
    energySolver->setAccuracy(0.0, eps/10);
	
    TFCorrection::Energy::parameter[0] = pow(2.0, 7.0/6.0)
                                       * pow(3.0, 2.0/3.0)
                                       * pow(M_PI, -5.0/3.0)
                                       * sqrt(T.getValue())*pow(V.getValue(), 2.0/3.0);
    TFCorrection::Energy::parameter[1] = V.getValue()*pow(T.getValue(), 2.0)/M_PI/M_PI/M_PI;

    value = energySolver->driverApply().y[4];
    E.setValue(value);

    calculatedEnergy = true;
    energySolver->reset();
    makeChargeful(E);
    makeChargeful(T);
    makeChargeful(V);
}

void ThomasFermiCorrection::calculateThermalE() {
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

void ThomasFermiCorrection::calculateEntropy() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    makeChargeless(V);
    double value;
    double* initials = new double[TFCorrection::Entropy::dimension];

    initials[0] = potential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = potentialCorrection->valueAt_1();
    initials[3] = initials[2];
    initials[4] = 1.0/(sqrt(2*T.getValue())*3.0*M_PI)*potentialCorrection->derivativeAt_0();

    entropySolver->setLimits(1.0, 0.0);
    entropySolver->setInitials(initials);
    entropySolver->setAccuracy(0.0, eps/10);
	
    TFCorrection::Entropy::parameter[0] = pow(2.0, 7.0/6.0)
                                        * pow(3.0, 2.0/3.0)
                                        * pow(M_PI, -5.0/3.0)
                                        * sqrt(T.getValue())*pow(V.getValue(), 2.0/3.0);
    TFCorrection::Entropy::parameter[1] = V.getValue()*T.getValue()/M_PI/M_PI/M_PI;

    value = entropySolver->driverApply().y[4];
    S.setValue(value);

    calculatedEntropy = true;
    entropySolver->reset();
    makeChargeful(S);
    makeChargeful(T);
    makeChargeful(V);
}

void ThomasFermiCorrection::calculateThermalS() {
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


void ThomasFermiCorrection::calculateChemPotential() {
    if (!readyToCalculate) getReady();
    makeChargeless(T);
    double value = sqrt(T.getValue()*0.5)/3.0/M_PI
                 * (0.5*FermiDirac::MHalf::function(potential->valueAt_1()) 
                 + potentialCorrection->valueAt_1()
                   );
    M.setValue(value, Atomic);
    calculatedChemPot = true;
    makeChargeful(M);
    makeChargeful(T);
}

void ThomasFermiCorrection::calculateThermalCP() {
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

double ThomasFermiCorrection::coldP() {
    coldPotential->setParameters(V.getValue(), coldT);
    coldPotential->setPrecision(eps);
    coldPotentialCorrection->setParameters(V.getValue(), coldT);
    coldPotentialCorrection->setPrecision(eps);
    
    double value = pow(coldT, 2.0)/3.0/M_PI/M_PI/M_PI
                 * (FermiDirac::Half::function
                        (
                                coldPotential->valueAt_1()
                        )
                    * coldPotentialCorrection->valueAt_1() 
                    + Y::function(coldPotential->valueAt_1())
                 );
    return value;
}

double ThomasFermiCorrection::coldE() {
    double value;
    double* initials = new double[TFCorrection::Energy::dimension];

    potential->setParameters(V.getValue(), coldT);
    potential->setPrecision(eps);
    potentialCorrection->setParameters(V.getValue(), coldT);
    potentialCorrection->setPrecision(eps);
    
    initials[0] = potential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = potentialCorrection->valueAt_1();
    initials[3] = initials[2];
    initials[4] = sqrt(coldT*0.5)/3.0/M_PI*potentialCorrection->derivativeAt_0() + 0.269900170;

    energySolver->setLimits(1.0, 0.0);
    energySolver->setInitials(initials);
    energySolver->setAccuracy(0.0, eps/10);
	
    TFCorrection::Energy::parameter[0] = pow(2.0, 7.0/6.0)
                                       * pow(3.0, 2.0/3.0)
                                       * pow(M_PI, -5.0/3.0)
                                       * sqrt(coldT)*pow(V.getValue(), 2.0/3.0);
    TFCorrection::Energy::parameter[1] = V.getValue()*pow(coldT, 2.0)/M_PI/M_PI/M_PI;

    value = energySolver->driverApply().y[4];
    energySolver->reset();
    return value;
}

double ThomasFermiCorrection::coldS() {
    double value;
    double* initials = new double[TFCorrection::Entropy::dimension];

    potential->setParameters(V.getValue(), coldT);
    potential->setPrecision(eps);
    potentialCorrection->setParameters(V.getValue(), coldT);
    potentialCorrection->setPrecision(eps);
    
    initials[0] = potential->valueAt_1();
    initials[1] = initials[0];
    initials[2] = potentialCorrection->valueAt_1();
    initials[3] = initials[2];
    initials[4] = 1.0/(sqrt(2*coldT)*3.0*M_PI)*potentialCorrection->derivativeAt_0();

    entropySolver->setLimits(1.0, 0.0);
    entropySolver->setInitials(initials);
    entropySolver->setAccuracy(0.0, eps/10);
	
    TFCorrection::Entropy::parameter[0] = pow(2.0, 7.0/6.0)
                                        * pow(3.0, 2.0/3.0)
                                        * pow(M_PI, -5.0/3.0)
                                        * sqrt(coldT)*pow(V.getValue(), 2.0/3.0);
    TFCorrection::Entropy::parameter[1] = V.getValue()*coldT/M_PI/M_PI/M_PI;

    value = entropySolver->driverApply().y[4];
    entropySolver->reset();
    return value;
}

double ThomasFermiCorrection::coldM() {
    coldPotential->setParameters(V.getValue(), coldT);
    coldPotential->setPrecision(eps);
    coldPotentialCorrection->setParameters(V.getValue(), coldT);
    coldPotentialCorrection->setPrecision(eps);
    
    double value = sqrt(coldT*0.5)/3.0/M_PI
                 * (0.5*FermiDirac::MHalf::function(coldPotential->valueAt_1()) 
                 + coldPotentialCorrection->valueAt_1()
                   );
    return value;
}