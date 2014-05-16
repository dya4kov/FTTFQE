#include "../header/ThomasFermiPotentialCorrection.h"


TFPotentialCorrection::TFPotentialCorrection(void) :
                        vTableSize(191), 
                        tTableSize(171), 
                        lgV0(-9.0), 
                        lgT0(-7.0),
                        lgVStep(0.1),
                        lgTStep(0.1) 
{
    solver = new ODESolver<TFCorrection::Potential>;
    potential = new TFPotential;
    setPsiTable();
    V = 0;
    T = 0;
}


TFPotentialCorrection::~TFPotentialCorrection(void) {
    delete solver;
    delete potential;
    for (int i = 0; i < vTableSize; ++i) {
        delete[] psiTable[i];
    }
    delete[] psiTable;
}

void TFPotentialCorrection::setPsiTable() {
    std::ifstream data("potentialData/TFPotentialCorrection.dat", std::ios::in);
    double currentValue;
    psiTable = new double*[vTableSize];
    for (int i = 0; i < vTableSize; ++i)
        psiTable[i] = new double[tTableSize];

    for (int t = 0; t < tTableSize; ++t) {
        for (int v = 0; v < vTableSize; ++v) {
            data >> currentValue;
            psiTable[v][t] = currentValue;
        }
    }
    data.close();
}

void TFPotentialCorrection::setInitialShotParameters(double& p1, double& p2) {
    double psi1_min = 0;
    double psi1_max = 0;
    int v;
    int t;
    bool v_is_calculated = false;
    bool t_is_calculated = false;

    v = 0;
    t = 0;

    v = (int) ceil((log10(V) - lgV0)*10);
    if (fabs(log10(V) - lgV0 - v*lgVStep) < 1e-10)
    {
        v_is_calculated = true;
    }
    else if (fabs(log10(V) - lgV0 - (v - 1)*lgVStep) < 1e-10)
    {
        v_is_calculated = true;
        v--;
    }

    t = (int) ceil((log10(T) - lgT0)*10);
    if (fabs(log10(T) - lgT0 - t*lgTStep) < 1e-10)
    {
        t_is_calculated = true;
    }
    else if (fabs(log10(T) - lgT0 - (t - 1)*lgTStep) < 1e-10)
    {
        t_is_calculated = true;
        t--;
    }

    if (v_is_calculated && t_is_calculated)
    {
        p1 = p2 = psiTable[v][t];
    }
    else 
    {
        if (v == 0 || v_is_calculated) { 
            psi1_max = psiTable[v][t];
            psi1_min = psiTable[v][t - 1];
        }
        else if (t == 0 || t_is_calculated) {
            psi1_max = psiTable[v][t];
            psi1_min = psiTable[v - 1][t];
        }
        else {
            psi1_max = psiTable[v][t];
            psi1_min = psiTable[v - 1][t - 1];
        }
        p1 = psi1_min;
        p2 = psi1_max;
    }
}

void TFPotentialCorrection::calculate() {
    Solution currentSolution(TFCorrection::Potential::dimension);
    double shotParameter_1;
    double shotParameter_2;
    double initials[TFCorrection::Potential::dimension] = { 0 };

    setInitialShotParameters(shotParameter_1, shotParameter_2);

    initials[0] = potential->valueAt_1();
    initials[1] = initials[0];

    if (fabs(shotParameter_1 - shotParameter_2)/fabs(shotParameter_1) < eps/10) {
        shotParameter_1 = shotParameter_1 - eps*fabs(shotParameter_1);
        shotParameter_2 = shotParameter_1 + eps*fabs(shotParameter_2);
    }

    solver->setLimits(1.0, 0);
    solver->setAccuracy(0.0, eps/10);

    while (fabs(shotParameter_2 - shotParameter_1)/fabs(shotParameter_1) > eps) {
        initials[2] = (shotParameter_1 + shotParameter_2)/2;
        initials[3] = initials[2];
        solver->setInitials(initials);
        currentSolution = solver->driverApply();
        if (currentSolution.y[2] > 0)
            shotParameter_2 -= (shotParameter_2 - shotParameter_1)/2;
        else shotParameter_1 += (shotParameter_2 - shotParameter_1)/2;
        solver->reset();
    }
    psi_1 = (shotParameter_1 + shotParameter_2)/2;
    calculated = true;
}

double TFPotentialCorrection::valueAt_1() {
    if (!calculated) calculate();
    return psi_1;
}

double TFPotentialCorrection::valueAt_x(double x) {
    if (!calculated) calculate();
    double initials[2] = { psi_1, psi_1 };
    double result;
    solver->setInitials(initials);
    solver->setLimits(1.0, x);
    result = solver->driverApply().y[0];
    solver->reset();
    return result;
}

double TFPotentialCorrection::derivativeAt_0() {
    if (!calculated) calculate();
    double initials[4] = {  potential->valueAt_1(), 
                            potential->valueAt_1(),
                            psi_1,
                            psi_1
                         };
    double result;
    solver->setInitials(initials);
    solver->setLimits(1.0, 0);
    result = solver->driverApply().y[3];
    solver->reset();
    return result;
}

void TFPotentialCorrection::setParameters(double V, double T) {
    if (fabs(log10(V) - log10(this->V)) > 1e-10
            || fabs(log10(T) - log10(this->T)) > 1e-10) 
    { 
        calculated = false;
        this->V = V;
        this->T = T;
        potential->setParameters(V, T);
        *TFCorrection::Potential::parameter = pow(2.0, 7.0/6.0)
                                            * pow(3.0, 2.0/3.0)
                                            * pow(M_PI, -5.0/3.0)
                                            * sqrt(T)*pow(V, 2.0/3.0);
    }
}

void TFPotentialCorrection::setPrecision(double eps) {
    potential->setPrecision(eps);
    this->eps = eps;
}
