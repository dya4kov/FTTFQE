#include "../header/ThomasFermiPotential.h"

TFPotential::TFPotential(void) : 
                vTableSize(200), 
                tTableSize(175), 
                lgV0(-9.9), 
                lgT0(-7.4),
                lgVStep(0.1),
                lgTStep(0.1) 
{
    solver = new ODESolver<TF::Potential>;
    setPhiTable();
    V = 0;
    T = 0;
}

TFPotential::~TFPotential(void) {

    delete solver;
    for (int i = 0; i < vTableSize; ++i) {
            delete[] phiTable[i];
    }
    delete[] phiTable;
}

void TFPotential::setPhiTable() {
    std::ifstream data("potentialData/TFPotential.dat", std::ios::in);
    double currentValue;
    phiTable = new double*[vTableSize];
    for (int i = 0; i < vTableSize; ++i)
        phiTable[i] = new double[tTableSize];

    for (int t = 0; t < tTableSize; ++t) {
        for (int v = 0; v < vTableSize; ++v) {
            data >> currentValue;
            phiTable[v][t] = currentValue;
        }
    }
    data.close();
}

void TFPotential::setInitialShotParameters(double& p1, double& p2) {
    double phi1_max;
    double phi1_min;
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
        p1 = p2 = phiTable[v][t];
    }
    else 
    {
        if (v == 0 || v_is_calculated) { 
            phi1_min = phiTable[v][t];
            phi1_max = phiTable[v][t - 1];
        }
        else if (t == 0 || t_is_calculated) {
            phi1_min = phiTable[v][t];
            phi1_max = phiTable[v - 1][t];
        }
        else {
            phi1_min = phiTable[v][t];
            phi1_max = phiTable[v - 1][t - 1];
        }
        p1 = phi1_min;
        p2 = phi1_max;
    }
}

void TFPotential::calculate() {
    Solution currentSolution(TF::Potential::dimension);
    double shotParameter_1;
    double shotParameter_2;
    double initials[2] = { 0 };

    setInitialShotParameters(shotParameter_1, shotParameter_2);

    if (fabs(shotParameter_1 - shotParameter_2)/fabs(shotParameter_1) < eps/10) {
        shotParameter_1 = shotParameter_1 - eps*fabs(shotParameter_1);
        shotParameter_2 = shotParameter_1 + eps*fabs(shotParameter_2);
    }
    solver->setLimits(1.0, 0);
    solver->setAccuracy(0.0, eps/10);
    while (fabs(shotParameter_2 - shotParameter_1)/fabs(shotParameter_1) > eps) {
        initials[0] = (shotParameter_1 + shotParameter_2)/2;
        initials[1] = initials[0];
        solver->setInitials(initials);
        currentSolution = solver->driverApply();
        if (currentSolution.y[0] - phi_0 > 0)
            shotParameter_2 -= (shotParameter_2 - shotParameter_1)/2;
        else shotParameter_1 += (shotParameter_2 - shotParameter_1)/2;
        solver->reset();
    }
    phi_1 = (shotParameter_1 + shotParameter_2)/2;
    calculated = true;
}

double TFPotential::valueAt_1() {
    if (!calculated) calculate();
    return phi_1;
}

double TFPotential::valueAt_x(double x) {
    if (!calculated) calculate();
    double initials[2] = { phi_1, phi_1 };
    double result;
    solver->setInitials(initials);
    solver->setLimits(1.0, x);
    result = solver->driverApply().y[0];
    solver->reset();
    return result;
}

double TFPotential::derivativeAt_0() {
    if (!calculated) calculate();
    double initials[2] = { phi_1, phi_1 };
    double result;
    solver->setInitials(initials);
    solver->setLimits(1.0, 0);
    result = solver->driverApply().y[1];
    solver->reset();
    return result;
}

void TFPotential::setParameters(double V, double T) {
    if (fabs(log10(V) - log10(this->V)) > 1e-10
        || fabs(log10(T) - log10(this->T)) > 1e-10) 
    { 
        calculated = false;
        this->V = V;
        this->T = T;
        *TF::Potential::parameter = pow(2.0, 7.0/6.0)
                                  * pow(3.0, 2.0/3.0)
                                  * pow(M_PI, -5.0/3.0)
                                  * sqrt(T)*pow(V, 2.0/3.0);
        phi_0 = pow(4.0*M_PI/3.0/V, 1.0/3.0)/T;
    }
}

void TFPotential::setPrecision(double eps) {
	/*double phiMin;
	double phiMax;
	double eps1;
	double eps2;
	setInitialShotParameters(phiMin, phiMax);
	eps1 = eps*fabs(2*FermiDirac::ThreeHalf::function(phiMin)/3.0/FermiDirac::Half::function(phiMin));
	eps2 = eps*fabs(2*FermiDirac::ThreeHalf::function(phiMax)/3.0/FermiDirac::Half::function(phiMax));
	if (eps1 > eps2) this->eps = eps2;
	else this->eps = eps1;*/
	this->eps = eps;
}
