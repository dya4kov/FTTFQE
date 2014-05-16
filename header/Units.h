#pragma once
typedef double Charge;
typedef double Mass;

enum Unit {
    dimensionless,
    Atomic,
    eV,
    Kelvin,
    cmc,
    mc,
    Pa,
    GPa,
    MBar,
    gcc,
    gOverCmc,
    kgOverMc,
    oneOverCmc,
    oneOverMc,
    undefinedUnit
};

enum Scaling {
    lin = 0,
    lg = 1
};

const double Avogadro = 6.02214129e+23;