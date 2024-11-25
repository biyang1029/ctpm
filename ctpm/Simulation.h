#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H

#include "PatientData.h"
#include "Initialization.h"

//  function for initialization 
void runSimulation(Patient* thePatient, InitializationData simData);

// stuct for boundary condition
struct BoundaryConditions {
    //parameters for calcilation
    double ae, aw, an, as, af, ab, spe, spw, spn, sps, spf, spb, sp, sce, scw, scn, scs, scf, scb, ap0, Qrad, Qrade, Qradw, Qradn, Qrads, Qradf, Qradb, Qsweat,qsweat, qsweate, qsweatw, qsweatn, qsweatf, qsweatb;
    double qmet;
    double qmetH, qmetT, qmetA, qmetL;
    double area;
    double qlamp,qlampSide;
};

// function for calculate boundary conditions
BoundaryConditions calculateBoundaryConditions(
    int i, int j, int k,
    const InitializationData& simData,
    const Patient* thePatient  // ÃÌº” Patient ÷∏’Î

);


#endif // SIMULATION_H
