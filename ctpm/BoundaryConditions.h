#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H
#include "Initialization.h"
#include "PatientData.h"
#include <vector>

// ����߽������Ľṹ��
struct BoundaryConditions {
    double ae, aw, an, as, af, ab, spe, spw, spn, sps, spf, spb, sp, sce, scw, scn, scs, scf, scb, ap0;
};
 


// �������������ڼ���߽�����
BoundaryConditions calculateBoundaryConditions(
    int i, int j, int k,
    const InitializationData& simData,
    const Patient* thePatient  // ��� Patient ָ��
   
);



#endif // BOUNDARY_CONDITIONS_H
