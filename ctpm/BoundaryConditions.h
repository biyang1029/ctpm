#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H
#include "Initialization.h"
#include "PatientData.h"
#include <vector>

// 定义边界条件的结构体
struct BoundaryConditions {
    double ae, aw, an, as, af, ab, spe, spw, spn, sps, spf, spb, sp, sce, scw, scn, scs, scf, scb, ap0;
};
 


// 函数声明，用于计算边界条件
BoundaryConditions calculateBoundaryConditions(
    int i, int j, int k,
    const InitializationData& simData,
    const Patient* thePatient  // 添加 Patient 指针
   
);



#endif // BOUNDARY_CONDITIONS_H
