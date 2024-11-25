#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "PatientData.h"
#include <vector>

// Define a struct to store initialization data
struct InitializationData {
    int gridSizeX ;
    int gridSizeY;
    int gridSizeZ;

    // Define a struct to store the 3D temperature field data
    std::vector<std::vector<std::vector<double>>> temperatureField;

    double WidthChest;
    double DiameterHead;

    // Matrix for storing dx, dy, dz
    std::vector<double> gridSpacingX;
    std::vector<double> gridSpacingY;
    std::vector<double> gridSpacingZ;

    // A 3D matrix for storing information such as specific heat and density.
    std::vector<std::vector<std::vector<double>>> specificHeatField;
    std::vector<std::vector<std::vector<double>>> densityField;
    std::vector<std::vector<std::vector<double>>> thermalConductivityField;
    //Define a 3D matrix to store partial derivatives such as ?x, ?y, and ?z
    std::vector<std::vector<std::vector<double>>> partialXe;
    std::vector<std::vector<std::vector<double>>> partialXw;
    std::vector<std::vector<std::vector<double>>> partialYn;
    std::vector<std::vector<std::vector<double>>> partialYs;
    std::vector<std::vector<std::vector<double>>> partialZf;
    std::vector<std::vector<std::vector<double>>> partialZb;
    
    double RadiationPower;
    double DeltaTime;
    int MaximumIteration;
    //metabolic rate
    double Qmet;
    //size of the body for future calculation
    double volumeHead,volumeThorax,volumeAbdomen,volumeLeg;
    double lengthY;
    bool conductivityReductionComplete = false;
    bool conductivityIncreaseComplete = false;
};

InitializationData initializeSimulation(Patient* thePatient);

#endif 
