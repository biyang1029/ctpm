#include "Initialization.h"
#include "PatientData.h"
#include <vector>
#include <iostream>

// Function definition: Use the passed patient pointer to initialize the data
InitializationData initializeSimulation(Patient* thePatient) {
    // the maximum number of grids in the xyz direction
    InitializationData simData;
    simData.gridSizeX = 16;
    simData.gridSizeY = 4;
    simData.gridSizeZ = 8;
    //initialization of the matrix of temperature field,specificHeatField,densityField,thermalConductivityField...
    simData.temperatureField.resize(simData.gridSizeX+1,
        std::vector<std::vector<double>>(simData.gridSizeY+1,
            std::vector<double>(simData.gridSizeZ+1, thePatient->Ta + 10)));
    simData.specificHeatField.resize(simData.gridSizeX+1,
        std::vector<std::vector<double>>(simData.gridSizeY+1,
            std::vector<double>(simData.gridSizeZ+1, 1.2)));
    simData.densityField.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));
    simData.thermalConductivityField.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));
    simData.partialXe.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));
    simData.partialXw.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));
    simData.partialYn.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));
    simData.partialYs.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));
    simData.partialZf.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));
    simData.partialZb.resize(simData.gridSizeX + 1,
        std::vector<std::vector<double>>(simData.gridSizeY + 1,
            std::vector<double>(simData.gridSizeZ + 1, 2)));

    // initialization of the temperature field
    for (int i = 0; i <= simData. gridSizeX; ++i) {
        for (int j = 0; j <= simData.gridSizeY; ++j) {
            for (int k = 0; k <= simData.gridSizeZ; ++k) {
                simData.temperatureField[i][j][k] =  thePatient->Ta ;
            }
        }
    }

    for (int i = 3; i <= simData.gridSizeX-3; ++i) {
        for (int j = 0; j <= simData.gridSizeY; ++j) {
            for (int k = 0; k <= simData.gridSizeZ; ++k) {
                simData.temperatureField[i][j][k] = 0.6 * thePatient->Ta + 0.7 * thePatient->Tcor;
            }
        }
    }
    for (int i = 5; i <= simData.gridSizeX - 5; ++i) {
        for (int j = 0; j <= simData.gridSizeY - 1; ++j) {
            for (int k = 0; k <= simData.gridSizeZ; ++k) {
                simData.temperatureField[i][j][k] = 0.5 * thePatient->Ta + 0.5 * thePatient->Tcor;
            }
        }
    }
    for (int i = 7; i <= simData.gridSizeX - 7; ++i) {
        for (int j = 0; j <= simData.gridSizeY - 2; ++j) {
            for (int k = 1; k <= simData.gridSizeZ-1; ++k) {
                simData.temperatureField[i][j][k] = 0.25 * thePatient->Ta + 0.75 * thePatient->Tcor;
            }
        }
    }
    for (int i = 8; i <= simData.gridSizeX - 8; ++i) {
        for (int j = 0; j <= simData.gridSizeY - 3; ++j) {
            for (int k = 2; k <= simData.gridSizeZ - 2; ++k) {
                simData.temperatureField[i][j][k] = 0.1*thePatient->Ta +0.9 *thePatient->Tcor;
            }
        }
    }
    simData.temperatureField[8][1][3] = simData.temperatureField[8][1][4]= thePatient->Tcor;

    // calculate the body size parameters
    simData.WidthChest = 0.386 + 0.0015 * (thePatient->weight / thePatient->height / thePatient->height - 21.75);
    simData.DiameterHead = (0.422 + 0.08673 * thePatient->height) / 4;
    //meshing
    
    simData.gridSpacingX = {
        (0.6 - simData.WidthChest - 0.06) / 8,
        (0.6 - simData.WidthChest - 0.06) / 8,
        (0.6 - simData.WidthChest - 0.06) / 8,
        (0.6 - simData.WidthChest - 0.06) / 8,
        0.03,
        (simData.WidthChest - simData.DiameterHead) / 4, 
        (simData.WidthChest - simData.DiameterHead) / 4,
        (simData.DiameterHead - 0.09) / 2,
        0.09,
        (simData.DiameterHead - 0.09) / 2,
        (simData.WidthChest - simData.DiameterHead) / 4,
        (simData.WidthChest - simData.DiameterHead) / 4,
        0.03,
        (0.6 - simData.WidthChest - 0.06) / 8,
        (0.6 - simData.WidthChest - 0.06) / 8,
        (0.6 - simData.WidthChest - 0.06) / 8,
        (0.6 - simData.WidthChest - 0.06) / 8
    };

    simData.gridSpacingY = { 
        (simData.DiameterHead - 0.09) / 2,
        0.09,
        (simData.DiameterHead - 0.09) / 2,
        thePatient->weight / 985 / (thePatient->height - simData.DiameterHead) / simData.WidthChest - 0.03,
        0.03 
    };

    simData.gridSpacingZ = {
        simData.DiameterHead,
        (thePatient->height - simData.DiameterHead) / 8,
        (thePatient->height - simData.DiameterHead) / 8,
        (thePatient->height - simData.DiameterHead) / 8,
        (thePatient->height - simData.DiameterHead) / 8,
        (thePatient->height - simData.DiameterHead) / 8,
        (thePatient->height - simData.DiameterHead) / 8,
        (thePatient->height - simData.DiameterHead) / 8,
        (thePatient->height - simData.DiameterHead) / 8 ,
    };

    //Body thickness
    simData.lengthY = 0.0;
    for (double spacing : simData.gridSpacingY) {
        simData.lengthY += spacing;
    }

    //Set a different specific heat density and thermal conductivity for each node
    for (int j = 0; j <= simData.gridSizeY; ++j) {
        for (int i = 4; i <= simData.gridSizeX - 4; ++i) {
            for (int k = 1; k <= 2; ++k) {
                simData.specificHeatField[i][j][k] = 3021;
                simData.densityField[i][j][k] = 980;
                simData.thermalConductivityField[i][j][k] = 0.387 + 0.13684;
            }
            for (int k = 3; k <= 4; ++k) {
                simData.specificHeatField[i][j][k] = 3018;
                simData.densityField[i][j][k] = 1045;
                simData.thermalConductivityField[i][j][k] = 0.423 + 0.13684;
            }
            for (int k = 5; k <= simData.gridSizeZ; ++k) {
                simData.specificHeatField[i][j][k] = 3028;
                simData.densityField[i][j][k] = 1057;
                simData.thermalConductivityField[i][j][k] = 0.407 + 0.13684;
            }
        }
    }

    for (int j = 0; j <= simData.gridSizeY; ++j) {
        for (int i = 0; i <= 3; ++i) {
            for (int k = 0; k <= simData.gridSizeZ; ++k) {
                simData.specificHeatField[i][j][k] = 1030;
                simData.densityField[i][j][k] = 1.225;
                simData.thermalConductivityField[i][j][k] = 0.026;
            }
        }
        for (int i = simData.gridSizeX - 3; i <= simData.gridSizeX; ++i) {
            for (int k = 0; k <= simData.gridSizeZ; ++k) {
                simData.specificHeatField[i][j][k] = 1030;
                simData.densityField[i][j][k] = 1.225;
                simData.thermalConductivityField[i][j][k] = 0.026;
            }
        }
    }
    for (int j = 0; j <= simData.gridSizeY; ++j) {
        for (int i = 0; i <= simData.gridSizeX; ++i) {
            simData.specificHeatField[i][j][0] = 1030;
            simData.densityField[i][j][0] = 1.225;
            simData.thermalConductivityField[i][j][0] = 0.000001;      
        }
    }
    //Set a small thermal conductivity for areas of the grid that are not involved in the calculation
    for (int k = 1; k <= simData.gridSizeZ; ++k) {
        simData.thermalConductivityField[0][2][k] = simData.thermalConductivityField[0][3][k] = simData.thermalConductivityField[0][4][k]
            = simData.thermalConductivityField[1][3][k] = simData.thermalConductivityField[1][4][k] = simData.thermalConductivityField[2][4][k] =
            simData.thermalConductivityField[simData.gridSizeX][2][k] = simData.thermalConductivityField[simData.gridSizeX][3][k] = simData.thermalConductivityField[simData.gridSizeX][4][k]
            = simData.thermalConductivityField[simData.gridSizeX - 1][3][k] = simData.thermalConductivityField[simData.gridSizeX - 1][4][k] = simData.thermalConductivityField[simData.gridSizeX - 2][4][k] =
            0.000001;
    }
    //Density thermal conductivity and thermal conductivity of the head grid
    for (int i = 7; i <= 9; ++i) {
        for (int j = 0; j <= 2; ++j) {
            simData.specificHeatField[i][j][0] = 3257;
            simData.densityField[i][j][0] = 1109;
            simData.thermalConductivityField[i][1][0] = 0.533 + 0.136849;
        }
    }

    //Initialize partial x partial y partial z
    for (int k = 0; k <= simData.gridSizeZ; ++k) {
        for (int j = 0; j <= simData.gridSizeY; ++j) {
            for (int i = 0; i <= simData.gridSizeX - 1; ++i) {
                simData.partialXe[i][j][k] = (simData.gridSpacingX[i] / 2.0) + (simData.gridSpacingX[i + 1] / 2.0);
            }
            simData.partialXe[simData.gridSizeX][j][k] = (simData.gridSpacingX[simData.gridSizeX] / 2.0);
        }
        simData.partialXe[13][4][k] = simData.partialXe[14][3][k] = simData.partialXe[15][2][k] = simData.partialXe[16][1][k] = 0.001;
    }
    simData.partialXe[9][0][0]= simData.partialXe[9][1][0] = simData.partialXe[9][2][0] = (simData.gridSpacingX[9] / 2.0);

    for (int k = 0; k <= simData.gridSizeZ; ++k) {
        for (int j = 0; j <= simData.gridSizeY; ++j) {
            for (int i = 1; i <= simData.gridSizeX; ++i) {
                simData.partialXw[i][j][k] = (simData.gridSpacingX[i] / 2.0) + (simData.gridSpacingX[i - 1] / 2.0);
            }
            simData.partialXw[0][j][k] = (simData.gridSpacingX[0] / 2.0);
        }
        simData.partialXw[0][1][k] = simData.partialXw[1][2][k] = simData.partialXw[2][3][k] = simData.partialXw[3][4][k] = 0.001;
    }
    simData.partialXw[7][0][0] = simData.partialXw[7][1][0] = simData.partialXw[7][2][0] = (simData.gridSpacingX[7] / 2.0);

    for (int k = 0; k <= simData.gridSizeZ; ++k) {
        for (int i = 0; i <= simData.gridSizeX; ++i) {
            for (int j = 0; j <= simData.gridSizeY-1; ++j) {
                simData.partialYn[i][j][k] = (simData.gridSpacingY[j] / 2.0) + (simData.gridSpacingY[j + 1] / 2.0);
            }
            simData.partialYn[i][simData.gridSizeY][k] = (simData.gridSpacingY[simData.gridSizeY] / 2.0);
        }
        simData.partialYn[0][1][k] = simData.partialYn[1][2][k] = simData.partialYn[2][3][k] = simData.partialYn[3][4][k]
            = simData.partialYn[16][1][k] = simData.partialYn[15][2][k] = simData.partialYn[14][3][k] = simData.partialYn[13][4][k]
            = 0.001;
    }
    simData.partialYn[7][2][0] = simData.partialYn[8][2][0] = simData.partialYn[9][2][0] = (simData.gridSpacingY[2] / 2.0);

    for (int k = 0; k <= simData.gridSizeZ; ++k) {
        for (int i = 0; i <= simData.gridSizeX; ++i) {
            for (int j = 1; j <= simData.gridSizeY; ++j) {
                simData.partialYs[i][j][k] = (simData.gridSpacingY[j] / 2.0) + (simData.gridSpacingY[j - 1] / 2.0);
            }
            simData.partialYs[i][0][k] = (simData.gridSpacingY[0] / 2.0);
        }
;    }
    simData.partialYs[7][0][0] = simData.partialYs[8][0][0] = simData.partialYs[9][0][0] = (simData.gridSpacingY[0] / 2.0);

    for (int i = 0; i <= simData.gridSizeX; ++i) {
        for (int j = 0; j <= simData.gridSizeY; ++j) {
            for (int k = 0; k <= simData.gridSizeZ-1; ++k) {
                simData.partialZf[i][j][k] = (simData.gridSpacingZ[k] / 2.0) + (simData.gridSpacingZ[k + 1] / 2.0);
            }
            simData.partialZf[i][j][simData.gridSizeZ] = (simData.gridSpacingZ[simData.gridSizeZ] / 2.0);
        }
    }

    for (int i = 0; i <= simData.gridSizeX; ++i) {
        for (int j = 0; j <= simData.gridSizeY; ++j) {
            for (int k = 1; k <= simData.gridSizeZ ; ++k) {
                simData.partialZb[i][j][k] = (simData.gridSpacingZ[k] / 2.0) + (simData.gridSpacingZ[k - 1] / 2.0);
            }
            simData.partialZb[i][j][0] = (simData.gridSpacingZ[0] / 2.0);
        }
    }
    simData.partialZb[7][0][0] = simData.partialZb[8][0][0] = simData.partialZb[9][0][0] =
        simData.partialZb[7][1][0] = simData.partialZb[8][1][0] = simData.partialZb[9][1][0]=
        simData.partialZb[7][2][0] = simData.partialZb[8][2][0] = simData.partialZb[9][2][0]
        = (simData.gridSpacingZ[0] / 2.0);


    simData.RadiationPower = 226.225;
    simData.DeltaTime = 0.1;
    simData.MaximumIteration = 108000;

    simData.Qmet = ((3.941 * thePatient->V_o2 / 100 * thePatient->Bre_r + 1.106 * thePatient->V_co2 / 100 * thePatient->Bre_r - 0.061 * 2.17) * 4.184 * 1000 / 60 - 1.76 * thePatient->Bre_r/ 1000 * (thePatient->Tcor - (thePatient->Ta - (0.370 * 0.81 * 4400))));

    simData.volumeHead = simData.DiameterHead * simData.DiameterHead* simData.DiameterHead;
    simData.volumeAbdomen = simData.WidthChest * simData.lengthY * (simData.gridSpacingZ[1] + simData.gridSpacingZ[2]);
    simData.volumeThorax = simData.WidthChest * simData.lengthY * (simData.gridSpacingZ[3] + simData.gridSpacingZ[4]);
    simData.volumeLeg = simData.WidthChest * simData.lengthY * (simData.gridSpacingZ[5] + simData.gridSpacingZ[6] + simData.gridSpacingZ[7] + simData.gridSpacingZ[8]);
    //It is used to determine whether the number of thermal adjustment steps is reached
    simData.conductivityReductionComplete = false;
    simData.conductivityIncreaseComplete = false;
    return simData;

}