#include "Simulation.h"
#include "Initialization.h"
#include "PatientData.h"
#include <iostream>
#include <vector>
#include <fstream>

 std::ofstream outFile("Coretemperature.txt");

 std::vector<std::vector<std::vector<double>>> temperatureField;//Used for temporary storage of temperature data for time node t
std::vector<std::vector<std::vector<double>>> MTtem;//Used for temporary storage of temperature data for time node t-1
double  thermalConductivityField[17][5][9];//Used for temporary storage of temperature data for time node t-1

double flag_Tij;//Decide if get out of the loop
double resi;//Residual error
double ConvectiveCoefficient;//Convective heat transfer coefficient is related to real time body temperature

BoundaryConditions calculateBoundaryConditions(
    int i, int j, int k,

    const InitializationData& simData,
    const Patient* thePatient
    ) {
    BoundaryConditions bc;
      // define the simluation parameters, and the radiation to the environment.

    bc.area = simData.lengthY * simData.WidthChest * 2 + simData.WidthChest * (thePatient->height - simData.DiameterHead)*2 + simData.lengthY * (thePatient->height - simData.DiameterHead) * 2 + 4 * simData.DiameterHead * simData.DiameterHead;
    bc.Qsweat = (16.2 * simData.Qmet + 92.6 * temperatureField[i][j][k] - 3117.16) / bc.area;
    if( ((i == 0 && j == 2) || (i == 0 && j == 3) || (i == 0 && j == 4) || (i == 1 && j == 3) || (i == 1 && j == 4) || (i == 2 && j == 4) ||
        (i == 16 && j == 2) || (i == 16 && j == 3) || (i == 16 && j == 4) || (i == 15 && j == 3) || (i == 15 && j == 4) || (i == 14 && j == 4)) ||((k==0)&&(i!=7&&j!=0)&& (i != 7 && j != 1)&& (i != 7 && j != 2)&& (i != 8 && j != 0)&& (i != 8 && j != 1)&& (i != 8 && j != 2)&& (i != 9 && j != 0) && (i != 9 && j != 1) && (i != 9 && j != 2))){

        bc.ae = bc.aw = bc.an = bc.as = bc.af = bc.ab = 0;
        bc.spe = bc.spw = bc.spn = bc.sps = bc.spf = bc.spb = 0;
        bc.sce = bc.scw = bc.scn = bc.scs = bc.scf = bc.scb = 0;
        bc.Qrade = bc.Qradw = bc.Qradn = bc.Qrads = bc.Qradf = bc.Qradb = 0;
        bc.qsweate = bc.qsweatw = bc.qsweatn  = bc.qsweatf = bc.qsweatb = 0;
        bc.qmetH = bc.qmetT = bc.qmetA = bc.qmetL = 0;
    }
    else {
        bc.ae = ((j == 0 || j == 1) && (i == simData.gridSizeX)) || (j == 2 && i == simData.gridSizeX - 1) || (j == 3 && i == simData.gridSizeX - 2) || (j == 4 && i == simData.gridSizeX - 3)||(i==9&&(j=0||j==1||j==2)&&k==0) ? 0 : (simData.gridSpacingZ[k] * simData.gridSpacingY[j] / simData.partialXe[i][j][k] * thermalConductivityField[i][j][k]);
        bc.aw = ((j == 0 || j == 1) && (i == 0)) || (j == 2 && i == 1) || (j == 3 && i == 2) || (j == 4 && i == 3) || (i == 7 && (j == 0 || j == 1 || j == 2) && k == 0) ? 0 : (simData.gridSpacingZ[k] * simData.gridSpacingY[j] / simData.partialXw[i][j][k] * thermalConductivityField[i][j][k]);
        bc.an = (j == 4) ||((i==7||i==8||i==9)&&j==2&&k==0) ? 0 : (simData.gridSpacingZ[k] * simData.gridSpacingX[i] / simData.partialYn[i][j][k] * thermalConductivityField[i][j][k]);
        bc.as = (j == 0) ? 0 : (simData.gridSpacingZ[k] * simData.gridSpacingX[i] / simData.partialYs[i][j][k] * thermalConductivityField[i][j][k]);
        bc.ab = (k == 1 && (i != 8 && j != 1))||(k==0&&(i==7||i==8||i==9)&&(j==0||j==1||j==2) )? 0 : (simData.gridSpacingY[j] * simData.gridSpacingX[i] / simData.partialZb[i][j][k] * thermalConductivityField[i][j][k]);
        bc.af = (k == 8) ? 0 : (simData.gridSpacingY[j] * simData.gridSpacingX[i] / simData.partialZf[i][j][k] * thermalConductivityField[i][j][k]);


        bc.spe =( (j == 0 || j == 1) && (i == simData.gridSizeX) )|| (j == 2 && i == simData.gridSizeX - 1) || (j == 3 && i == simData.gridSizeX - 2) || (j == 4 && i == simData.gridSizeX - 3) || (i == 9 && (j = 0 || j == 1 || j == 2) && k == 0) ? (1.0 / simData.gridSpacingX[i]) / (1 / ConvectiveCoefficient + simData.partialXe[i][j][k] / thermalConductivityField[i][j][k]) : 0;
        bc.spe = ((j == 0 || j == 1) && (i == 0)) || (j == 2 && i == 1) || (j == 3 && i == 2) || (j == 4 && i == 3) || (i == 7 && (j == 0 || j == 1 || j == 2) && k == 0) ? (1.0 / simData.gridSpacingX[i]) / (1 / ConvectiveCoefficient + simData.partialXw[i][j][k] / thermalConductivityField[i][j][k]) : 0;
        bc.spn = (j == 4) || ((i == 7 || i == 8 || i == 9) && j == 2 && k == 0) ? (1.0 / simData.gridSpacingY[j]) / (1 / ConvectiveCoefficient + simData.partialYn[i][j][k] / thermalConductivityField[i][j][k]) : 0;
        bc.sps = 0;
        bc.spb = (k == 1 && (i != 8 && j != 1)) || (k == 0 && (i == 7 || i == 8 || i == 9) && (j == 0 || j == 1 || j == 2)) ? (1.0 / simData.gridSpacingZ[k]) / (1 / ConvectiveCoefficient + simData.partialZb[i][j][k] / thermalConductivityField[i][j][k]) : 0;
        bc.spf = (k == 8) ? (1.0 / simData.gridSpacingZ[k]) / (1 / ConvectiveCoefficient + simData.partialZf[i][j][k] / thermalConductivityField[i][j][k]) : 0;

        bc.Qrade = ((j == 0 || j == 1) && (i == simData.gridSizeX) )|| (j == 2 && i == simData.gridSizeX - 1) || (j == 3 && i == simData.gridSizeX - 2) || (j == 4 && i == simData.gridSizeX - 3) || (i == 9 && (j = 0 || j == 1 || j == 2) && k == 0) ? ((pow(temperatureField[i][j][k] + 273.15, 4) - pow(19.9 + 273.15, 4)) * 5.56 / pow(10, 8)) / simData.gridSpacingX[i] : 0;
        bc.Qradw = ((j == 0 || j == 1) && (i == 0)) || (j == 2 && i == 1) || (j == 3 && i == 2) || (j == 4 && i == 3) || (i == 7 && (j == 0 || j == 1 || j == 2) && k == 0) ? ((pow(temperatureField[i][j][k] + 273.15, 4) - pow(19.9 + 273.15, 4)) * 5.56 / pow(10, 8)) / simData.gridSpacingX[i] : 0;
        bc.Qradn = (j == 4) || ((i == 7 || i == 8 || i == 9) && j == 2 && k == 0) ? ((pow(temperatureField[i][j][k] + 273.15, 4) - pow(19.9 + 273.15, 4)) * 5.56 / pow(10, 8)) / simData.gridSpacingY[j] : 0;
        bc.Qrads = 0;
        bc.Qradb = (k == 1 && (i != 8 && j != 1)) || (k == 0 && (i == 7 || i == 8 || i == 9) && (j == 0 || j == 1 || j == 2)) ? ((pow(temperatureField[i][j][k] + 273.15, 4) - pow(19.9 + 273.15, 4)) * 5.56 / pow(10, 8)) / simData.gridSpacingZ[k] : 0;
        bc.Qradf = (k == 8) ? ((pow(temperatureField[i][j][k] + 273.15, 4) - pow(19.9 + 273.15, 4)) * 5.56 / pow(10, 8)) / simData.gridSpacingZ[k] : 0;

        bc.qsweate = ((i == 12) && ((j == 0)||(j==1) ||(j == 2)||(j == 3)|| (j == 4))) ? bc.Qsweat / simData.gridSpacingX[i] : 0;
        bc.qsweatw = ((i == 4) && ((j == 0) || (j == 1) || (j == 2) || (j == 3) || (j == 4))) ? bc.Qsweat / simData.gridSpacingX[i] : 0;
        bc.qsweatn = (j == 4) || ((i == 7 || i == 8 || i == 9) && j == 2 && k == 0) ? bc.Qsweat / simData.gridSpacingY[j] : 0;
        bc.qsweatf = (k == 1 && (i != 8 && j != 1)) ? bc.Qsweat / simData.gridSpacingZ[k] : 0;
        bc.qsweatb = (k == 8) || (k == 0 && (i == 7 || i == 8 || i == 9) && (j == 0 || j == 1 || j == 2)) ? bc.Qsweat / simData.gridSpacingZ[k] : 0;
       
        bc.qmetH = ((k == 0 && (i == 7 || i == 8 || i == 9) && (j == 0 || j == 1 || j == 2)))?simData.Qmet * 0.33907 / simData.volumeHead:0;
        bc.qmetT = ((k == 1 || k == 2) && i >= 4 && i <= 12) ? simData.Qmet * 0.099898 / simData.volumeThorax : 0;
        bc.qmetA = ((k == 3 || k == 4) && i >= 4 && i <= 12) ? simData.Qmet * 0.519621 / simData.volumeAbdomen : 0;
        bc.qmetL = ((k >= 5 && k <= simData.gridSizeZ) && i >= 4 && i <= 12) ? simData.Qmet * 0.041411 / simData.volumeLeg : 0;
        bc.qlamp = (i == 8 && (k == 3 || k == 4) && j == 4) ? 262.665 / 4/8/ simData.gridSpacingX[i] / simData.gridSpacingZ[k] / simData.gridSpacingY[j] : 0;
        bc.qlampSide = (i == 8 && j == 4 && (k == 2 || k == 5) || ((i == 7 || i == 9) && j == 4 && (k == 3 || k == 4))) ? 262.665 / 12/8 / simData.gridSpacingX[i] / simData.gridSpacingZ[k] / simData.gridSpacingY[j] : 0;
    }

    bc.Qrad = bc.Qrade + bc.Qradw + bc.Qradn + bc.Qradf + bc.Qradb;
    bc.qsweat = (bc.qsweate + bc.qsweatw + bc.qsweatn + bc.qsweatf + bc.qsweatb)*4.2/3600;
    bc.sp = -bc.spe - bc.spw - bc.spn - bc.sps - bc.spf - bc.spb;
    bc.ap0 = simData.densityField[i][j][k] * simData.specificHeatField[i][j][k] * simData.gridSpacingX[i] * simData.gridSpacingY[j] * simData.gridSpacingZ[k] / simData.DeltaTime;
    bc.qmet = (bc.qmetH + bc.qmetT + bc.qmetA + bc.qmetL);

    return bc;
}
// run the simultion
void runSimulation(Patient* thePatient, InitializationData simData ) {

    temperatureField.resize(simData.gridSizeX+1,
        std::vector<std::vector<double>>(simData.gridSizeY+1,
            std::vector<double>(simData.gridSizeZ+1, 1)));

    MTtem.resize(simData.gridSizeX+1,
        std::vector<std::vector<double>>(simData.gridSizeY+1,
            std::vector<double>(simData.gridSizeZ+1, 1)));
    for (int i = 0; i <= simData.gridSizeX; i++) {
        for (int j = 0; j <= simData.gridSizeY; ++j) {
            for (int k = 0; k <= simData.gridSizeZ; ++k) {
                temperatureField[i][j][k] = simData.temperatureField[i][j][k];
                MTtem[i][j][k] = simData.temperatureField[i][j][k];
               thermalConductivityField[i][j][k] = simData.thermalConductivityField[i][j][k];
            }
        }
    }

    for (int iteration = 0; iteration < simData.MaximumIteration; ++iteration) {

        for (int i = 0; i <= simData.gridSizeX; i++) {
            for (int j = 0; j <= simData.gridSizeY; ++j) {
                for (int k = 0; k <= simData.gridSizeZ; ++k) {

                    MTtem[i][j][k] = temperatureField[i][j][k];
                }
            }
        }

        double Tcore = (temperatureField[8][1][3] + temperatureField[8][1][4]) / 2;
        std::cout << "" << iteration << std::endl;;
        if (1 <= iteration && iteration <= 108000 && (iteration % 600 == 0)) {
            outFile << "time =" << iteration/600 << "min, and Tcore= "<<Tcore << std::endl;;
        }


        //std::cout << "resi =" << resi << std::endl;;
        int reductioncounter = 1;
        int increasecounter = 1;

        //If the core temperature exceeds the limit, the thermal conductivity of the skin is changed
        if (Tcore < 35.5 && !simData.conductivityReductionComplete) {
                // thermal regulation
                reductioncounter++; 
                for (int k = 1; k <= simData.gridSizeZ; ++k) {
                    for (int j = 0; j <= simData.gridSizeY; ++j) {
                        thermalConductivityField[4][j][k] *= 0.99915;
                        thermalConductivityField[12][j][k] *= 0.99915;
                    }
                    for (int i = 5; i < 11; ++i) {
                        flag_Tij = 1;
                        thermalConductivityField[i][4][k] *= 0.99915;
                    }
                    if (reductioncounter >= 600) {
                        simData.conductivityReductionComplete = true; //regulation is done
                    }
                }
            }
        if (Tcore  >= 35.5 && Tcore <= 38) {

            for (int i = 0; i <= simData.gridSizeX; i++) {
                for (int j = 0; j <= simData.gridSizeY; ++j) {
                    for (int k = 0; k <= simData.gridSizeZ; ++k) {

                        thermalConductivityField[i][j][k] = simData.thermalConductivityField[i][j][k];
                    }
                }
            }
        }

        if (Tcore > 38 && !simData.conductivityIncreaseComplete) {
                // thermal regulation
                increasecounter -= 1; 
                for (int k = 1; k < simData.gridSizeZ; ++k) {
                    for (int j = 0; j < simData.gridSizeY; ++j) {
                        thermalConductivityField[4][j][k] *= 0.99915;
                        thermalConductivityField[12][j][k] *= 0.99915;
                    }
                    for (int i = 5; i < 11; ++i) {
                        thermalConductivityField[i][4][k] *= 0.99915;
                    }
                    if (reductioncounter >= 600) {
                        simData.conductivityReductionComplete = true; // regulation is done
                    }
                }
            }
        flag_Tij = 1;
        while (flag_Tij) {
            resi = 0;
            for (int i = 0; i <= simData.gridSizeX; ++i) {
                for (int j = 0; j <= simData.gridSizeY; ++j) {
                    for (int k = 0; k <= simData.gridSizeZ; ++k) {
                        double Tp = temperatureField[i][j][k];
                        double Ttem = Tp;

                        //The convectivecoefficient is different from different part of the body.
                        if (k == 0 && (i == 8||i==7||i==9) && (j == 1||j==0||j==2)) {
                            ConvectiveCoefficient = (temperatureField[i][j][k] - thePatient->Ta)*0.0542+3.5662;
                        }

                        else if (k >= 1 && k <= 2) {
                            ConvectiveCoefficient =0.0096*(temperatureField[i][j][k] - thePatient->Ta)+3.8071;

                        }

                        else if (k >= 3 && k <= 4) {
                            ConvectiveCoefficient =0.0226*(temperatureField[i][j][k] - thePatient->Ta)+3.75;
                        }

                        else if (k >= 5 && k <= 8) {
                            ConvectiveCoefficient = 0.1369*(temperatureField[i][j][k] - thePatient->Ta)+1.4042;
                        }

                        else {
                            ConvectiveCoefficient = 0.00000000001;
                        }

                        BoundaryConditions bc = calculateBoundaryConditions(i, j, k, simData, thePatient);

                        double ap = bc.aw + bc.ae + bc.an + bc.as + bc.af + bc.ab + bc.ap0 - bc.sp * simData.gridSpacingX[i] * simData.gridSpacingY[j] * simData.gridSpacingZ[k];
                        double sou = bc.qmet - bc.Qrad - bc.qsweat + bc.qlamp;
                        double sc = bc.sce + bc.scw + bc.scn + bc.scs + bc.scf + bc.scb + sou;
                        double b = sc * simData.gridSpacingX[i] * simData.gridSpacingY[j] * simData.gridSpacingZ[k] + bc.ap0 * MTtem[i][j][k];
                        double T_left = (i == 0)||(i==7&&(j==0||j==1||j==2)&&k==0) ? temperatureField[i][j][k] : temperatureField[i - 1][j][k];
                        double T_right = (i == simData.gridSizeX) || (i == 9 && (j == 0 || j == 1 || j == 2) && k == 0) ? temperatureField[i][j][k] : temperatureField[i + 1][j][k];
                        double T_down = (j == 0) ? temperatureField[i][j][k] : temperatureField[i][j - 1][k];
                        double T_up = (j == simData.gridSizeY)||((i==7||i==8||i==9)&&j==2&&k==0) ? temperatureField[i][j][k] : temperatureField[i][j + 1][k];
                        double T_back = (k == 0) ? temperatureField[i][j][k] : temperatureField[i][j][k - 1];
                        double T_front = (k == simData.gridSizeZ) ? temperatureField[i][j][k] : temperatureField[i][j][k + 1];

                        // Update the current grid temperature
                        Tp = (bc.aw * T_left +
                            bc.ae * T_right +
                            bc.as * T_down +
                            bc.an * T_up +
                            bc.ab * T_back +
                            bc.af * T_front +
                            b) / ap;
                        resi = resi + fabs(Ttem - Tp);

                        temperatureField[i][j][k] = Tp; // update temperature field

                    }
                }
            }
            flag_Tij = (resi / 795 <= 0.0000001) ? 0 : 1;
        }
    }
// write the file
        outFile.close();

        std::cout << "Temperature data has been successfully written to " << "temperature.txt" << std::endl;
}
