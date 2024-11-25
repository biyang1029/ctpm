#include <iostream>
#include <fstream>
#include <vector>
#include "PatientData.h"
#include "Initialization.h"
#include "Simulation.h"
using namespace std;

int main() {


    // Claim patient data
    Patient* thePatient = &patient1;

    // initializeSimulation
    InitializationData simData = initializeSimulation(thePatient);

    // run the simulation
    runSimulation(thePatient, simData);


    


    return 0;
}
