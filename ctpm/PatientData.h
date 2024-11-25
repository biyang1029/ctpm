#ifndef PATIENT_DATA_H
#define PATIENT_DATA_H
//struct for all inputs
struct Patient {
    double height;
    double weight;
    double V_o2;
    double V_co2;
    double Bre_r;
    double Ta;
    double Tcor;
};

//different patient information
    extern Patient patient1;
    extern Patient patient2;
    extern Patient patient3;
#endif
