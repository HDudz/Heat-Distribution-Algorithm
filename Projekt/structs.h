#pragma once
#include <vector>
#include<string>
#include <fstream>
#include <iostream>

using namespace std;


struct node {
    double x = -100, y = -100;
    bool onSide = false;
    string wypisz();
};


struct side {
    int ilePC;
    node* pktCalk;
    long double* N[4];

    side();

    string wypisz();
};

struct element {
    int ID[4];
    side boki[4];
    string wypisz();
    string wypisz_boki();
};

struct grid {
    int nN;
    int nE;

    element* ele = new element[nE];
    node* nod = new node[nN];

    grid(int N, int E);

    void odczytaj();
    void ele_gen();
    void czyBok();
    void side_gen();
};

struct GlobalData
{
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int nN;
    int nE;


    void odczytaj();

    GlobalData();

};

struct Jakobian {
    double J[2][2];
    double detJ;
    double invJ[2][2];
};