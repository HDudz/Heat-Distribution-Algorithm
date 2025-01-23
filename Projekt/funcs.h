#pragma once
#include <vector>
#include <iostream>
#include "config.h"
#include "structs.h"
using namespace std;


vector<node> pktCalk(int n);


double* getW(int n);

void liczdN(double xi, double eta,
    double& dN1_dxi, double& dN1_deta,
    double& dN2_dxi, double& dN2_deta,
    double& dN3_dxi, double& dN3_deta,
    double& dN4_dxi, double& dN4_deta);

vector<double> liczN(double xi, double eta);

Jakobian liczJakobian(double xi, double eta, node* wezly);

void wypiszJakobian(Jakobian jac);

vector<vector<double>> liczHbc(double alfa, vector<double> len, element elem, double* w);

vector<double> liczP(GlobalData dane, vector<double> len, element elem, double* w);

vector<vector<double>> liczH(GlobalData dane, vector<node> punktyCalkowania, grid siatka, double* w, int nrElem);

void liczGlobalne(vector<vector<double>>& H_global, vector<double>& P_global, vector<vector<double>>& C_global, GlobalData dane, vector<node> punktyCalkowania, grid& siatka, double* w);

vector<double> rozwiazUkladRownan(vector<vector<double>> H, vector<double> P);

vector<double> liczT(
    const vector<vector<double>> H_global,
    const vector<vector<double>> C_global,
    const vector<double> P_global,
    vector<double>& T0,
    GlobalData dane);