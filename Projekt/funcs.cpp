#include "funcs.h"
using namespace std;

vector<node> pktCalk(int n)
{
    vector<node> pktCal;
    if (n == 2)
    {
        double p = 1.0 / sqrt(3.0);
        pktCal = {
            {-p,-p}, {p,-p},
            {-p, p}, {p, p}
        };
    }
    else if (n == 3)
    {
        double p = sqrt(3.0 / 5.0);
        pktCal = {
        {-p, -p}, {0.0, -p}, {p, -p},
           {-p, 0.0}, {0.0, 0.0}, {p, 0.0},
        {-p, p}, {0.0, p}, {p, p}
        };
    }
    else if (n == 4)
    {
        double p1 = -0.861136;
        double p2 = -0.339981;
        double p3 = 0.339981;
        double p4 = 0.861136;

        pktCal = {
            {p1, p1}, {p2, p1}, {p3, p1}, {p4, p1},
            {p1, p2}, {p2, p2}, {p3, p2}, {p4, p2},
            {p1, p3}, {p2, p3}, {p3, p3}, {p4, p3},
            {p1, p4}, {p2, p4}, {p3, p4}, {p4, p4}

        };
    }
    return pktCal;
}

double* getW(int n)
{
    double* w = new double[n];
    if (n == 2)
    {
        w[0] = 1;
        w[1] = 1;

    }
    else if (n == 3)
    {
        w[0] = 5.0 / 9.0;
        w[1] = 8.0 / 9.0;
        w[2] = 5.0 / 9.0;
    }
    else if (n == 4) {
        w[0] = 0.347855;
        w[1] = 0.652145;
        w[2] = 0.652145;
        w[3] = 0.347855;
    }
    return w;
}

void liczdN(double xi, double eta, 
    double& dN1_dxi, double& dN1_deta,
    double& dN2_dxi, double& dN2_deta,
    double& dN3_dxi, double& dN3_deta,
    double& dN4_dxi, double& dN4_deta) {
    dN1_dxi = 0.0 - 0.25 * (1 - eta);
    dN1_deta = 0.0 - 0.25 * (1 - xi);

    dN2_dxi = 0.25 * (1 - eta);
    dN2_deta = 0.0 - 0.25 * (1 + xi);

    dN3_dxi = 0.25 * (1 + eta);
    dN3_deta = 0.25 * (1 + xi);

    dN4_dxi = 0.0 - 0.25 * (1 + eta);
    dN4_deta = 0.25 * (1 - xi);
}

vector<double> liczN(double xi, double eta)
{
    vector<double> N(4, 0.0);

    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);
    return N;

}

Jakobian liczJakobian(double xi, double eta, node* wezly) {
    double dN1_dxi, dN1_deta, dN2_dxi, dN2_deta, dN3_dxi, dN3_deta, dN4_dxi, dN4_deta;
    liczdN(xi, eta, dN1_dxi, dN1_deta, dN2_dxi, dN2_deta, dN3_dxi, dN3_deta, dN4_dxi, dN4_deta);

    Jakobian jac;



    jac.J[0][0] = dN1_dxi * wezly[0].x + dN2_dxi * wezly[1].x + dN3_dxi * wezly[2].x + dN4_dxi * wezly[3].x;
    jac.J[0][1] = dN1_dxi * wezly[0].y + dN2_dxi * wezly[1].y + dN3_dxi * wezly[2].y + dN4_dxi * wezly[3].y;
    jac.J[1][0] = dN1_deta * wezly[0].x + dN2_deta * wezly[1].x + dN3_deta * wezly[2].x + dN4_deta * wezly[3].x;
    jac.J[1][1] = dN1_deta * wezly[0].y + dN2_deta * wezly[1].y + dN3_deta * wezly[2].y + dN4_deta * wezly[3].y;

    jac.detJ = jac.J[0][0] * jac.J[1][1] - jac.J[0][1] * jac.J[1][0];

    if (jac.detJ != 0) {
        jac.invJ[0][0] = jac.J[1][1] / jac.detJ;
        jac.invJ[0][1] = -jac.J[0][1] / jac.detJ;
        jac.invJ[1][0] = -jac.J[1][0] / jac.detJ;
        jac.invJ[1][1] = jac.J[0][0] / jac.detJ;
    }

    return jac;
}

void wypiszJakobian(Jakobian jac) {
    cout << "\nMacierz Jacobiego:\n";
    cout << jac.J[0][0] << " " << jac.J[0][1] << "\n";
    cout << jac.J[1][0] << " " << jac.J[1][1] << "\n";
    cout << "\nWyznacznik Jacobiego: " << jac.detJ << "\n";

    if (jac.detJ != 0) {
        cout << "\nOdwrotność Jacobiego:\n";
        cout << jac.invJ[0][0] << " " << jac.invJ[0][1] << "\n"
            << jac.invJ[1][0] << " " << jac.invJ[1][1] << "\n";
    }
    else {
        cout << "Macierz Jacobiego jest osobliwa (nieodwracalna)\n";
    }
}

vector<vector<double>> liczHbc(double alfa, vector<double> len, element elem, double* w) {
    vector<vector<double>> Hbc(4, vector<double>(4, 0.0));


    
    for (int i = 0; i < 4; i++)
    {
        vector<vector<double>> localHbc(4, vector<double>(4, 0.0));
        vector<vector<double>> sumHbc(4, vector<double>(4, 0.0));
        double detJ = len[i] / 2;
        if (elem.boki[i].pktCalk[0].x != -100) {
            side bok = elem.boki[i];
            for (int j = 0; j < punktyCalk; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    for (int l = 0; l < 4; l++)
                    {
                        Hbc[k][l] += alfa * w[j] * bok.N[k][j] * bok.N[l][j] * detJ;
                    }
                }
            }

        }
    }



    return Hbc;
}

vector<double> liczP(GlobalData dane, vector<double> len, element elem, double* w) {

    vector<double> P(4, 0.0);

    for (int i = 0; i < 4; i++)
    {
        double detJ = len[i] / 2.0;
        if (elem.boki[i].pktCalk[0].x != -100) {
            vector<double> sumP(4, 0.0);
            side bok = elem.boki[i];
            for (int j = 0; j < punktyCalk; j++)
            {
                for (int k = 0; k < 4; k++)
                {

                    sumP[k] += dane.Alfa * w[j] * bok.N[k][j] * dane.Tot * detJ;

                }
            }
            P[0] += sumP[0];
            P[1] += sumP[1];
            P[2] += sumP[2];
            P[3] += sumP[3];

        }
    }


    return P;
}



vector<vector<double>> liczH(GlobalData dane, vector<node> punktyCalkowania, grid siatka, double* w, int nrElem) {
    vector<vector<double>> Hbc(4, vector<double>(4, 0.0));
    vector<double> P(4, 0.0);
    vector<vector<double>> H(4, vector<double>(4, 0.0));

    int k = 0;
    int licz = 0;
    double* wLoc = new double[pow(punktyCalk, 2)];
    for (int i = 0; i < punktyCalk; i++)
    {
        for (int j = 0; j < punktyCalk; j++)
        {
            wLoc[licz] = w[i] * w[j];
            licz++;
        }

    }

    node wezly[4];
    wezly[0] = siatka.nod[siatka.ele[nrElem].ID[0] - 1];
    wezly[1] = siatka.nod[siatka.ele[nrElem].ID[1] - 1];
    wezly[2] = siatka.nod[siatka.ele[nrElem].ID[2] - 1];
    wezly[3] = siatka.nod[siatka.ele[nrElem].ID[3] - 1];



    for (const auto& punkt : punktyCalkowania) {
        vector<vector<double>> LocalH(4, vector<double>(4, 0.0));

        
        double xi = punkt.x;
        double eta = punkt.y;
        if(debug)
        {
            cout << endl;
            for (int i = 0; i < 70; i++)
            {
                cout << "=";
            }
            cout << "\nLiczenie dla punktu całkowania " << k + 1 << ":\txi " << xi << ", eta " << eta << "\n";


            for (int i = 0; i < 70; i++)
            {
                cout << "=";
            }
            cout << "\n";
        }

        k++;

        Jakobian jac = liczJakobian(xi, eta, wezly);

        if(debug){
            wypiszJakobian(jac);
        }
        if (jac.detJ == 0) continue; // Pomijamy przypadki, gdzie Jakobian jest osobliwy


        double dN_dx[4];
        double dN_dy[4];

        double dN1_dxi, dN1_deta, dN2_dxi, dN2_deta, dN3_dxi, dN3_deta, dN4_dxi, dN4_deta;
        liczdN(xi, eta, dN1_dxi, dN1_deta, dN2_dxi, dN2_deta, dN3_dxi, dN3_deta, dN4_dxi, dN4_deta);



        dN_dx[0] = jac.invJ[0][0] * dN1_dxi + jac.invJ[0][1] * dN1_deta;
        dN_dy[0] = jac.invJ[1][0] * dN1_dxi + jac.invJ[1][1] * dN1_deta;

        dN_dx[1] = jac.invJ[0][0] * dN2_dxi + jac.invJ[0][1] * dN2_deta;
        dN_dy[1] = jac.invJ[1][0] * dN2_dxi + jac.invJ[1][1] * dN2_deta;

        dN_dx[2] = jac.invJ[0][0] * dN3_dxi + jac.invJ[0][1] * dN3_deta;
        dN_dy[2] = jac.invJ[1][0] * dN3_dxi + jac.invJ[1][1] * dN3_deta;

        dN_dx[3] = jac.invJ[0][0] * dN4_dxi + jac.invJ[0][1] * dN4_deta;
        dN_dy[3] = jac.invJ[1][0] * dN4_dxi + jac.invJ[1][1] * dN4_deta;

        

        if(debug)
        {
            cout << "\ndN/dx równa się:\n";
            cout << dN_dx[0] << ", " << dN_dx[1] << ", " << dN_dx[2] << ", " << dN_dx[3] << endl;
            cout << "dN/dy równa się:\n";
            cout << dN_dy[0] << ", " << dN_dy[1] << ", " << dN_dy[2] << ", " << dN_dy[3] << "\n\n";
            cout << "\nMacierz H lokalna:\n";
        }


        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                LocalH[i][j] = (dane.Conductivity * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]) * jac.detJ) * wLoc[k - 1];

                if (debug) { cout << LocalH[i][j] << " \t"; }
                H[i][j] += LocalH[i][j];
            }
            if(debug) cout << "\n";
        }


    }

    


    return H;
}

vector<vector<double>> liczC(GlobalData dane, vector<node> punktyCalkowania, grid siatka, double* w, int nrElem) {
    vector<vector<double>> C(4, vector<double>(4, 0.0));
    int k = 0;
    int licz = 0;
    double* wLoc = new double[pow(punktyCalk, 2)];
    for (int i = 0; i < punktyCalk; i++)
    {
        for (int j = 0; j < punktyCalk; j++)
        {
            wLoc[licz] = w[i] * w[j];
            licz++;
        }

    }

    node wezly[4];
    wezly[0] = siatka.nod[siatka.ele[nrElem].ID[0] - 1];
    wezly[1] = siatka.nod[siatka.ele[nrElem].ID[1] - 1];
    wezly[2] = siatka.nod[siatka.ele[nrElem].ID[2] - 1];
    wezly[3] = siatka.nod[siatka.ele[nrElem].ID[3] - 1];



    for (const auto& punkt : punktyCalkowania) {
        vector<vector<double>> LocalC(4, vector<double>(4, 0.0));


        double xi = punkt.x;
        double eta = punkt.y;
        if (debug)
        {
            cout << endl;
            for (int i = 0; i < 70; i++)
            {
                cout << "=";
            }
            cout << "\nLiczenie dla punktu całkowania " << k + 1 << ":\txi " << xi << ", eta " << eta << "\n";


            for (int i = 0; i < 70; i++)
            {
                cout << "=";
            }
            cout << "\n";
        }

        k++;

        Jakobian jac = liczJakobian(xi, eta, wezly);

        if (debug) {
            wypiszJakobian(jac);
        }
        if (jac.detJ == 0) continue;


        vector<double> N = liczN(xi, eta);


        


        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                C[i][j] += (dane.SpecificHeat * dane.Density * jac.detJ * (N[i] * N[j])) * wLoc[k - 1];

                if (debug) { cout << LocalC[i][j] << " \t"; }


                C[i][j] = round(C[i][j] * 10000.0) / 10000.0;

            }
            if (debug) cout << "\n";
        }


    }



    return C;
}




void liczGlobalne(vector<vector<double>> &H_global, vector<double> &P_global, vector<vector<double>>& C_global,GlobalData dane, vector<node> punktyCalkowania, grid& siatka, double* w) {
    vector<double> len(4);


    for (int k = 0; k < siatka.nE; k++) {
        cout << "\n\n\n\n";
        for (int i = 0; i < 120; i++)
        {
            cout << "=";
        }
        cout << "\nELEMENT " << k + 1 << endl;
        cout << "Jego boki:" << "\n" << siatka.ele[k].wypisz_boki();
        for (int i = 0; i < 120; i++)
        {
            cout << "=";
        }
        cout << "\n";
        int* local_ids = siatka.ele[k].ID;
        for (int i = 0; i < 4; i++)
        {
            int id1 = i % 4;
            int id2 = (i + 1) % 4;
            len[i] = sqrt(pow(siatka.nod[local_ids[id1] - 1].x - siatka.nod[local_ids[id2] - 1].x, 2) + pow(siatka.nod[local_ids[id1] - 1].y - siatka.nod[local_ids[id2] - 1].y, 2));
        }
        

        vector<vector<double>> H_local = liczH(dane, punktyCalkowania, siatka, w, k);
        vector<vector<double>> Hbc = liczHbc(dane.Alfa, len, siatka.ele[k], w);
        vector<vector<double>> C = liczC(dane, punktyCalkowania, siatka, w, k);
        vector<double> P = liczP(dane, len, siatka.ele[k], w);

        cout << "\n\nMacierz H:\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << H_local[i][j] << "\t";
            }
            cout << "\n";
        }


        cout << "\n\nMacierz Hbc elementu:\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << Hbc[i][j] << " \t";
            }
            cout << "\n";
        }


        cout << "\n\nMacierz H+Hbc:\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                H_local[i][j] += Hbc[i][j];
                cout << H_local[i][j] << "\t";
            }
            cout << "\n";
        }


        cout << "\n\nMacierz C:\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                cout << C[i][j] << "\t";
            }
            cout << "\n";
        }


        cout << "\n\nWektor P elementu:\n";
        for (double wartosc : P) {
            cout << wartosc << " \t";
        }
        cout << "\n";
        

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                H_global[local_ids[i] - 1][local_ids[j] - 1] += H_local[i][j];
                C_global[local_ids[i] - 1][local_ids[j] - 1] += C[i][j];
            }
            P_global[local_ids[i] - 1] += P[i];
        }

    }
    

    return;
}

vector<double> rozwiazUkladRownan(vector<vector<double>> H, vector<double> P) {
    int n = H.size();

    //eliminacja gausa
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double factor = H[j][i] / H[i][i];
            for (int k = i; k < n; ++k) {
                H[j][k] -= factor * H[i][k];
            }
            P[j] -= factor * P[i];
        }
    }

    //podstawiam wstecznie
    vector<double> T(n);
    for (int i = n - 1; i >= 0; --i) {
        T[i] = P[i];
        for (int j = i + 1; j < n; ++j) {
            T[i] -= H[i][j] * T[j];
        }
        T[i] /= H[i][i];
    }

    return T;
}

vector<double> liczT(
    const vector<vector<double>> H_global,
    const vector<vector<double>> C_global,
    const vector<double> P_global,
    vector<double>& T0,
    GlobalData dane) {

    vector<vector<double>> A(dane.nN, vector<double>(dane.nN, 0.0));
    double min, max;
    if (debug) cout << "\nObliczone [H]+[C]/dT: \n";

    for (int i = 0; i < dane.nN; i++) {
        for (int j = 0; j < dane.nN; j++) {
            A[i][j] = H_global[i][j] + C_global[i][j] / dane.SimulationStepTime;
            if (debug)cout << A[i][j] << "  ";
        }
        if (debug)cout << "\n";
    }

    if (!debug)cout << "Czas:\tmin:\t\tmax:\n";
    for (int czas = 0; czas <= dane.SimulationTime; czas+= dane.SimulationStepTime)
    {
        if (debug){
            for (int i = 0; i < 60; i++)
            {
                cout << "===";
            }
            cout << "\n\nDla czasu = " << czas << ":\n";
        }

        vector<double> b(dane.nN, 0.0);
        if (debug) cout << "\n\nObliczone {P}+{[C]/dT}*{T0}: \n";
        for (int i = 0; i < dane.nN; i++) {
            for (int j = 0; j < dane.nN; j++) {
                b[i] += (C_global[i][j] / dane.SimulationStepTime) * T0[j];
            }
            b[i] += P_global[i];
            if(debug) cout << b[i] << "  ";
        }



        if(debug){
            cout << "\nTemperatura: \n";

            for (int i = 0; i < dane.nN; i++)
            {
                cout << T0[i] << "  ";
            }
        }


        min = *min_element(T0.begin(), T0.end());
        max = *max_element(T0.begin(), T0.end());

        if(debug){
            cout << "\nMin Temp:\t" << min << "\n";
            cout << "Max Temp:\t" << max << "\n\n";
        }


        if (!debug) cout << czas << "\t" << min << "\t" << max << "\n";


        T0 = rozwiazUkladRownan(A, b);


    }
    return T0;
}
