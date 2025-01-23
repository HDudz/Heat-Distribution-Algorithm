#include "funcs.h"
#include "structs.h"

int main() {

    cout << fixed;
    cout.precision(5);

    GlobalData dane;
    grid siatka(dane.nN, dane.nE);

    cout << "Nodes:\n";
    for (int i = 0; i < dane.nN; i++)
    {
        cout << i + 1 << ".\t" << siatka.nod[i].x << ", \t" << siatka.nod[i].y << " \t" << siatka.nod[i].onSide << endl;
    }

    cout << "Elements:\n";
    for (int i = 0; i < dane.nE; i++)
    {
        cout << i + 1 << ".\t(" << siatka.ele[i].wypisz() << ")\n";

    }
  


    double* w = getW(punktyCalk);

    vector<node> punktyCalkowania = pktCalk(punktyCalk);




    vector<double> P_global(dane.nN, 0.0);
    vector<vector<double>> H_global(dane.nN, vector<double>(dane.nN, 0.0));
    vector<vector<double>> C_global(dane.nN, vector<double>(dane.nN, 0.0));
    vector<double> T0(dane.nN, dane.InitialTemp);
    vector<double> b(dane.nN, 0.0);
    liczGlobalne(H_global, P_global, C_global , dane, punktyCalkowania, siatka, w);

    if(data_path != "dane/dane_31x31.txt"){
        cout << "\n\n\nMacierz globalna H:\n";
        for (auto& wiersz : H_global) {
            for (double wartosc : wiersz) {
                cout << wartosc << "   ";
            }
            cout << "\n";
        }

        cout << "\n\n\nMacierz globalna C:\n";
        for (auto& wiersz : C_global) {
            for (double wartosc : wiersz) {
                cout << wartosc << "   ";
            }
            cout << "\n";
        }


        cout << "\n\nWektor globalny P:\n";

        for (double wartosc : P_global) {
            cout << wartosc << "  ";
        }
        cout << "\n";

        cout << "\n\n\n";
    }

    vector<double> stacjonarne(dane.nN);

    stacjonarne = rozwiazUkladRownan(H_global, P_global);

    cout << "\n\nWektor Temperatur rozwiązania stacjonarnego:\n";

    for (double wartosc : stacjonarne) {
        cout << wartosc << "  ";
    }
    cout << "\n\n";

   

    liczT(H_global, C_global, P_global, T0, dane);

    return 0;
}