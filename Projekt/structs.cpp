#include "config.h"
#include "structs.h"
#include "funcs.h"

string node::wypisz()
{
	string txt = "(" + to_string(x) + ", " + to_string(y) + ")";
	return txt;
}

side::side() 
{

    ilePC = punktyCalk;
    pktCalk = new node[ilePC];
    N[0] = new long double[ilePC];
    N[1] = new long double[ilePC];
    N[2] = new long double[ilePC];
    N[3] = new long double[ilePC];
    
}

string side::wypisz()
{
    string txt = "";
    if (pktCalk->x == -100)
    {
        return "Wewnętrzny.\n";
    }
    for (int i = 0; i < ilePC; i++)
    {
        txt += pktCalk[i].wypisz() + "\t" + to_string(N[0][i]) + "\t" + to_string(N[1][i]) + "\t" + to_string(N[2][i]) + "\t" + to_string(N[3][i]) + "\n";
    }
    return txt;
}

string element::wypisz()
{
    string txt = to_string(ID[0]) + ", " + to_string(ID[1]) + ", " + to_string(ID[2]) + ", " + to_string(ID[3]);
    return txt;
}

string element::wypisz_boki() {
    string txt = "";
    for (int i = 0; i < 4; i++)
    {
        txt += "Bok " + to_string(i + 1) + ":\n" + boki[i].wypisz();
    }
    return txt;
}


grid::grid(int N, int E) : nN(N), nE(E) {
    odczytaj();
    ele_gen();
    czyBok();
    side_gen();
}

void grid::odczytaj() {
    string line;
    ifstream file1(data_path);
    int licznik = 0;

    for (int i = 0; i < 10 && getline(file1, line); i++) {}

    while (getline(file1, line) && licznik < this->nN) {

        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        size_t first_comma = line.find(',');
        size_t second_comma = line.find(',', first_comma + 1);

        if (first_comma != string::npos && second_comma != string::npos) {
            this->nod[licznik].x = stod(line.substr(first_comma + 1, second_comma - first_comma - 1));
            this->nod[licznik].y = stod(line.substr(second_comma + 1));

            licznik++;
        }
    }
}


void grid::ele_gen()
{
    int licz = 0;
    for (int i = 0; i < this->nE; i++)
    {
        licz++;
        if (licz % int(sqrt(this->nN)) == 0)
        {
            licz++;
        }
        this->ele[i].ID[0] = licz;
        this->ele[i].ID[1] = licz + 1;
        this->ele[i].ID[2] = licz + sqrt(this->nN) + 1;
        this->ele[i].ID[3] = licz + sqrt(this->nN);


    }
}



void grid::czyBok() {
    int sideLength = static_cast<int>(sqrt(nN));


    node LD = nod[0];
    node LG = nod[nN - sideLength];
    node PD = nod[sideLength - 1];
    node PG = nod[nN - 1];


    for (int i = 0; i < nN; ++i) {
        nod[i].onSide = false;

        if (nod[i].x >= min(LD.x, LG.x) && nod[i].x <= max(LD.x, LG.x))//czy na lewym boku
        {
            nod[i].onSide = true;
        }
        else if (nod[i].x >= min(PD.x, PG.x) && nod[i].x <= max(PD.x, PG.x))//czy na prawym boku
        {
            nod[i].onSide = true;
        }
        else if (nod[i].y >= min(LD.y, PD.y) && nod[i].y <= max(LD.y, PD.y))//czy na dolnym boku
        {
            nod[i].onSide = true;
        }
        else if (nod[i].y >= min(LG.y, PG.y) && nod[i].y <= max(LG.y, PG.y))//czy na górnym boku
        {
            nod[i].onSide = true;
        }

    }
}

void grid::side_gen()
{
    element elem;
    vector<node>pCal = pktCalk(punktyCalk);
    for (int i = 0; i < nE; i++)
    {
        elem = ele[i];
        for (int j = 0; j < 4; j++)
        {
            int nextJ = (j + 1) % 4;
            for (int k = 0; k < punktyCalk; k++)
            {
                if (nod[elem.ID[j] - 1].onSide && nod[elem.ID[nextJ] - 1].onSide && (j + 1 % 2)) {

                    elem.boki[j].pktCalk[k].x = pCal[k].x;
                    elem.boki[j].pktCalk[k].y = -1 + (j % 3);

                }
                if (nod[elem.ID[j] - 1].onSide && nod[elem.ID[nextJ] - 1].onSide && (j % 2))
                {
                    elem.boki[j].pktCalk[k].x = -1 + (2 * (j % 3));
                    elem.boki[j].pktCalk[k].y = pCal[k * punktyCalk].y;
                }

                vector<double> N = liczN(elem.boki[j].pktCalk[k].x, elem.boki[j].pktCalk[k].y);
                for (int m = 0; m < 4; m++) {
                    elem.boki[j].N[m][k] = N[m];
                }

            }
        }
    }
}

GlobalData::GlobalData() {
    odczytaj();
}

void GlobalData::odczytaj()
{
    ifstream file(data_path);
    string line;
    string temp;
    int value;
    int index;
    int licznik = 0;

    string nazwy[10] = { "SimulationTime", "SimulationStepTime", "Conductivity", "Alfa", "Tot", "InitialTemp", "Density", "SpecificHeat", "liczba węzłów", "liczba elementów" };
    while (getline(file, line) && licznik < 10) {
        index = line.rfind(" ") + 1;
        temp = line.substr(index, 16);
        value = stoi(temp);

        switch (licznik) {
        case 0: SimulationTime = value; break;
        case 1: SimulationStepTime = value; break;
        case 2: Conductivity = value; break;
        case 3: Alfa = value; break;
        case 4: Tot = value; break;
        case 5: InitialTemp = value; break;
        case 6: Density = value; break;
        case 7: SpecificHeat = value; break;
        case 8: nN = value; break;
        case 9: nE = value; break;
        }


        cout << nazwy[licznik] << ": " << value << endl;
        licznik++;
    }

}