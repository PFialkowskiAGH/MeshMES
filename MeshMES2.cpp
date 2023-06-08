#include <vector>
#include <iostream>

using namespace std;

vector<double> w = {1, 1};
vector<double> e = { -0.5773502692, 0.5773502692 };
vector<double> n1 = { 0.5 * (1 - e[0]), 0.5 * (1 - e[1]) };
vector<double> n2 = { 0.5 * (1 + e[0]), 0.5 * (1 + e[1]) };
double dTau = 10;
//int time = 2;

vector<double> GaussElimination(const vector<vector<double>>& _A, const vector<double>& _B) {
    int n = _A.size();
    vector<double> x(n);

    vector<vector<double>> tmpA(n, vector<double>(n + 1));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmpA[i][j] = _A[i][j];
        }
        tmpA[i][n] = _B[i];
    }

    double tmp;

    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            tmp = tmpA[i][k] / tmpA[k][k];
            for (int j = k; j < n + 1; j++) {
                tmpA[i][j] -= tmp * tmpA[k][j];
            }
        }
    }

    for (int k = n - 1; k >= 0; k--) {
        tmp = 0;
        for (int j = k + 1; j < n; j++) {
            tmp += tmpA[k][j] * x[j];
        }
        x[k] = (tmpA[k][n] - tmp) / tmpA[k][k];
    }

    return x;
}

struct GlobalData
{
    //0, 0.08, 300, 100, 1200, 1000, 700, 25, 7800, 9, 8
    double rMin; //promieñ minimalny wsadu
    double rMax; //promieñ maksymalny wsadu
    double alfaAir; //wspó³czynnik konwekcyjnej wymiany ciep³a
    double tempBegin; //temperatura pocz¹tkowa
    double tempAir; //temperatura otoczenia
    double tauMax; //czas procesu s
    double c; //efektywne ciep³o w³aœciwe
    double k; //wspó³czynnik przewodzenia ciep³a
    double rho; //gêstoœæ metalu
    int nN; //liczba wezlow siatki MES
    int nE; //liczba elementow siatki MES
    double dR;

    GlobalData(double rMin, double rMax, double alfaAir, double tempBegin, double tempAir, double tauMax, double c, double k, double rho, int nN, int nE)
    {
        this->rMin = rMin;
        this->rMax = rMax;
        this->alfaAir = alfaAir;
        this->tempBegin = tempBegin;
        this->tempAir = tempAir;
        this->tauMax = tauMax;
        this->c = c;
        this->k = k;
        this->rho = rho;
        this->nN = nN;
        this->nE = nE;
        this->dR = (rMax-rMin) / nE;
    }
};

struct Node
{
    double t; //temperatura na wezle (wynikowa)
    double x; //obliczyc L zamiast x - dla elementu
    int bc; //warunki brzegowe 1 - konwekcja, 2 - staly strumien ciepla

    Node(double t, double x, int bc)
    {
        this->t = t;
        this->x = x;
        this->bc = bc;
    }
};

struct Element
{
    double k; //wspolczynnik przewodzenia ciepla
    int id[2]{ 0,0 }; // E1 - N1, N2 - id[1, 2]
    double le; //dlugosc elementu
    vector<vector<double>> h; //lokalne macierz
    vector<double> p; //lokalne wektor

    Element(double k, int id[2], double le)
    {
        this->k = k;
        this->id[0] = id[0];
        this->id[1] = id[1];
        this->le = le;
        h = vector<vector<double>>(3, vector<double>(3));
        p = vector<double>(3);
    }

    void calculate(GlobalData globalData, vector<Node> nodes)
    {
        for (int i = 0; i < w.size(); ++i)
        {
            double intePos = (n1[i] * nodes[id[0]].x) + (n2[i] * nodes[id[1]].x);
            double intelTemperature = (n1[i] * nodes[id[0]].t) + (n2[i] * nodes[id[1]].t);
            double pom1 = (globalData.k * intePos * w[i] * (1 / globalData.dR));
            double pom2 = (globalData.c * globalData.rho * globalData.dR * intePos * w[i]) / dTau;

            //macierz H
            h[0][0] += (pom1 + (pom2 * n1[i] * n1[i]));
            h[0][1] += (-pom1 + (pom2 * n1[i] * n2[i]));
            h[1][0] += (-pom1 + (pom2 * n1[i] * n2[i]));
            h[1][1] += (pom1 + (pom2 * n2[i] * n2[i]));

            //wektor P
            p[0] += -((globalData.c * globalData.rho * globalData.dR * intelTemperature * intePos * w[i] * n1[i]) / dTau);
            p[1] += -((globalData.c * globalData.rho * globalData.dR * intelTemperature * intePos * w[i] * n2[i]) / dTau);

            //konwekcja, dodajmey raz !
            if (nodes[id[0]].bc == 1 && i == 1) {
                h[0][0] += 2 * globalData.alfaAir * globalData.rMax;
                p[0] += 2 * globalData.alfaAir * globalData.rMax * globalData.tempAir;
            }

            if (nodes[id[1]].bc == 1 && i == 1) {
                h[1][1] += 2 * globalData.alfaAir * globalData.rMax;
                p[1] += -(2 * globalData.alfaAir * globalData.rMax * globalData.tempAir);
            }

            cout << "H" << id[0] << "," << i << ":" << endl;
            cout << h[0][0] << " " << h[0][1] << endl;
            cout << h[1][0] << " " << h[1][1] << endl;
            cout << endl;

            cout << "P" << id[0] << "," << i << ":" << endl;
            cout << p[0] << " " << p[1] << endl;
            cout << endl;
        }
    }
};

struct siatkaMES
{
    vector<Element> elementVector; //E = new element[nE]
    vector<Node> nodeVector;

    siatkaMES(GlobalData globalData, int side) {

        //Tworzenie nodów
        int bc;
        for (int i = 0; i < globalData.nN; i++) {
            if (i == globalData.nN - 1) bc = 1;
            else bc = 0;
            nodeVector.push_back(Node(globalData.tempBegin, globalData.dR * i, bc));
        }

        //Tworzenie elementów siatki
        for (int i = 0; i < globalData.nE; i++)
        {
            int id[2] = { i, i + 1 };
            double le = nodeVector[i + 1].x - nodeVector[i].x;

            elementVector.push_back(Element(globalData.k, id, globalData.dR));
        }
    }

    void calculate(GlobalData globalData)
    {
        //for (int j = 0; j < time; j++)
        //{
            for (int i = 0; i < elementVector.size(); i++)
            {
                elementVector[i].calculate(globalData, nodeVector);
            }
        //}
    }
};

struct SOE //System of Equations - uklad rownan
{
    vector<vector<double>> hg; //nN x nN
    vector<double> pg; //nN x 1
    vector<double> t; //nN x 1

    void calculate(siatkaMES siatka)
    {
        int siatkaMeshSize = siatka.nodeVector.size();
        hg = vector<vector<double>>(siatkaMeshSize, vector<double>(siatkaMeshSize));
        pg = vector<double>(siatka.nodeVector.size());
        t = vector<double>(siatka.nodeVector.size());

        for (int i = 0; i < siatka.elementVector.size(); i++)
        {
            hg[siatka.elementVector[i].id[0]][siatka.elementVector[i].id[0]] += siatka.elementVector[i].h[0][0];
            hg[siatka.elementVector[i].id[0]][siatka.elementVector[i].id[1]] += siatka.elementVector[i].h[0][1];
            hg[siatka.elementVector[i].id[1]][siatka.elementVector[i].id[0]] += siatka.elementVector[i].h[1][0];
            hg[siatka.elementVector[i].id[1]][siatka.elementVector[i].id[1]] += siatka.elementVector[i].h[1][1];

            pg[siatka.elementVector[i].id[0]] += siatka.elementVector[i].p[0];
            pg[siatka.elementVector[i].id[1]] += siatka.elementVector[i].p[1];
        }

        t = GaussElimination(hg, pg);
        cout << endl << "Temperatura w danych Node: " << endl;
        for (int i = 0; i < siatkaMeshSize; i++) cout << t[i] << endl;
    }
};

int main() {
    GlobalData globalData = GlobalData(0, 0.08, 300, 100, 1200, 1000, 700, 25, 7800, 9, 8);

    siatkaMES siatka = siatkaMES(globalData, 1);
    SOE soe = SOE();

    siatka.calculate(globalData);
    soe.calculate(siatka);

    return 0;
};
