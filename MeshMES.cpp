﻿#include <vector>
#include <iostream>

using namespace std;

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
    double tot; //temperatura otoczenia
    double alfa; //warunek brzegowy konwekcji
    double q; //gestosc strumienia ciepla
    double k; //wspolczynnik przewodzenia ciepla
    double l; //dlugosc preta
    double s; //pole przekroju preta
    int nN; //liczba wezlow siatki MES
    int nE; //liczba elementow siatki MES

    GlobalData(double tot, double alfa, double q, double k, double l, double s, int nN, int nE) 
    {
        this->tot = tot;
        this->alfa = alfa;
        this->q = q;
        this->k = k;
        this->l = l;
        this->s = s;
        this->nN = nN;
        this->nE = nE;
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
        //wyliczenie macierzy H
        if (nodes[id[0]].bc == 2) h[0][0] = (globalData.s * globalData.k * (1 / le)) + (globalData.alfa * globalData.s);
        else h[0][0] = (globalData.s * globalData.k * (1 / le));

        h[0][1] = globalData.s * globalData.k * (-1 / le);
        h[1][0] = globalData.s * globalData.k * (-1 / le);
        
        if (nodes[id[1]].bc == 2) h[1][1] = (globalData.s * globalData.k * (1 / le)) + globalData.alfa * globalData.s;
        else h[1][1] = globalData.s * globalData.k * (1 / le);

        //wyliczenie macierzy P
        if (nodes[id[0]].bc == 1) p[0] = globalData.q * globalData.s;

        if (nodes[id[1]].bc == 1) p[1] = globalData.q * globalData.s;

        if (nodes[id[0]].bc == 2) p[0] = -1 * globalData.alfa * globalData.tot * globalData.s;

        if (nodes[id[1]].bc == 2) p[1] = -1 * globalData.alfa * globalData.tot * globalData.s;
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
            if ((i == 0 && side == 1) || (i == globalData.nN - 1 && side == 2)) bc = 1;
            else if ((i == 0 && side == 2) || (i == globalData.nN - 1 && side == 1)) bc = 2;
            else bc = 0;

            nodeVector.push_back(Node(0, i * (globalData.l / (globalData.nN - 1)), bc));
        }

        //Tworzenie elementów siatki
        for (int i = 0; i < globalData.nN - 1; i++)
        {
            int id[2] = { i, i + 1 };
            double le = nodeVector[i + 1].x - nodeVector[i].x;

            elementVector.push_back(Element(globalData.k, id, le));
        }
    }

    void calculate(GlobalData globalData)
    {
        for (int i = 0; i < elementVector.size(); i++)
        {
            elementVector[i].calculate(globalData, nodeVector);
        }
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
    GlobalData globalData = GlobalData(400, 10, -150, 50, 5, 2, 3, 2);

    siatkaMES siatka = siatkaMES(globalData, 1);
    SOE soe = SOE();

    siatka.calculate(globalData);
    soe.calculate(siatka);

    return 0;
}
