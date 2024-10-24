#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <TRandom3.h>

using namespace std;


//  ATTENZIONE: il file param.dat contiene i parametri da fornire per la simulazione 
//              1° riga     --> beta minimo
//              2° riga     --> beta massimo
//              3° riga     --> numero di particelle


/***********************************************
*        Lettura dei parametri di input        *
***********************************************/

//Funzione per leggere i parametri di input
void parametriSimulativi(string nome, vector<double> &contenitore){
    
    int conta = 0;
    double appo = 0; string pippo;

    ifstream filein;
    filein.open(nome);

    //Controllo corretta apertura del file
    if(!filein) {
        cerr << "Errore in apertura file che fornisce la configurazione di input." << endl;
        cerr << "Termino esecuzione del programma" << endl;
        exit(1);
    }

    while(!filein.eof()){
        getline(filein, pippo);
        conta++;
    }

    if(conta != 4){
        cout << "File dei parametri con input errati!" << endl;
        cout << "Formato richiesto: " << endl;
        cout << "1° riga    --> Beta minimo" << endl; 
        cout << "2° riga    --> Beta massimo" << endl; 
        cout << "3° riga    --> Numero particelle" << endl; 
        cout << "4° riga    --> Vuota" << endl; 
        exit(-1);
    }

    else{
        filein.clear();
        filein.seekg(0,ios::beg);

        while(filein >> appo){
            contenitore.push_back(appo); 
        }
    }

    filein.close();
}


// Stampo i parametri della simulazione
void stampaPar(const vector<double> &contenitore){
    if(size(contenitore) == 0){
        cout << "Contenitore dei parametri vuoto, non abbiamo alcun parametro simulativo!" << endl;
    }

    else{
        cout << endl;
        cout << "Parametri per la simulazione" << endl;
        cout << "Beta minimo: " << contenitore[0] << endl;
        cout << "Beta massimo: " << contenitore[1] << endl;
        cout << "Numero di particelle: " << contenitore[2] << endl;
        cout << endl;
    }
}


/******************************************************
*        Studio delle permutazioni del sistema        *
******************************************************/

// Funzione per valutare i pesi delle permutazioni
void weightCalc(const double &beta, const int &Npart, vector<double> &weight){
    for(int i=1; i<=Npart; i++){
        weight.push_back(1/pow(1-exp(-i * beta), 3));
    }
}

// Funzione per valutare le derivate dei pesi delle permutazioni
void derWeightCalc(const double &beta, const int &Npart, const vector<double> weight, vector<double> &derWeight){
    for(int i=1; i<=Npart; i++){
        derWeight.push_back(-3 * weight[i-1] * i * exp(-i * beta)/(1-exp(-i * beta)));
    }
}

// Funzione per il calcolo ricorsivo della funzione di partizione
double canonicRecursion(const vector<double> &weight){
    vector<double> fPart; double appo = 0;
    fPart.push_back(1);
    
    for(int m=1; m<=int(size(weight)); m++) {
        for(int i=0; i<m; i++){
            appo += fPart[i] * weight[m-i-1];
        } 
        fPart.push_back(appo/m); appo = 0;
    } 

    return fPart[int(size(weight))];
}

// Funzione per il calcolo ricorsivo della funzione di partizione, energia e 
// del numero di particelle nel condensato
vector<double> canonicRecursionObs(const vector<double> &weight, const vector<double> &derWeight){
    vector<double> fPart; double appo = 0;
    vector<double> dfPart; double pippo = 0;
    fPart.push_back(1);
    dfPart.push_back(0);
    
    for(int m=1; m<=int(size(weight)); m++) {
        for(int i=0; i<m; i++){
            appo += fPart[i] * weight[m-i-1];
            pippo += fPart[i] * derWeight[m-i-1] + dfPart[i] * weight[m-i-1];
        } 
        fPart.push_back(appo/m); appo = 0;
        dfPart.push_back(pippo/m); pippo = 0;
    } 

//    cout << "La funzione di partizione è pari a: " << endl;
//    cout << "Z0: " << fPart[0] << endl;
//    cout << "Z1: " << fPart[1] << endl;
//    cout << "Z2: " << fPart[2] << endl;
//    cout << "Z3: " << fPart[3] << endl;
//    cout << "Z4: " << fPart[4] << endl;
//    cout << "Z5: " << fPart[5] << endl;
//
//    cout << endl << endl;
//    cout << "La derivata della funzione di partizione è pari a: " << endl;
//    cout << "Z0: " << dfPart[0] << endl;
//    cout << "Z1: " << dfPart[1] << endl;
//    cout << "Z2: " << dfPart[2] << endl;
//    cout << "Z3: " << dfPart[3] << endl;
//    cout << "Z4: " << dfPart[4] << endl;
//    cout << "Z5: " << dfPart[5] << endl;

    // Calcolo energia e frazione di condensato
    appo = -dfPart[int(size(weight))]/fPart[int(size(weight))];
    pippo = 0;

    for(int i=0; i<int(size(weight)); i++){
        pippo += fPart[i];
    }

    pippo = pippo/fPart[int(size(weight))];


    fPart.push_back(appo);
    fPart.push_back(pippo);

    return fPart;
}

// Funzione per stampare gli osservabili di nostro interesse a file
void stampaOss(const vector<double> usedBeta, const vector<double> energy, const vector<double> confrac, string name) {

    if(int(size(energy)) != int(size(confrac))) {
        cout << "Dimensioni differenti dei vettori delle energie e delle frazioni di condensato." << endl;
        cout << "Termine esecuzione programma" << endl;
        exit(-2);
    }

    ofstream fileout;
    fileout.open(name);

    int dim = int(size(energy));
    for(int i=0; i<dim; i++){
        fileout << usedBeta[i] << "    " << energy[i] << "    " << confrac[i] << endl;
    }

    fileout.close();
}


int main(int argc, char* argv[]){

    if(argc != 2) { 
        cout << "Utilizzo programma: " << argv[0] << "<file parametri>" << endl;
        return -1;
    }

    string nome_par = argv[1];

    vector<double> paramIn;
    TRandom3* generator = new TRandom3();


    // Leggo i parametri iniziali
    parametriSimulativi(nome_par, paramIn);
    stampaPar(paramIn);

    // Parametri per la simulazione
    double betamin = paramIn[0];
    double betamax = paramIn[1];
    int Npart = int(paramIn[2]);


    // Calcolo i pesi delle permutazioni e la loro derivata
    vector<double> weight;
    vector<double> derWeight;
    vector<double> appo; double beta;

    vector<double> energy;
    vector<double> confrac;
    vector<double> usedBeta;

    double dbeta = 0.01;
    int Nmax = int((betamax - betamin)/dbeta);

    // Cicliamo a varie temperature per vedere il comportamento del sistema 
    // al variare della temperatura
    for(int i = 0; i<=Nmax; i++) {
        beta = betamax - dbeta * i;
        weightCalc(beta, Npart, weight);
        derWeightCalc(beta, Npart, weight, derWeight);

        // Ricorsione per determinazione della temperatura
        appo = canonicRecursionObs(weight, derWeight);

        usedBeta.push_back(beta);
        energy.push_back(appo[Npart+1]);
        confrac.push_back(appo[Npart+2]/Npart);

        if(i%100 == 0 and i > 0){
            cout << "Fatte prime " << i << " temperature!" << endl;
        }

        weight.clear();
        derWeight.clear();

    }

    
    // Stampo a file i parametri d'interesse della nostra indagine, ossia
    // frazione di condensato ed energia del sistema bosonico
    stampaOss(usedBeta, energy, confrac, "osservabili.dat");

    delete generator;
    return 0;
}