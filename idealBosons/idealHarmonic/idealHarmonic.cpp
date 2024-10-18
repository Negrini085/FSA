#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <TRandom3.h>

using namespace std;



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

    if(conta != 3){
        cout << "File dei parametri con input errati!" << endl;
        cout << "Formato richiesto: " << endl;
        cout << "1° riga    --> Beta" << endl; 
        cout << "2° riga    --> Numero di particelle" << endl; 
        cout << "3° riga    --> Vuota" << endl; 
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
        cout << "Beta: " << contenitore[0] << endl;
        cout << "Numero di particelle: " << contenitore[1] << endl;
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
    double beta = paramIn[0];
    int Npart = int(paramIn[1]);


    // Calcolo i pesi delle permutazioni
    vector<double> weight;
    weightCalc(beta, Npart,weight);

    beta = canonicRecursion(weight);

    // Calcolo ricorsivamente la funzione di partizione
    cout << "La funzione di partizione ha valore: " << canonicRecursion(weight) << endl;

    delete generator;
    return 0;
}