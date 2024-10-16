#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <TRandom3.h>

using namespace std;

//  ATTENZIONE: il file param.dat contiene i parametri da fornire per la simulazione 
//              1° riga     --> beta
//              2° riga     --> numero di completezze
//              3° riga     --> numero di iterazioni


// Leggo la configurazione iniziale
void inizializzaCammino(string nome, vector<double> &contenitore){
    
    ifstream filein;
    double appo = 0;
    filein.open(nome);

    //Controllo corretta apertura del file
    if(!filein) {
        cerr << "Errore in apertura file che fornisce la configurazione di input." << endl;
        cerr << "Termino esecuzione del programma" << endl;
        exit(1);
    }

    while(!filein.eof()){
        filein >> appo;
        contenitore.push_back(appo);
    }

    filein.close();
}

// Leggo parametri simulativi
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

    filein.clear();
    filein.seekg(0,ios::beg);

    if(conta != 4){
        cout << "File dei parametri con input errati!" << endl;
        cout << "Formato richiesto: " << endl;
        cout << "Beta" << endl; 
        cout << "Numero completezze" << endl; 
        cout << "Numero iterazioni" << endl; 
        exit(-1);
    }

    else{
        while(filein >> appo){
            contenitore.push_back(appo); 
        }
    }

    filein.close();
}

// Stampo la configurazione desiderata
void stampaConfig(const vector<double> &contenitore){
    if(size(contenitore) == 0){
        cout << "Vettore vuoto: non contiene alcun elemento!" << endl;
    }

    else{
        cout << endl;
        cout << "Elementi del vettore: " << endl;
        for(int i=0; i<int(size(contenitore)); i++) {
            cout << contenitore[i] << endl;
        }
        cout << endl;
    }
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
        cout << "Numero di completezze: " << int(contenitore[1]) << endl;
        cout << "Numero di iterazioni: " << int(contenitore[2]) << endl;
        cout << endl;
    }
}

int main(int argc, char* argv[]){

    if(argc != 3) { 
        cout << "Utilizzo programma: " << argv[0] << " <file configurazione iniziale> <file parametri>" << endl;
        return -1;
    }

    string nome_confIn = argv[1];
    string nome_par = argv[2];

    vector<double> contieni;
    vector<double> paramIn;
    TRandom3* generator = new TRandom3();

    // Leggo la configurazione iniziale
    inizializzaCammino(nome_confIn, contieni);
    stampaConfig(contieni);

    // Leggo i parametri iniziali
    parametriSimulativi(nome_par, paramIn);
    stampaPar(paramIn);

    delete generator;
    return 0;
}