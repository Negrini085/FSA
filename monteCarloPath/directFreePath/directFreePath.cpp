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
//              4° riga     --> dimensioni box


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

    if(conta != 5){
        cout << "File dei parametri con input errati!" << endl;
        cout << "Formato richiesto: " << endl;
        cout << "Beta" << endl; 
        cout << "Numero completezze" << endl; 
        cout << "Numero iterazioni" << endl; 
        cout << "Dimensioni box" << endl; 
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
        cout << "Numero di completezze: " << int(contenitore[1]) << endl;
        cout << "Numero di iterazioni: " << int(contenitore[2]) << endl;
        cout << "Dimensioni box: " << contenitore[3] << endl;
        cout << endl;
    }
}

// Funzione gaussiana
double gaussiana(double x, double x1, double sigma) {
    double coeff = 1.0 / (sigma * sqrt(2.0 * M_PI)); 
    double exponent = - (pow((x - x1), 2)) / (2 * pow(sigma, 2));  
    return coeff * exp(exponent);
}

// Metodo per costruire cammino di evoluzione libera, punto iniziale randomico fra 0 ed L
void mossaCammino(double dt, double L, int Nid, TRandom* generatore, vector<double> &config){
    
    // Scelgo punto iniziale per il cammino (estratto uniformemente in 0 -> L)
    double start = generatore -> Uniform(0, L);
    double end = start; double dt1 = 0;
    config[0] = start;

    double mean, sigma;

    for(int i=1; i<int(size(config)); i++) {
        dt1 = int((size(config) - i)) * dt;

        // Calcolo valore medio e varianza della gaussiana
        mean = (dt1 * config[i-1] + dt * end)/(dt + dt1);
        sigma = 1/sqrt(1/dt1 + 1/dt);

        // Campiono la mossa
        config[i] = mean + generatore -> Gaus(0, sigma);
    }
}

// Metodo per stampare la configurazione
void stampaConfig(ofstream &fileOut, const vector<double> &config) {

    for(int i = 0; i<int(size(config)); i++) {
        fileOut << config[i] << "    ";
    }
    fileOut << endl;

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
    int Ncompl = int(paramIn[1]);
    int Niter = int(paramIn[2]);
    double L = paramIn[3];
    vector<double> config(Ncompl, 0);

    // Contenitori & variabili per istogramma e studio dei cammini
    vector<int> histo(80, 0);
    ofstream fileOut; fileOut.open("camminiFree.dat");

    double dt = beta/Ncompl;    //Intervallo di tempo immaginario
    for(int i=0; i<Niter; i++){
        mossaCammino(dt, L, Ncompl, generator, config);
        stampaConfig(fileOut, config);
    }


    fileOut.close();
    delete generator;
    return 0;
}