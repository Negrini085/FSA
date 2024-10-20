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
//              4° riga     --> starting point 


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
        cout << "Starting point catena" << endl; 
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
        cout << "Starting point: " << contenitore[3] << endl;
        cout << endl;
    }
}

// Algoritmo per fare sampling del cammino di singola particella che vogliamo descrivere
// Levy harmonic algorithm, in questo caso consideriamo starting ed ending point come coincidenti
void mossaCammino(double dt, int Nid, TRandom* generatore, vector<double> &config){
    
    //Scelgo punto d'inizio per il cammino
    double start = config[1];
    double end = config[1];
    config[0] = start;

    //Calcolo elementi necessari per gaussiana
    double y1, y2;
    double mean, sigma;

    for(int i=1; i<int(size(config)); i++){
        y1 = 1/tanh(dt) + 1/tanh((Nid - i) * dt);
        y2 = start/sinh(dt) + end/sinh((Nid - i) * dt);

        mean = y2/y1; sigma = 1/sqrt(y1);
        config[i] = mean + generatore -> Gaus(0, sigma);

        start = config[i];
    }
}

// Metodo per la formazione dell'istogramma
void istoPosizioni(vector<int> &histo, double pos) {

    double appo = pos + 4;
    int index = int(appo * 10);
    
    if(index <= 79){ histo[index] += 1; }
}

// Metodo per stampare l'istogramma
void stampaHisto(string nome, const vector<int> &histo, int Niter) {

    ofstream fileout;
    fileout.open(nome);

    double appo = 0;
    for(int i = 0; i<int(size(histo)); i++) {
        appo = double(histo[i])/Niter;
        fileout << appo << endl;
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
    double beta = paramIn[0];
    int Ncompl = int(paramIn[1]);
    int Niter = int(paramIn[2]);
    vector<double> config(Ncompl, 0); config[0] = paramIn[3];

    // Contenitore per istogramma e studio accettazione delle mosse
    vector<int> histo(80, 0);

    double dt = beta/Ncompl;    //Intervallo di tempo immaginario
    for(int i=0; i<Niter; i++){
        mossaCammino(dt, Ncompl, generator, config);
        istoPosizioni(histo, config[0]);
    }

    cout << endl;
    cout << "Stampa dell'istogramma a file: histo.dat" << endl;
    stampaHisto("histo.dat", histo, Niter);

    delete generator;
    return 0;
}