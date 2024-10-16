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
//              4° riga     --> delta (per intervallo mosse)
//              5° riga     --> starting point 


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

    if(conta != 6){
        cout << "File dei parametri con input errati!" << endl;
        cout << "Formato richiesto: " << endl;
        cout << "Beta" << endl; 
        cout << "Numero completezze" << endl; 
        cout << "Numero iterazioni" << endl; 
        cout << "Delta" << endl; 
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
        cout << "Delta (intervallo mosse): " << contenitore[3] << endl;
        cout << "Starting point: " << contenitore[4] << endl;
        cout << endl;
    }
}

// Funzione gaussiana
double gaussiana(double x, double x1, double sigma) {
    double coeff = 1.0 / (sigma * sqrt(2.0 * M_PI)); 
    double exponent = - (pow((x - x1), 2)) / (2 * pow(sigma, 2));  
    return coeff * exp(exponent);
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

    if(argc != 3) { 
        cout << "Utilizzo programma: " << argv[0] << " <file configurazione iniziale> <file parametri>" << endl;
        return -1;
    }

    string nome_confIn = argv[1];
    string nome_par = argv[2];

    vector<double> config;
    vector<double> paramIn;
    TRandom3* generator = new TRandom3();


    // Leggo i parametri iniziali
    parametriSimulativi(nome_par, paramIn);
    stampaPar(paramIn);

    // Parametri per la simulazione
    double beta = paramIn[0];
    int Ncompl = int(paramIn[1]);
    int Niter = int(paramIn[2]);
    double delta = paramIn[3];

    // Contenitore per istogramma e studio accettazione delle mosse
    int accRate = 0;
    vector<int> histo(80, 0);

    double dt = beta/Ncompl;    //Intervallo di tempo immaginario
    for(int i=0; i<Niter; i++){
        accRate += 1; 
        istoPosizioni(histo, config[0]);
    }

    cout << endl;
    cout << "Acceptance rate: " << endl;
    cout << double(accRate)/Niter << endl;

    cout << endl;
    cout << "Stampa dell'istogramma a file: histo.dat" << endl;
    stampaHisto("histo.dat", histo, Niter);

    delete generator;
    return 0;
}