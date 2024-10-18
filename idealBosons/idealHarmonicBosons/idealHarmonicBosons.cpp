#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include <TRandom3.h>

using namespace std;


//  ATTENZIONE: il file param.dat contiene i parametri da fornire per la simulazione 
//              1° riga     --> beta
//              2° riga     --> numero di particelle


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
        cout << "1° riga    --> Beta minimo" << endl; 
        cout << "2° riga    --> Numero particelle" << endl; 
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

// Funzione per valutare le derivate dei pesi delle permutazioni
void derWeightCalc(const double &beta, const int &Npart, const vector<double> weight, vector<double> &derWeight){
    for(int i=1; i<=Npart; i++){
        derWeight.push_back(-3 * weight[i-1] * i * exp(-i * beta)/(1-exp(-i * beta)));
    }
}

// Funzione per il calcolo ricorsivo della funzione di partizione
vector<double> canonicRecursion(const vector<double> &weight){
    vector<double> fPart; double appo = 0;
    fPart.push_back(1);
    
    for(int m=1; m<=int(size(weight)); m++) {
        for(int i=0; i<m; i++){
            appo += fPart[i] * weight[m-i-1];
        } 
        fPart.push_back(appo/m); appo = 0;
    } 

    return fPart;
}


/******************************************************
*    Campionamento delle permutazioni del sistema     *
******************************************************/

// Funzione per fare tower sampling
int towerSample(const vector<double> &prob, TRandom3* generatore) {
    vector<double> cumulativa;
    cumulativa.push_back(0);
    for(int i=0; i<int(size(prob)); i++){
        cumulativa.push_back(cumulativa[i] + prob[i]);
    }

    int ind = 0;
    double val = generatore -> Uniform(0, cumulativa[int(size(cumulativa))-1]);
    for(int i=0; i<int(size(cumulativa))-1; i++){
        if(cumulativa[i] <= val and val <= cumulativa[i+1]){ 
            ind = i; 
            return ind;
        }
    }

    cout << "Errore in calcolo della permutazione da fare" << endl;
    cout << "Termino esecuzione del programma!" << endl;
    exit(-1);
}

// Funzione per la determinazione diretta dei cicli considerati
vector<int> cicliDiretti(const vector<double> &weight, const vector<double> &fPart, TRandom3* generatore){
    vector<int> lCicli(int(size(weight)), 0);
    int M = int(size(weight)); int ind = 1;

    vector<double> prob;
    while(M>0){

        // Creo il vettore delle probabilità & trovo indice per tower sampling
        for(int i=0; i<M; i++){ prob.push_back(weight[i] * fPart[M-1-i]); }
        ind = towerSample(prob, generatore);
        
        // Aggiorno lunghezza ciclo e vario lunghezza permutazione
        M = M - (ind+1);
        lCicli[ind] += 1;
        prob.clear();
    }

    return lCicli;
}


/******************************************************
*   Campionamento delle configurazioni del sistema    *
******************************************************/

// Algoritmo per fare sampling del cammino di singola particella che vogliamo descrivere
// Levy harmonic algorithm, in questo caso consideriamo starting ed ending point come coincidenti
void mossaCammino(double dt, int kPart, TRandom* generatore, double first, vector<double> &config){
    
    //Scelgo punto d'inizio per il cammino
    double start = first;
    double end = start;
    config[0] = first;

    //Calcolo elementi necessari per gaussiana
    double y1, y2;
    double mean, sigma;

    for(int i=1; i<kPart; i++){
        y1 = 1/tanh(dt) + 1/tanh((kPart - i) * dt);
        y2 = start/sinh(dt) + end/sinh((kPart - i) * dt);

        mean = y2/y1; sigma = 1/sqrt(y1);
        config[i] = mean + generatore -> Gaus(0, sigma);

        start = config[i];
    }

    config.push_back(end);

}

// Bosoni armonici diretti
vector<double> bosoniDiretti(double beta, const vector<double> &weight, const vector<double> &fPart, TRandom3* generatore){
    double first = 0;
    vector<double> appo;
    vector<double> cordx;
    vector<double> cordy;
    vector<double> cordz;
    vector<double> conftot;

    // Determino i cicli con i quali stiamo lavorando
    vector<int> lCicli = cicliDiretti(weight, fPart, generatore);

    // Riempio il vettore delle configurazioni
    for(int i=0; i<int(size(lCicli)); i++){
        // Considero solo quei cicli di permutazioni non vuote
        if(lCicli[i] != 0){

            // Lavoro tutte le volte che è necessario
            for(int j=0; j<lCicli[i]; j++){

                // Coordinata x
                first = generatore -> Gaus(0, 1/sqrt(beta));
                mossaCammino(beta, i, generatore, first, appo);

                cordx.insert(cordx.end(), appo.begin(), appo.end());
                appo.clear();


                // Coordinata y
                first = generatore -> Gaus(0, 1/sqrt(beta));
                mossaCammino(beta, i, generatore, first, appo);

                cordy.insert(cordy.end(), appo.begin(), appo.end());
                appo.clear();


                // Coordinata z
                first = generatore -> Gaus(0, 1/sqrt(beta));
                mossaCammino(beta, i, generatore, first, appo);

                cordz.insert(cordz.end(), appo.begin(), appo.end());
                appo.clear();
            }
        }
    }    
    cout << "daje" << endl;

    // Contenitore per le varie coordinate
    conftot.insert(conftot.end(), cordx.begin(), cordx.end());
    conftot.insert(conftot.end(), cordy.begin(), cordy.end());
    conftot.insert(conftot.end(), cordz.begin(), cordz.end());


    return conftot;
}


void stampaConf(string name, const vector<double> &conf){
    ofstream fileout;
    fileout.open(name);

    for(int i=0; i<int(size(conf)); i++) {
        fileout << conf[i] << endl;
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
    vector<double> conf;
    TRandom3* generator = new TRandom3();


    // Leggo i parametri iniziali
    parametriSimulativi(nome_par, paramIn);
    stampaPar(paramIn);

    // Parametri per la simulazione
    double beta = paramIn[0];
    int Npart = int(paramIn[1]);


    // Calcolo i pesi delle permutazioni e la loro derivata
    vector<double> weight;
    vector<double> fPart;
    
    cout << "Determinazione dei pesi" << endl;
    weightCalc(beta, Npart, weight);
    fPart = canonicRecursion(weight);

    // Studio le configurazioni del sistema
    cout << "Campionamento della configurazione" << endl;
    conf = bosoniDiretti(beta, weight, fPart, generator);
    stampaConf("bosConf.dat", conf);


    delete generator;
    return 0;
}