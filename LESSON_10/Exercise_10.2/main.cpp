/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 10
			 ESERCIZIO 10.2

*****************************************************************/

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "main.h"


Random* rnd;

using namespace std;

int main (int argc, char *argv[]){

MPI::Init(argc,argv);
input();

if(N_in==0) {inizialization_square();}
else {inizialization_circle();} 

int size = MPI::COMM_WORLD.Get_size();	// numero di nodi
int rank = MPI::COMM_WORLD.Get_rank();	// numero del nodo considerato in quel momento

ofstream Best_l;
Best_l.open("Output/Best_l_"  + std::to_string(rank) + ".dat");

ofstream Data;
Data.open("Output/Parametri_"  + std::to_string(rank) + ".dat");

pop_old.resize(size);
pop_new.resize(size);

conf_0(pop_old[rank].path);		// non voglio partire dalle città ordinate da 1 a N, inserisco una mutazione iniziale
double h1=rnd->Rannyu();
double h2=rnd->Rannyu();
	if (h1<0.800) { permutation(pop_old[rank].path);}
	if (h2<0.200) { permutation(pop_old[rank].path);}
pop_old[rank].length=chromo_length(pop_old[rank].path); 		

// GENERAZIONE CASUALE DI BETA INIZIALE E BETA FINALE
if (rank==0) {
beta=1;
beta_fin=300;
}
else if (rank==1) {
beta=floor(rnd->Rannyu(1,10));
beta_fin=floor(rnd->Rannyu(230,300));
}
else if (rank==2) {
beta=floor(rnd->Rannyu(1,5));
beta_fin=floor(rnd->Rannyu(250,280));
}
else {
beta=floor(rnd->Rannyu(5,10));
beta_fin=floor(rnd->Rannyu(200,300));
}
// GENERAZIONE CASUALE DELL'INCREMENTO SU BETA
if (rank==0) { }
else if (rank==1) {incr=rnd->Rannyu(0.1,0.2); }
else if (rank==2) {incr=rnd->Rannyu(0.2,0.3);}
else {incr=rnd->Rannyu(0.3,0.4);}
double beta_0=beta;
N_step=(int)(beta_fin-beta)/incr;

for (int i=0; i<N_step; i++) {					// numero di cicli totali, dipende da beta iniziale, finale e incremento	
	for (int j=0; j<N_metro; j++) {				// numero di passi metropolis per ciclo i-esimo
		pop_old[rank].length=chromo_length(pop_old[rank].path);
		pop_new[rank]=pop_old[rank];
		double r=rnd->Rannyu();
		if(r<0.500) { permutation(pop_new[rank].path); }// al 50% effettuo permutazioni
		else {shift(pop_new[rank].path); }		// al 50% effettuo shift
		pop_new[rank].length=chromo_length(pop_new[rank].path);
		alpha = fmin( 1.0 , exp(-beta*(pop_new[rank].length-pop_old[rank].length)) );
		p = rnd->Rannyu();
		if (p<=alpha) {					// Metropolis: accetto il passo o no? 
			pop_old[rank]=pop_new[rank];
			}	
	}
	
	Best_l <<  i+1 <<  "  " << pop_old[rank].length << endl;
	beta+=incr;
}

Data << rank << "	" << beta_0 << "	" << beta_fin << "	" << incr << "	" << pop_old[rank].length << endl; // stampo su output

cout << "Nodo " << rank << " ha finito." << endl;		// controllo su terminale dei nodi che hanno finito

print_conf(pop_old[rank].path, rank); 
Best_l.close();
Data.close();

MPI::Finalize();
delete rnd;
return 0;
}


//////////   FUNCTIONS   //////////////

void input(void){
	ifstream ReadInput;
	ReadInput.open("input.dat");

	ReadInput >> N_cities;
	ReadInput >> N_metro;
	ReadInput >> incr;
	ReadInput >> N_in;

	ReadInput.close();
}

void inizialization_square(void){   	// Assegna una posizione casuale alle città nel quadrato [0,1]x[0,1]
	rnd = new Random();
	x_cities.resize(N_cities);
	y_cities.resize(N_cities);
	for (int i=0; i<N_cities; i++) {
		x_cities[i]=rnd->Rannyu();
		y_cities[i]=rnd->Rannyu();
	}	
}

void inizialization_circle(void){   	// Assegna una posizione casuale alle città sulla circonferenza unitaria  
	rnd = new Random();
	x_cities.resize(N_cities);
	y_cities.resize(N_cities);
	for (int i=0; i<N_cities; i++) {
		theta=rnd->Rannyu(0,2*M_PI);
		x_cities[i]=cos(theta);
		y_cities[i]=sin(theta);
	}		
}

void conf_0(vector<int> &v) {          // Assegna al cromosoma la configurazione iniziale [1, 2, ...., n]
	v.resize(N_cities);
	for(int i=0; i<N_cities; i++) {
		v[i]=i+1;
	}
	check(v);
}




double chromo_length(vector<int> &v){  // restituisce la lunghezza del percorso (double), data la configurazione
	v.resize(N_cities);	       // da usare per assegnare il valore alla variabile length della struct 
	double l=0;
	double x_length;
	double y_length;

	for (int i=0; i<N_cities; i++) {
		if(i==(N_cities-1)) {
			x_length = x_cities[v[0]-1] - x_cities[v[i]-1];
			y_length = y_cities[v[0]-1] - y_cities[v[i]-1];
		}
		else {
			x_length = x_cities[v[i+1]-1] - x_cities[v[i]-1];
			y_length = y_cities[v[i+1]-1] - y_cities[v[i]-1];
		}
		l += x_length*x_length + y_length*y_length;
	}
return l;
}

void check(vector<int> &v){            // controlla che nel singolo cromosoma non siano presenti due città uguali  
	v.resize(N_cities);
	for (int i=0; i<N_cities; i++) {
		for (int j=0; j<N_cities; j++) {
			if(j!=i) {
				if(v[j]==v[i]) {cout << "Error! Chromosome has the same city twice" << endl;}
			}
		}
	}
}

void permutation(vector<int> &v){         // permutazione, scambia una coppia di città in modo casuale 
	v.resize(N_cities);

	int a=floor(rnd->Rannyu(0,N_cities));
	int b=floor(rnd->Rannyu(0,N_cities));
	while (b==a) {
		b=floor(rnd->Rannyu(0,N_cities));
	}
	double v_ = v[a];
	v[a] = v[b];
	v[b] = v_;

	check(v);
}

void shift(vector<int> &v){         
	v.resize(N_cities);
	int a=floor(rnd->Rannyu(0,N_cities)); 
	rotate(v.rbegin(), v.rbegin()+a, v.rend());
	check(v);
}


void print_conf(vector<int> &v, int n) {     // Stampa la configurazione del singolo cromosoma sul terminale 
	v.resize(N_cities);	             // in caso aggiungere opzione per stampare su file output

	ofstream Conf; 
	Conf.open("Output/Configuration_" + std::to_string(n) + ".dat");

	for(int i=0; i<N_cities; i++) {
		Conf << x_cities[v[i]-1] << " " << y_cities[v[i]-1] << endl;
	}
	Conf << endl;
	Conf.close();
}

/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 10
			 ESERCIZIO 10.2

*****************************************************************/

	


