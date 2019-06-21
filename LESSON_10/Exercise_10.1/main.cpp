/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 10
			 ESERCIZIO 10.1

*****************************************************************/

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

input();

if(N_in==0) {inizialization_square();}
else {inizialization_circle();}			//N_in=1         

ofstream Best_l;
Best_l.open("Output/Best_l.dat");

conf_0(pop_old.path);
pop_old.length=chromo_length(pop_old.path);   		

vector<double> acceptance (N_step, 0);

for (int i=0; i<N_step; i++) {					// numero di cicli totali, dipende da beta iniziale, finale e incremento	
	for (int j=0; j<N_metro; j++) {				// numero di passi metropolis per ciclo i-esimo
		pop_old.length=chromo_length(pop_old.path);
		pop_new=pop_old;
		double r=rnd->Rannyu();
		if(r<0.500) { permutation(pop_new.path); }	// al 50% effettuo permutazioni
		else {shift(pop_new.path); }			// al 50% effettuo shift
		pop_new.length=chromo_length(pop_new.path);
		alpha = fmin( 1.0 , exp(-beta*(pop_new.length-pop_old.length)) );
		p = rnd->Rannyu();
		if (p<=alpha) {					// Metropolis: accetto il passo o no?
			pop_old=pop_new;
			acceptance[i]++;
			}	
	}
	
	acceptance[i]/=N_metro;
	Best_l <<  i+1 <<  "  " << pop_old.length << endl;
	beta+=incr;
}


/*
for (int i=0; i<N_step; i++) {				// scommentare per un controllo su terminale dell'accettazione di Metropolis
	cout << acceptance[i]*100 << " % " << endl;	// voglio che sia intorno al 50%
}
*/

print_conf(pop_old.path); 
	
Best_l.close();
delete rnd;
return 0;
}


//////////   FUNCTIONS   //////////////

void input(void){
	ifstream ReadInput;
	ReadInput.open("input.dat");

	ReadInput >> N_cities;
	ReadInput >> N_metro;
	ReadInput >> beta;
	ReadInput >> beta_fin;
	ReadInput >> incr;
	ReadInput >> N_in;
	N_step=floor((beta_fin-beta)/incr)+1;

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

void shift(vector<int> &v){        	// shifta di una quantità intera a<N_cities la lista di città 
	v.resize(N_cities);
	int a=floor(rnd->Rannyu(0,N_cities)); 
	rotate(v.rbegin(), v.rbegin()+a, v.rend());
	check(v);
}


void print_conf(vector<int> &v) {            // Stampa la configurazione del singolo cromosoma sul terminale 
	v.resize(N_cities);	             // in caso aggiungere opzione per stampare su file output

	ofstream Conf; 
	Conf.open("Output/Configuration.dat");

	for(int i=0; i<N_cities; i++) {
		Conf << x_cities[v[i]-1] << " " << y_cities[v[i]-1] << endl;
	}
	Conf << endl;
	Conf.close();
}

/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 10
			 ESERCIZIO 10.1

*****************************************************************/

