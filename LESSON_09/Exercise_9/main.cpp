/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 9
			 ESERCIZIO 9

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
else {inizialization_circle();}         	// N_in=1

ofstream Best_l, Ave_l;
Best_l.open("Output/Best_l.dat");
Ave_l.open("Output/Average_l.dat");

pop_old.resize(N_pop);
pop_new.resize(N_pop);

for(int i=0; i<N_pop; i++) {			// setto la popolazione iniziale: N_pop cromosomi [1,...,N_cities]
	conf_0(pop_old[i].path);   		// in base ai percorsi iniziali calcolo le rispettive lunghezze L^2
	pop_old[i].length=chromo_length(pop_old[i].path);
}

for(int i=0; i<N_iter; i++) {

	sort(pop_old.begin(), pop_old.end());

	for(int j=0; j<N_pop; j++) {		// riempio la popolazione nuova selezionando elementi tramite "selection" dalla vecchia
		pop_new[j]=pop_old[selection(pop_old)];
		pop_new[j].length=chromo_length(pop_new[j].path);
	}
	
	sort(pop_new.begin(), pop_new.end());
	
	for(int j=0; j<N_pop; j++) {		// su una frazione della popolazione nuova faccio crossover
		cross=rnd->Rannyu();			
		if (cross<p_crossover) {
			crossover(pop_new[selection(pop_new)].path, pop_new[selection(pop_new)].path);
		}
	}

	for(int j=0; j<N_pop; j++) {		// su una frazione della popolazione nuova faccio permutazioni
		mut=rnd->Rannyu();
		if (mut<p_permutation) {
			permutation(pop_new[j].path);
		}
	}

	for(int j=0; j<N_pop; j++) {		// su una frazione della popolazione nuova faccio shift
		sh=rnd->Rannyu();
		if (sh<p_shift) {
			shift(pop_new[j].path);	
		}
	}

	for(int j=0; j<N_pop; j++) {		// su una frazione della popolazione nuova faccio una inversione
		inv=rnd->Rannyu();
		if (sh<p_inversion) {
			inversion(pop_new[j].path);	
		}
	}

	for(int j=0; j<N_pop; j++) {		// aggiorno le lunghezze dei percorsi e poi ordino
		pop_new[j].length=chromo_length(pop_new[j].path);
	}
	sort(pop_new.begin(), pop_new.end());

	for(int j=0; j<N_pop; j++) {		// popolazione vecchia = popolazione nuova (per ripartire nelle iterazioni)
		pop_old[j]=pop_new[j];
	}
	sort(pop_old.begin(), pop_old.end());
	
	sum=0;
	for(int j=0; j<(N_pop/2); j++) {	// media delle lunghezze sulla prima metà (la migliore) della popolazione
		sum+=pop_old[j].length;
	}

	Ave_l << i+1 << " " << sum/((double)N_pop/2) << endl;
	Best_l << i+1 << " " << pop_old[0].length << endl;
	//cout << pop_old[0].length << endl;	// scommentare per un controllo rapido da terminale della lunghezza minore di ogni generazione
}


	print_conf(pop_old[0].path);

Best_l.close();
Ave_l.close();
delete rnd;
return 0;
}


//////////   FUNCTIONS   //////////////

void input(void){
	ifstream ReadInput;
	ReadInput.open("input.dat");

	ReadInput >> N_cities;
	ReadInput >> N_pop;
	ReadInput >> N_iter;
	ReadInput >> N_mut;
	ReadInput >> N_in;
	ReadInput >> p_crossover;
	ReadInput >> p_permutation;
	ReadInput >> p_shift;
 	ReadInput >> p_inversion;

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

void mutation(vector<int> &v){         
	v.resize(N_cities);
	double prob1=rnd->Rannyu();
	if(prob1<p_permutation) {permutation(v);}
	double prob2=rnd->Rannyu();
	if(prob2<p_shift) {shift(v);}
	double prob3=rnd->Rannyu();
	if(prob3<p_inversion) {inversion(v);}
}

void permutation(vector<int> &v){         // permutazione, scambia N_mut coppie di città in modo casuale (N_mut=1 tendenzialmente)
	v.resize(N_cities);
	for (int i=0; i<N_mut; i++) {
		int a=floor(rnd->Rannyu(0,N_cities));
		int b=floor(rnd->Rannyu(0,N_cities));
		while (b==a) {
			b=floor(rnd->Rannyu(0,N_cities));
		}
		double v_ = v[a];
		v[a] = v[b];
		v[b] = v_;
	}
	check(v);
}

void shift(vector<int> &v){		// shifta di una quantità intera a<N_cities la lista di città
	v.resize(N_cities);
	int a=floor(rnd->Rannyu(0,N_cities)); 
	rotate(v.rbegin(), v.rbegin()+a, v.rend());
	check(v);
}

void inversion(vector<int> &v){     	// effettua una "inversione a specchio" sulla lista di città    
	v.resize(N_cities);
	vector<int> v_app (N_cities);
	for(int i=0; i<N_cities; i++) {
		v_app[i]=v[i];
	}
	for(int i=0; i<N_cities; i++) {
		v[i]=v_app[fabs(N_cities-(i+1))];
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

int selection(vector<Chromosome> &v) {        // Da usare con il vettore dei cromosomi GIA' ORDINATO!!! :
	v.resize(N_pop);		      // seleziona (con una probabilità proporzionale ad una potenza dell'inverso della lunghezza) 
	vector<double> w (N_pop);	      // l'indice di uno dei cromosomi della popolazione
	double norm=0.0;		       
	for (int i=0; i<N_pop; i++) {
		w[i]=1/(v[i].length*v[i].length*v[i].length*v[i].length);  // peso 1/l^4 
		norm+=w[i];
	}
	for (int i=0; i<N_pop; i++) {
		w[i] /= norm;
	}
	double p = rnd->Rannyu();
	double l=0;
	for (int i=0; i<N_pop; i++) {
		l+=w[i];
		if (p<l){
			selected=i;
			break;
		}
	}
return selected;
}

void crossover(vector<int> &v, vector<int> &w){  // Crossover, estraggo un intero casuale m, l'indice dal quale "taglio" i cromosomi
	v.resize(N_cities);			 // alla fine faccio un controllo: i cromosomi hanno tutte le città una volta sola ?
	w.resize(N_cities);
	vector<int> w_app (N_cities);
	int m=floor(rnd->Rannyu(1,(N_cities))); 
	for (int i=0; i<N_cities; i++) {
		w_app[i]=w[i];
	}
	int l=0;
	int h=0;
	for (int i=0; i<N_cities; i++) {
		int c=0;
		for (int j=0; j<m; j++) {	// Se trovo un numero uguale metto c=1
			if(v[i]==w[j]) {
				c=1;
			}	
		}
		if (c==0) {			// se non ho trovato numeri uguali allora c=0, quindi devo inserire v[i] in w[m+...]
			w[m+l]=v[i];
			l++;
		}			
	} 
	for (int i=0; i<N_cities; i++) {
		int d=0;
		for (int j=0; j<m; j++) {   
			if(w_app[i]==v[j]) {
				d=1;
			}	
		}
		if (d==0) {		  
			v[m+h]=w_app[i];
			h++;
		}			
	} 
	check(v);
	check(w); 
}

void print_conf(vector<int> &v) {		// Stampa la configurazione del singolo cromosoma sul terminale 
	v.resize(N_cities);			// in caso aggiungere opzione per stampare su file output

	ofstream Conf; 
	Conf.open("Output/Configuration.dat",ios::app);		// rinominare questi output (best_path)

	for(int i=0; i<N_cities; i++) {
		Conf << x_cities[v[i]-1] << " " << y_cities[v[i]-1] << endl;
	}
	Conf << endl;
	Conf.close();
}

/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 9
			 ESERCIZIO 9

*****************************************************************/

