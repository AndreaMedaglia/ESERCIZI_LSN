/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 8
			 ESERCIZIO 8.1

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "main.h"

using namespace std;
 
int main (int argc, char *argv[]){

input(); //inizialization 
ofstream out ("Output/Minimization.txt");
ofstream out2 ("Output/Energy.txt");
ofstream out3 ("Output/Isto.txt");

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
	 input >> property;
	 if( property == "RANDOMSEED" ){
	    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	    rnd.SetRandom(seed,p1,p2);
	 }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	 
double ave[M], ave2[M], sum_prog[M], sum2_prog[M], err_prog[M], e[N];

Var.resize(q);

for(int l=0; l<q; l++) {		// ciclo grande sulla coppia di parametri mu sigma
	x=x_0;
	Var[l].energy.resize(M);
	Var[l].error.resize(M);
	Var[l].isto.resize(N);
	conta=0;
	for(int i=0; i<M; i++) {	// inizializzazione vettori	
		ave[i]=ave2[i]=sum_prog[i]=sum2_prog[i]=err_prog[i]=0.0;
	}
	for(int i=0; i<N; i++) {	// inizializzazione vettori
		e[i]=0.0;
	}

	Var[l].par_sigma=rnd.Rannyu(sigma_min, sigma_max);
	Var[l].par_mu=rnd.Rannyu(mu_min, mu_max);

	for (int i=0; i<N; i++) {	// Ciclo relativo all'algoritmo di Metropolis

		x_new=x;

		b=rnd.Rannyu();
		g=rnd.Rannyu();
		d=rnd.Rannyu(_min,_max); 

		if(b<0.5000) {x_new+=-d;}
		else {x_new+=d;}

		alpha = fmin ( 1.0 , (distribution(x_new, Var[l].par_sigma, Var[l].par_mu)/distribution(x, Var[l].par_sigma, Var[l].par_mu)) );

		if (g<=alpha) { 
			x=x_new;
			conta ++;
		}
		Var[l].isto[i]=x;
		e[i]=pot(x)+kin(x, Var[l].par_sigma, Var[l].par_mu);	
	}

	for(int i=0; i<M; i++) {	// Ciclo per il blocking
		double sum=0;
		for (int j=0; j<L; j++) {
			int k=j+i*L;
			sum += e[k];
		}
			ave[i] = sum/L;
			ave2[i] = ave[i] * ave[i]; 
	}
		
	for(int i=0; i<M; i++) {	// Ciclo per il blocking
		for (int j=0; j<(i+1); j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] /= (i+1);
		Var[l].energy[i]=sum_prog[i];
		sum2_prog[i] /= (i+1);
		err_prog[i] = error(sum_prog[i], sum2_prog[i], i) ;
		Var[l].error[i]=err_prog[i];
	}

	Var[l].en_best=sum_prog[M-1];
	
}
	
sort(Var.begin(), Var.end());		// ordino il vettore di struct in base al valore di en_best

for(int l=0; l<q; l++) {		// stampo output
	out << Var[l].en_best << "   " << Var[l].par_sigma << "   " << Var[l].par_mu << endl;	
}

for(int l=0; l<M; l++) {		// stampo output
	out2 << l+1 << " " << Var[0].energy[l] << "   " << Var[0].error[l] << endl;	
}

for(int l=0; l<N; l++) {		// stampo output
	out3 << Var[0].isto[l] << endl;	
}

out.close();
out2.close();
out3.close();

rnd.SaveSeed();
return 0;
}


void input(void){
	ifstream ReadInput;
	ReadInput.open("input.dat");
	ReadInput >> N;
	ReadInput >> M;
	L=N/M;
	ReadInput >> q;
	ReadInput >> x_0;
	x=x_0;
	ReadInput >> _min;
	ReadInput >> _max;
	ReadInput >> sigma_min;
	ReadInput >> sigma_max;
	ReadInput >> mu_min;
	ReadInput >> mu_max;
	ReadInput.close();
}


double error(double AV, double AV2, int n) {
	if(n==0) { return 0; }
	else { return sqrt((AV2-AV*AV)/n); }
}


double distribution(double _x, double _sigma, double _mu) {
	double num1=(_x-_mu)*(_x-_mu)/(2*_sigma*_sigma);
	double num2=(_x+_mu)*(_x+_mu)/(2*_sigma*_sigma);
	return pow((exp(-num1)+exp(-num2)),2); 
}


double pot(double _x) {
	return (pow(_x,4)-5.0*pow(_x,2)/2.0); 
}


double kin(double _x, double _sigma, double _mu) {
	double num1=(_x-_mu)*(_x-_mu)/(2*_sigma*_sigma);
	double num2=(_x+_mu)*(_x+_mu)/(2*_sigma*_sigma);
	double num3=(_x-_mu)/(_sigma*_sigma);
	double num4=(_x+_mu)/(_sigma*_sigma);
	return -0.5*(exp(-num1)*(num3*num3-1.0/(_sigma*_sigma))+exp(-num2)*(num4*num4-1.0/(_sigma*_sigma)));
}

/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 8
			 ESERCIZIO 8.1

*****************************************************************/
