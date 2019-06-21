/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 5
			 ESERCIZIO 5.1

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

double error(double AV, double AV2, int n);
double distribution(double _x, double _y, double _z);
double distribution2(double _x, double _y, double _z);

using namespace std;
 
int main (int argc, char *argv[]){

   int N=100000;  // NUMERO DI PASSI TOTALI DI METROPOLIS
   int M=100 ;    // NUMERO DI BLOCCHI
   int L=N/M ;    // ELEMENTI PER BLOCCO

   ofstream out1 ("Output/Coordinates.txt"); // output da rinominare in base allo stato che si sta campionando
   ofstream out2 ("Output/Distance.txt");


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



double x[3]={0.0, 0.0, 0.0};
double x_new[3], r[N], ave[M], ave2[M], sum_prog[M], sum2_prog[M], err_prog[M];
double a, b, g, alpha, dist;
double d;
int conta=0;
int conta2=0;
char h;
int n;


/////        Passi Metropolis per l'equilibrazione del sistema        /////
cout << "Do you want to make an equilibration of the system? [Y/N]" << endl;
cin >> h;
while (h != 'Y' && h != 'N') {
	cout << "Invalid input, try again with [Y/N] " << endl;
	cin >> h;
}
while (h == 'Y') { 
	dist = 0.0;
	cout << "For how many Metropolis steps you want to simulate? [integer]" << endl;
	cin >> n; 
	for (int i=0; i<n; i++) {

		for(int j=0; j<3; j++) {
			x_new[j]=x[j];
		}

		a=rnd.Rannyu();
		b=rnd.Rannyu();
		g=rnd.Rannyu();
		d=rnd.Gauss(3,1); 		// da usare nel caso multivariato per l'orbitale 2p per avere 50% di accettanza 
		//d=rnd.Gauss(1.3,0.5); 	// da usare nel caso multivariato per il Ground State per avere 50% di accettanza
		//d=rnd.Rannyu(2.5,3.5); 	// da usare nel caso uniforme per l'orbitale 2p per avere 50% di accettanza 
		//d=rnd.Rannyu(0.75,1.55); 	// da usare nel caso uniforme per il Ground State per avere 50% di accettanza

		if(a<0.3333) {			// asse x
			if(b<0.5000) {x_new[0]+=-d;}
			else {x_new[0]+=d;}
		}
		else if(a>=0.3333 && a<0.6666) { // asse y
			if(b<0.5000) {x_new[1]+=-d;}
				else {x_new[1]+=d;}
		}
		else {			         // asse z
			if(b<0.5000) {x_new[2]+=-d;}
			else {x_new[2]+=d;}
		}

//		alpha = fmin ( 1.0 , (distribution(x_new[0], x_new[1], x_new[2])/distribution(x[0], x[1], x[2])) );     // ground state
		alpha = fmin ( 1.0 , (distribution2(x_new[0], x_new[1], x_new[2])/distribution2(x[0], x[1], x[2])) );   // 2p

		if (g<=alpha) { 
			for(int j=0; j<3; j++) {
				x[j]=x_new[j];
			}
		conta2 ++;
		}
	
		dist += sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
	}

	cout << "With acceptance " << (double)conta2/(double)n*100 << "%, the averege distance in the last " << n << " steps is:   " << dist/n << "   Do you want to equilibrate again? [Y/N] " << endl;
	cin >> h;
	while (h != 'Y' && h != 'N') {
		cout << "Invalid input, try again with [Y/N] " << endl;
		cin >> h;
	}
	conta2=0;
}
///////////////////////////////////////////////////////////////////////////



//   Finita l'equilibrazione, inizia la simulazione Metropolis "vera"   //
for (int i=0; i<N; i++) {

	for(int j=0; j<3; j++) {
		x_new[j]=x[j];
	}

	a=rnd.Rannyu();
	b=rnd.Rannyu();
	g=rnd.Rannyu();
	d=rnd.Gauss(3,1); 		// da usare nel caso multivariato per l'orbitale 2p per avere 50% di accettanza 
	//d=rnd.Gauss(1.3,0.5); 	// da usare nel caso multivariato per il Ground State per avere 50% di accettanza
	//d=rnd.Rannyu(2.5,3.5); 	// da usare nel caso uniforme per l'orbitale 2p per avere 50% di accettanza 
	//d=rnd.Rannyu(0.75,1.55); 	// da usare nel caso uniforme per il Ground State per avere 50% di accettanza

	if(a<0.3333) {			 // asse x
		if(b<0.5000) {x_new[0]+=-d;}
		else {x_new[0]+=d;}
	}
	else if(a>=0.3333 && a<0.6666) { // asse y
		if(b<0.5000) {x_new[1]+=-d;}
		else {x_new[1]+=d;}
	}
	else {				 // asse z
		if(b<0.5000) {x_new[2]+=-d;}
		else {x_new[2]+=d;}
	}

//	alpha = fmin ( 1.0 , (distribution(x_new[0], x_new[1], x_new[2])/distribution(x[0], x[1], x[2])) );    // ground state
	alpha = fmin ( 1.0 , (distribution2(x_new[0], x_new[1], x_new[2])/distribution2(x[0], x[1], x[2])) );  // 2p

	if (g<=alpha) { 
		for(int j=0; j<3; j++) {
			x[j]=x_new[j];
		}
		conta ++;
	}
	
	r[i] = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
	out1 << x[0] << ";" << x[1] << ";" << x[2] << endl ;
}
////////////////////////////////////////////////////////////////////

// 	Datablocking	 //
for(int i=0; i<M; i++) {
	double sum=0;
	for (int j=0; j<L; j++) {
		int k=j+i*L;
		sum += r[k];
	}
		ave[i] = sum/L;
		ave2[i] = ave[i] * ave[i]; 
}

for(int i=0; i<M; i++) {
	for (int j=0; j<(i+1); j++) {
		sum_prog[i] += ave[j];
		sum2_prog[i] += ave2[j];
	}
	sum_prog[i] /= (i+1);
	sum2_prog[i] /= (i+1);
	err_prog[i] = error(sum_prog[i], sum2_prog[i], i) ;
}


for (int i=0; i<M; i++) {
	out2 << (i+1) << ";" << sum_prog[i] << ";" << err_prog[i] << endl;
}


cout << "Percentuale passi accettati: " << (double)conta/(double)N*100 << " %" << endl;
out1.close();
out2.close();

rnd.SaveSeed();
return 0;
}



double error(double AV, double AV2, int n) {
if(n==0) { return 0; }
else { return sqrt((AV2-AV*AV)/n); }
}

double distribution(double _x, double _y, double _z) {
double _r = sqrt(_x*_x+_y*_y+_z*_z);
return pow(exp(-_r)/sqrt(M_PI),2); 
}
double distribution2(double _x, double _y, double _z) {
double cost = sqrt(2/M_PI)*(1./8.);
double _r = sqrt(_x*_x+_y*_y+_z*_z);
return pow(cost * _r * exp(-_r/2) * (_z/_r),2); 
}



/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 5
			 ESERCIZIO 5.1

*****************************************************************/
