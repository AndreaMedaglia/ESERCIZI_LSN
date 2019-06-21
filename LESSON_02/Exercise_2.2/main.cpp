/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 2
			 ESERCIZIO 2.2

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

double error(double AV, double AV2, int n);

using namespace std;
 
int main (int argc, char *argv[]){

   int M=1000;  // NUMERO DI RAMDOM WALK 
   int N=100;   // NUMERO DI STEP PER OGNI RANDOM WALK
  
   ofstream out ("Output/out.txt"); // output generico, rinominare per fare i grafici

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


double r[N], r2[N],e[N];
double theta, cos_phi, a, p;

double x[3][M];   // x[asse][numero di RW]
for(int i=0; i<M; i++) {
	x[0][i]=0;
	x[1][i]=0;
	x[2][i]=0;
}


// Commentare/scommentare in base a quale Random Walk si vuole simulare (continuo o discreto) 

for (int i=0; i<N; i++) {
	double sum=0;
	double sum2=0;
	for(int j=0; j<M; j++) {

		// RANDOM WALK NEL DISCRETO  // 
		/*
		a=rnd.Rannyu();
		p=rnd.Rannyu();
		if(a<0.33333) {			  		// asse x
			if(p<0.5) {x[0][j]+=-1;}
			else {x[0][j]+=1;}
		}
		else if(a>=0.33333 && a<0.66666) {		// asse y
			if(p<0.5) {x[1][j]+=-1;}
			else {x[1][j]+=1;}
		}
		else {						// asse z
			if(p<0.5) {x[2][j]+=-1;}
			else {x[2][j]+=1;}
		}
		*/ 
		//////////////////////////////


		// RANDOM WALK NEL CONTINUO  //
		theta=rnd.Rannyu(0, 2*M_PI);
		cos_phi=rnd.Rannyu(-1, 1);
		x[0][j]+=cos(theta)*sin(acos(cos_phi));
		x[1][j]+=sin(theta)*sin(acos(cos_phi));
		x[2][j]+=cos_phi;
		///////////////////////////////

		sum+=sqrt(x[0][j]*x[0][j]+x[1][j]*x[1][j]+x[2][j]*x[2][j]);
		sum2+=sqrt(x[0][j]*x[0][j]+x[1][j]*x[1][j]+x[2][j]*x[2][j])*sqrt(x[0][j]*x[0][j]+x[1][j]*x[1][j]+x[2][j]*x[2][j]);
	}
	r[i]=sum/M;
	r2[i]=sum2/M;
	e[i] = error(r[i], r2[i], M);
}




for (int i=0; i<N; i++) {
	out << (i+1)  << ";" << r[i] << ";" << e[i] << endl;
}

out.close();

   rnd.SaveSeed();
   return 0;
}



double error(double AV, double AV2, int n) {
if(n==0) { return 0; }
else { return sqrt((AV2-AV*AV)/n); }
}



/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 2
			 ESERCIZIO 2.2

*****************************************************************/

