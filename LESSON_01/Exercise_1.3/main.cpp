/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 1
			 ESERCIZIO 1.3

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

double error(double AV, double AV2, int n);

using namespace std;
 
int main (int argc, char *argv[]){

   int M=100000;  // numeri generati
   int N=100;	  // numero di blocchi
   int L=M/N;     // elementi per blocco
   int Q=10;      // intervalli in cui viene diviso [0,1]
   double d=0.1;  // larghezza intervalli
   double l=0.07; // lunghezza ago di Buffon


   double x_1[M], x_2[M], pi[N], pi2[N];
   int conta[N];
   double sum_prog[N], sum2_prog[N], err_prog[N];

   ofstream out ("Output.txt");


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


 
	double x, y;

	for(int i=0; i<M; i++){
		x_1[i] = rnd.Rannyu();
		do {
		y = rnd.Rannyu();
		x = rnd.Rannyu();
		} while ((x*x+y*y)>1);
		x_2[i] = x_1[i]+ l * cos(atan(y/x));
		if (x_2[i]<0) {
			x_2[i]=0;	
		}
		if (x_2[i]>1) {
			x_2[i]=1;	
		}
	}


	for (int i=0; i<N; i++) {
		conta[i]=0;
	}


for (int i=0; i<N; i++) {

	for (int j=0; j<L; j++) {
		int k=j+i*L;
		for(int w=0; w<Q; w++) {
			double a=(double)w*1/10;
			double b=(double)(w+1)*1/10;
			if (x_1[k]>=a && x_1[k]<=b) {
				if (x_2[k]<=a || x_2[k]>=b) {
					conta[i]++;
				}
			} 
		}
	}
		
	pi[i]=(2*l/d)*((double)L/conta[i]);
	pi2[i]=pi[i]*pi[i];

}

	for(int i=0; i<N; i++) {
		for (int j=0; j<(i+1); j++) {
			sum_prog[i] += pi[j];
			sum2_prog[i] += pi2[j];
		}
		sum_prog[i] /= (i+1);
		sum2_prog[i] /= (i+1);
		err_prog[i] = error(sum_prog[i], sum2_prog[i], i) ;
	}


	for (int i=0; i<N; i++) {
		out << (i+1) << ";" << sum_prog[i] << ";" << err_prog[i] << endl;
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
			   LEZIONE 1
			 ESERCIZIO 1.3

*****************************************************************/

