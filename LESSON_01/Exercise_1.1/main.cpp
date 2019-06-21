/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 1
			 ESERCIZIO 1.1

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

double error(double AV, double AV2, int n);

using namespace std;
 
int main (int argc, char *argv[]){

   int M=1000000; // quantità di numeri pseudo-casuali generati
   int N=100;	  // numero di blocchi
   int L=M/N;     // numero di elementi per blocco
   int Q=L/N;     // numero di sottointervalli di [0,1] // per la stima di chi^2

   ofstream out1 ("Output/Stima_r.txt");
   ofstream out2 ("Output/Stima_sigma.txt");
   ofstream out3 ("Output/Stima_chi.txt");


   Random rnd;  // generatore di numeri pseudo-casuali
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


   // variabili per il datablocking
   double r[M];
   double ave_r[N], ave2_r[N], sum_prog_r[N], sum2_prog_r[N], err_prog_r[N];
   double ave_s[N], ave2_s[N], sum_prog_s[N], sum2_prog_s[N], err_prog_s[N];
   double ave_c[N], ave2_c[N], sum_prog_c[N], sum2_prog_c[N], err_prog_c[N];
   int conta[Q];


	for(int i=0; i<M; i++){
		r[i] = rnd.Rannyu();
	}

	for(int i=0; i<N; i++) {
		double sum_r=0;
		double sum_s=0;
		double sum_c=0;
		for(int w=0; w<Q; w++) {
			conta[w]=0;
		}
		for (int j=0; j<L; j++) {
			int k=j+i*L;
			for(int w=0; w<Q; w++) {              
				double a=((double)w)/Q;
				double b=((double)w+1)/Q;
				if (r[k]>a && r[k]<b) {conta[w]++;}
			}			     
			sum_r += r[k];                          // somma per la stima di r
			sum_s += (r[k] - 0.5) * (r[k] - 0.5);   // somma per la stima di sigma^2 
		}
		for (int j=0; j<Q; j++) {			// somma per la stima di chi^2
			sum_c += (conta[j] - Q) * (conta[j] - Q) ;
		}
			ave_r[i] = sum_r/L;
			ave2_r[i] = ave_r[i] * ave_r[i]; 
			ave_s[i] = sum_s/L;
			ave2_s[i] = ave_s[i] * ave_s[i]; 
			ave_c[i] = sum_c/Q;
			ave2_c[i] = ave_c[i] * ave_c[i]; 
	}

	for(int i=0; i<N; i++) {
		for (int j=0; j<(i+1); j++) {
			sum_prog_r[i] += ave_r[j];
			sum2_prog_r[i] += ave2_r[j];
			sum_prog_s[i] += ave_s[j];
			sum2_prog_s[i] += ave2_s[j];
			sum_prog_c[i] += ave_c[j];
			sum2_prog_c[i] += ave2_c[j];
		}
		sum_prog_r[i] /= (i+1);
		sum2_prog_r[i] /= (i+1);
		err_prog_r[i] = error(sum_prog_r[i], sum2_prog_r[i], i) ;
		sum_prog_s[i] /= (i+1);
		sum2_prog_s[i] /= (i+1);
		err_prog_s[i] = error(sum_prog_s[i], sum2_prog_s[i], i) ;
		sum_prog_c[i] /= (i+1);
		sum2_prog_c[i] /= (i+1);
		err_prog_c[i] = error(sum_prog_c[i], sum2_prog_c[i], i) ;
	}


	for (int i=0; i<N; i++) {
		out1 << (i+1) << ";" << sum_prog_r[i] << ";" << err_prog_r[i] << endl;
		out2 << (i+1) << ";" << sum_prog_s[i] << ";" << err_prog_s[i] << endl;
		out3 << (i+1) << ";" << sum_prog_c[i] << ";" << err_prog_c[i] << endl;
 	}

	out1.close();
	out2.close();
	out3.close();

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
			 ESERCIZIO 1.1

*****************************************************************/

