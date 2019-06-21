/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 3
			 ESERCIZIO 3.1

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

double error(double AV, double AV2, int n);

using namespace std;
 
int main (int argc, char *argv[]){

   ofstream outC ("Output/Output_C.txt"); // rinominare gli output per visualizzare nel jupyter-notebook
   ofstream outP ("Output/Output_P.txt");

	int M=10000;
	int N=100;
	int L=M/N;

	double S_0=100;
	double T=1;
	double K=100;
	double r=0.1;
	double sigma=0.25;
	
  

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


	double S[M], C[M], P[M];
	double delta=T/100.;
	double s_discrete;
	double W;


	// Commentare/scommentare in base al caso che si vuole simulare (diretto, discretizzato)//

	/*
	//       Put e Call caso diretto      //
	for (int i=0; i<M; i++) {
		S[i]=S_0 * exp((r-sigma*sigma*0.5)*T+sigma*rnd.Gauss(0,T));
		P[i]=exp(-r*T)*fmax(0,(K-S[i]));
		C[i]=exp(-r*T)*fmax(0,(S[i]-K));
	}
	////////////////////////////////////////
	*/

	//       Put e Call caso discreto     //
	for (int i=0; i<M; i++) {
		s_discrete = S_0;
		for (int j=1; j<100; j++){
			W = rnd.Gauss(0,T);
			s_discrete *= exp((r-sigma*sigma*0.5)*delta+sigma*W*sqrt(delta));
		}
		S[i]=s_discrete;
		C[i]=exp(-r*T)*fmax(0,S[i]-K);
		P[i]=exp(-r*T)*fmax(0,K-S[i]);

	}
	////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////

   double aveC[N], ave2C[N], sum_progC[N], sum2_progC[N], err_progC[N];
   double aveP[N], ave2P[N], sum_progP[N], sum2_progP[N], err_progP[N];

for(int i=0; i<N; i++) { // inizializzazione vettori (per evitare problemi)
    sum_progC[i]=0;
    sum2_progC[i]=0;
    err_progC[i]=0;
    aveC[i]=0;
    ave2C[i]=0;

    aveP[i]=0;
    ave2P[i]=0;
    sum_progP[i]=0;
    sum2_progP[i]=0;
    err_progP[i]=0;
}


	for(int i=0; i<N; i++) {
		double sumC=0;
		double sumP=0;
		for (int j=0; j<L; j++) {
			int k=j+i*L;
			sumC += C[k];
			sumP += P[k];
		}
			aveC[i] = sumC/L;
			ave2C[i] = aveC[i] * aveC[i];
			aveP[i] = sumP/L;
			ave2P[i] = aveP[i] * aveP[i]; 
	}

	for(int i=0; i<N; i++) {
		for (int j=0; j<(i+1); j++) {
			sum_progC[i] += aveC[j];
			sum2_progC[i] += ave2C[j];
			sum_progP[i] += aveP[j];
			sum2_progP[i] += ave2P[j];
		}
		sum_progC[i] /= (i+1);
		sum2_progC[i] /= (i+1);
		err_progC[i] = error(sum_progC[i], sum2_progC[i], i) ;
		sum_progP[i] /= (i+1);
		sum2_progP[i] /= (i+1);
		err_progP[i] = error(sum_progP[i], sum2_progP[i], i) ;
	}


	for (int i=0; i<N; i++) {
		outC << (i+1) << ";" << sum_progC[i] << ";" << err_progC[i] << endl;
		outP << (i+1) << ";" << sum_progP[i] << ";" << err_progP[i] << endl;
 	}


	outC.close();
	outP.close();

   rnd.SaveSeed();
   return 0;
}



double error(double AV, double AV2, int n) {

if(n==0) { return 0; }
else { return sqrt((AV2-AV*AV)/n); }

}



/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 3
			 ESERCIZIO 3.1

*****************************************************************/

