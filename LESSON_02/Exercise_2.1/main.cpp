/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 2
			 ESERCIZIO 2.1

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

double error(double AV, double AV2, int n);

using namespace std;
 
int main (int argc, char *argv[]){

   int M=10000;    // numeri casuali
   int N=100;	   // numero di blocchi
   int L=M/N;      // numero di elementi per blocco

   ofstream out ("Output/out.txt"); // l'output è uno solo, rinominare per fare i grafici in jupyter-notebook

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


   double r[M], ave[N], ave2[N], sum_prog[N], sum2_prog[N], err_prog[N];
   double x, p; 



///// commentare/scommentare in base alla tecnica e alla distribuzione che si vuole usare /////


//     Distribuzione uniforme     //
	for(int i=0; i<M; i++){
		r[i] = rnd.Rannyu();
	}
///////////////////////////////////


/*
//  Importance sampling tecnique  // 
	for(int i=0; i<M; i++){
		do {
			x=rnd.Rannyu();
			p=rnd.Rannyu();
//		} while( p>( (1-x) ) ); 								// RETTA
		} while( p>( ((M_PI/2-tan(x))/(M_PI/2+log(cos(1))))/(M_PI/(M_PI+2*log(cos(1)))) ) );    // PARABOLA 

		r[i] = x;
	}
///////////////////////////////////
*/

///////////////////////////////////////////////////////////////////////////////////////////////


// Datablocking //
// commentare/scommentare sum+= in base alla tecnica usata in precedenza //
	for(int i=0; i<N; i++) {
		double sum=0;
		for (int j=0; j<L; j++) {
			int k=j+i*L;
			sum += (M_PI/2)*cos(M_PI*r[k]/2);						// UNIFORME
//			sum +=  M_PI/2*cos(M_PI/2*r[k]) / (2*(1-r[k])) ;  				// RETTA
//			sum += M_PI/2*cos(M_PI/2*r[k]) / ((M_PI/2-tan(r[k]))/(M_PI/2+log(cos(1))));	// PARABOLA   
		}
			ave[i] = sum/L;
			ave2[i] = ave[i] * ave[i];
	}


	for(int i=0; i<N; i++) {		
		for (int j=0; j<(i+1); j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] /= (double)(i+1);
		sum2_prog[i] /= (double)(i+1);
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
			   LEZIONE 2
			 ESERCIZIO 2.1

*****************************************************************/


