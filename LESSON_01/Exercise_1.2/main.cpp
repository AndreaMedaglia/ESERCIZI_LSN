/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 1
			 ESERCIZIO 1.2

*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"


using namespace std;
 
int main (int argc, char *argv[]){

   int N= 10000 ; // numero di medie 
   int L= 100 ;   // elementi per ogni media		<--------cambiare per ottenere i diversi istogrammi (L=1,2,10,100)
   int M= N*L ;   // numero totale di Random Variabili

   double l= 1 ; 
   double gamma=1;  
   double x_0=0;
  
   ofstream out ("Output/out.txt"); // l'output è solo uno, rinominare in modo opportuno ogni volta che cambi L 


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


double r[M];
double ave[N];

	for(int i=0; i<M; i++){					// commentare/scommentare in base alla distribuzione che si vuole usare:
//		r[i] = rnd.Exp(l);				// esponenziale
//		r[i] = rnd.CauchyLorentz(x_0, gamma);		// di cauchy-lorentz	
		r[i] = rnd.Rannyu(); 				// uniforme
	}

	for(int i=0; i<N; i++) {
		double sum=0;
		for (int j=0; j<L; j++) {
			int k=j+i*L;
			sum += r[k];
		}
		ave[i] = sum/L;
	}

	for(int i=0; i<N; i++) {
		out << ave[i] << endl ;
	}


	out.close();

	rnd.SaveSeed();

   return 0;
}



/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 1
			 ESERCIZIO 1.2

*****************************************************************/

