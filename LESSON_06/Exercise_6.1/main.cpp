/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 6
			 ESERCIZIO 6.1

*****************************************************************/


#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() {        

Input();

///// 			equilibrazione	 		/////
cout << endl << endl << "Do you want to make an equilibration of the system? [Y/N]" << endl;
cin >> q;
while (q != 'Y' && q != 'N') {
	cout << "Invalid input, try again with [Y/N] " << endl;
	cin >> q;
}
while (q == 'Y') { 

	cout << "For how many Metropolis/Gibbs steps you want to simulate? [integer]" << endl;
	cin >> n; 

	for(int istep=1; istep <= n; ++istep) {
		Move(metro);
	}

	cout << "Do you want to equilibrate again? [Y/N] " << endl;
	cin >> q;
	while (q != 'Y' && q != 'N') {
		cout << "Invalid input, try again with [Y/N] " << endl;
		cin >> q;
	}
}


///// 			simulazione	 		/////
if(sim==0) {
	ofstream E_int, Susc, Heat, Magn;
	E_int.open("Output/Internal_energy.dat");
	Susc.open("Output/Susceptivity.dat");
	Heat.open("Output/Heat_capacity.dat");
	Magn.open("Output/Magnetization.dat");
					
	while(temp<T_fin){
		Inizialization();
		for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
			Reset(iblk);   			//Reset block averages
			for(int istep=1; istep <= nstep; ++istep) {
			Move(metro);
			Measure();
			Accumulate(); 			//Update block averages
			}
		Averages(iblk);   			//Print results for current block
		}
		E_int << temp <<  ";" << glob_av[iu]/(double)nblk << ";" << err_u << endl;
		Susc << temp <<  ";" << glob_av[ix]/(double)nblk << ";" << err_x << endl;
		Heat << temp <<  ";" << glob_av[ic]/(double)nblk << ";" << err_c << endl; 
		Magn << temp << ";" << glob_av[im]/(double)nblk << ";" << err_m << endl; 
		Reset_conf(); 				// Reset configuration to restart with another Temperature
		temp = temp + incr;
		beta = 1.0/temp;
	}
	E_int.close();
	Susc.close();
	Heat.close();
	Magn.close();
}

else{
	for(int iblk=1; iblk <= nblk; ++iblk) {		//Simulation
		Reset(iblk);   				//Reset block averages
		for(int istep=1; istep <= nstep; ++istep) {
			Move(metro);
			Measure();
			Accumulate();			//Update block averages
		}
		Averages(iblk);   			//Print results for current block
	}
	ConfFinal(); 					//Write final configuration
}


return 0;
}

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> T_fin;
  cout << "Final Temperature = " << T_fin << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> incr;

  ReadInput >> sim;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  if(sim==1) cout << "The program is simulated at the same temperature" << endl;
  else { cout << "The program is simulated at the different temperature with increment "<< incr << endl; }

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double r, alpha, s_new;


  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    attempted++;

    if(metro==1){ //Metropolis     
	r=rnd.Rannyu();
	s_new = s[o]*(-1);
	alpha = fmin( 1, exp(-beta*(Boltzmann(s_new, o)-Boltzmann(s[o], o))) );
	if(r<=alpha) {
		s[o]=s_new;
		accepted++;
	}
    }
    else { //Gibbs sampling 
	r=rnd.Rannyu();    
	alpha = 1/( 1 + exp(-2*beta*J*( s[Pbc(o-1)] + s[Pbc(o+1)] ) - 2*beta*h) ) ;
	if(r<=alpha) { s[o]=+1; }
	else {s[o]=-1;} 
	accepted++;
    }
  }
}

double Boltzmann(double sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
//  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;    // 0 energia
  walker[ic] = u*u;  // 1 capacità termica
  walker[im] = m;    // 2 magnetizzazione
  walker[ix] = m*m;  // 3 suscettività

}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

    stima_u = blk_av[iu]/blk_norm/(double)nspin; 
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);

    stima_c = beta * beta * (blk_av[ic]/blk_norm - blk_av[iu]/blk_norm*blk_av[iu]/blk_norm)/(double)nspin; 
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);

    stima_m = blk_av[im]/blk_norm/(double)nspin; 
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);

    stima_x = beta * blk_av[ix]/blk_norm/(double)nspin; 
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);

	if (sim==1) {    
		ofstream Ene, Heat, Mag, Chi;
		const int wd=12;
	    
		cout << "Block number " << iblk << endl;
		cout << "Acceptance rate " << accepted/attempted << endl << endl;

		Ene.open("Output/output.ene.0",ios::app);
		Heat.open("Output/output.heat.0",ios::app);
		Mag.open("Output/output.mag.0",ios::app);
		Chi.open("Output/output.chi.0",ios::app);

		Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
		Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
		Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl; 
		Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;

		Mag.close();
		Heat.close();
		Ene.close();
		Chi.close();

		cout << "----------------------------" << endl << endl;
	}
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("Output/config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Reset_conf(void) {
for (int i=0; i<nspin; ++i) {
	s[i] = 0;
}
}

void Inizialization(void) {
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }

cout << "Current temperature " << temp << endl;

}


/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 6
			 ESERCIZIO 6.1

*****************************************************************/


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
