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
			   LEZIONE 4
			 ESERCIZIO 4.1

*****************************************************************/

#include <stdlib.h>    
#include <iostream>     
#include <fstream>     
#include <cmath>        
#include "MolDyn_NVE.h"
#include <cstdlib>
#include <iomanip>
#include <string>
#include "random.h"

using namespace std;

int main(int argc, char *argv[]){ 

Input();
cout << n_data;
	
///// 				EQUILIBRAZIONE 				/////
cout << endl << endl << "****** PRELIMINARY PART ******" << endl << endl;
cout << "If you read the old configuration from old.0, do you want to rescale ? [Y/N] " << endl;
cout << "[Remember that the old configuration generated from the random velocities is ALREADY rescaled]" << endl;
cin >> b;
while (b != 'Y' && b != 'N') {
	cout << "Invalid input, try again with [Y/N] " << endl;
	cin >> b;
}
if(b=='Y') { Rescale();}
cout << "Do you want to make an equilibration of the system? [Y/N]" << endl;
cin >> b;
int cont=0;
while (b == 'Y') { 
	if (cont!=0) {Rescale();} 
	cout << "For how many steps you want to simulate (and then rescale)? [integer]" << endl;
	cin >> n; 
	for (int i=0; i<n; i++) { Move(); }
	cout << "Now the temperature is:   " << Get_temp() << "   Do you want to rescale? [Y/N] " << endl;
	cin >> b;
	cont++;
}
///// 			        FINE EQUILIBRAZIONE 			/////


/////				SIMULAZIONE        			/////
cout << endl << endl << "****** SIMULATION ******" << endl << endl;
int nconf = 1;
for(int iblk=1; iblk <= n_block; ++iblk) {
	Reset(iblk);				
	for(int istep=1; istep <= nstep; ++istep) {
		Move();
		if(istep%10 == 0){
			Measure();
			Accumulate();			
			Print(nstep*(iblk-1)+istep);	
		} 
		if(istep%10 == 0){
			//ConfXYZ(nconf);	 
			nconf += 1;
		}
	}
	Averages(iblk);			
	if (iblk%10==0) cout << "block number: " << iblk << endl;
}
/////				FINE SIMULAZIONE        		/////

ConfFinal(); 

return 0;
}




///// 				FUNCTIONS				/////


// Input: prepara tutti i parametri per la simuazione, prepara le osservabili, legge configurazione iniziale e legge/genera quella vecchia
void Input(void){ 

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

	ifstream ReadInput,ReadConf,ReadOld;
	//double ep, ek, pr, et, vir;
	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl << endl;

	  
	ReadInput.open("Input/input.dat"); 
	ReadInput >> temp;
	ReadInput >> npart;
	cout << "****** PARAMETERS ******" << endl << endl;
	cout << "Number of particles = " << npart << endl;
	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;
	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	n_data=(int)nstep/10;
	ReadInput >> n_block;
	l_block=n_data/n_block;
	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl;
	cout << "Number of block = " << n_block << endl << endl;
	ReadInput.close();

	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	ip = 4; //Pressure
	n_props = 5; //Number of observables

	//measurement of g(r)
	igofr = 5;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box/2.0)/(double)nbins;

	vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
	ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));


	//Read initial configuration
	cout << "****** ACTUAL AND OLD CONFIGURATION ******" << endl << endl;
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("Input/config.0");
	for (int i=0; i<npart; ++i){
		x[i]=0; y[i]=0; z[i]=0;
	}

	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}

	ReadConf.close();

	//Read or generate old configuration
	char a;
	cout << "Old Configuration" << endl;
	cout << "Do you want to start with the configuration saved in file old.0 (in any) ? [O] " << endl;
	cout << "Or do you want to start with an old configuration calculate from random velocities ? [R] " << endl;
	cin >> a;
	while (a != 'O' && a != 'R') {
		cout << "Invalid input, try again with [O/R] " << endl;
		cin >> a;
	}

	if (a == 'O') {
		ReadOld.open("Input/old.0");

		if(ReadOld.is_open()) {
			cout << "Read old configuration from file old.0 " << endl << endl;
			for (int i=0; i<npart; ++i){
				ReadOld >> xold[i] >> yold[i] >> zold[i];
				xold[i] = xold[i] * box;
				yold[i] = yold[i] * box;
				zold[i] = zold[i] * box;	
			}
		}
		else {
			cout << "ATTENTION!!! There is no file old.0" << endl;
			cout  << "Prepare old configuration from random velocities with center of mass velocity equal to zero " << endl << endl;
			double sumv[3] = {0.0, 0.0, 0.0};
			for (int i=0; i<npart; ++i){
				vx[i]=0; vy[i]=0; vz[i]=0;
				xold[i]=0; yold[i]=0; zold[i]=0;
			}
			for (int i=0; i<npart; ++i){
				vx[i] = rnd.Rannyu(-0.5, 0.5);
				vy[i] = rnd.Rannyu(-0.5, 0.5);
				vz[i] = rnd.Rannyu(-0.5, 0.5);

				sumv[0] += vx[i];
				sumv[1] += vy[i];
				sumv[2] += vz[i];
			}
			for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
			double sumv2 = 0.0, fs;
			for (int i=0; i<npart; ++i){
				vx[i] = vx[i] - sumv[0];
				vy[i] = vy[i] - sumv[1];
				vz[i] = vz[i] - sumv[2];

				sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
			}
			sumv2 /= (double)npart;

			fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
			for (int i=0; i<npart; ++i){
				vx[i] *= fs;
				vy[i] *= fs;
				vz[i] *= fs;

				xold[i] = Pbc( x[i] - vx[i] * delta);
				yold[i] = Pbc( y[i] - vy[i] * delta);
				zold[i] = Pbc( z[i] - vz[i] * delta);
			}
		}
	} 

	else {
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i=0; i<npart; ++i){
			vx[i]=0; vy[i]=0; vz[i]=0;
			xold[i]=0; yold[i]=0; zold[i]=0;
		}
		for (int i=0; i<npart; ++i){
			vx[i] = rnd.Rannyu(-0.5, 0.5);
			vy[i] = rnd.Rannyu(-0.5, 0.5);
			vz[i] = rnd.Rannyu(-0.5, 0.5);

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		double sumv2 = 0.0, fs;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
			
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc( x[i] - vx[i] * delta);
			yold[i] = Pbc( y[i] - vy[i] * delta);
			zold[i] = Pbc( z[i] - vz[i] * delta);
		}
	}


	rnd.SaveSeed();

return;
}
// fine Input()


// Move: effettua un passo Verlet
void Move(void){ 
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){ 
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){ 

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		_xold[i]=xold[i];
		_yold[i]=yold[i];
		_zold[i]=zold[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}

	return;
}



//Force: calcola la forza come f=-Grad_ip V(r)	
double Force(int ip, int idir){ 	
	double f=0.0;
	double dvec[3], dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );  
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );

			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); 
			}
		}
	}
	  
	return f;
}



// Misura le osservabili in esame 
void Measure(){ 

	double v, t, vij, w, wij;
	double dx, dy, dz, dr;


	v = 0.0; //reset observables
	w = 0.0;
	t = 0.0;

	//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){

			dx = Pbc( x[i] - x[j] );
			dy = Pbc( y[i] - y[j] );
			dz = Pbc( z[i] - z[j] );

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

		 	for (int l=0; l<nbins; l++) {
				if ( (dr)>=l*bin_size && (dr)<(l+1)*bin_size ) {
					walker[igofr+l]+=2;
				}
			}

			if(dr < rcut){
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6); 
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

				v += vij;
				w += wij;
			}
		}          
	}

	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); 
	   
	walker[iv] = 4 * v;			//Potential energy
	walker[ik] = t;				//Kinetic energy
	walker[it] = (2.0 / 3.0) * t;		//Temperature
	walker[ie] = walker[iv] + walker[ik];	//Total energy
	walker[ip] = 48.0 * w / 3.0;		//Pressure

	return;
}



//Reset: inizializza a zero le variabili necessarie per le medie
void Reset(int iblk) 
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



//Accumulate: aggiorna le medie sui blocchi
void Accumulate(void) 
{

   for(int i=0; i<(n_props); ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}



//Averages: stampa i valori delle osservabili relativi al blocco in esame
void Averages(int iblk) 
{
    
   ofstream Gofr, Gave, Epot, Pres, Ekin, Etot, Temp;
   const int wd=12;
   int nbins=100;
   double stima_g[nbins], err_g[nbins];

    
//    cout << "Block number " << iblk << endl;
//    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("Output/blocking.epot.dat",ios::app);
    Pres.open("Output/blocking.pres.dat",ios::app);
    Ekin.open("Output/blocking.ekin.dat",ios::app);
    Etot.open("Output/blocking.etot.dat",ios::app);
    Temp.open("Output/blocking.temp.dat",ios::app);
    Gofr.open("Output/blocking.gofr.dat",ios::app);
    Gave.open("Output/blocking.gave.dat");
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail;				//Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_temp = blk_av[it]/blk_norm/(double)npart;					//Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    stima_kin = blk_av[ik]/blk_norm/(double)npart;					//Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
    
    stima_press = rho * stima_temp + (blk_av[ip]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[ip] += stima_press;
    glob_av2[ip] += stima_press*stima_press;
    err_press=Error(glob_av[ip],glob_av2[ip],iblk);

    stima_etot = stima_pot + stima_kin;					//Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    for (int i=0; i<nbins; i++) {
	stima_g[i] = blk_av[i+igofr]/blk_norm/(rho*(double)npart*(4.0*pi/3.0)* (pow((bin_size*(i+1)),3)-pow((bin_size*i),3)) ) ;
	glob_av[i+igofr] += stima_g[i];
	glob_av2[i+igofr] += stima_g[i]*stima_g[i];
	err_g[i]=Error(glob_av[i+igofr],glob_av2[i+igofr],iblk);
    }

    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
    Pres << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_press << endl;
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//g(r)
	Gofr << setw(wd) << iblk << endl << endl;
    for (int i=0; i<nbins; i++) {
	Gofr  <<  i+1 << " " << stima_g[i] << " " << glob_av[i+igofr]/(double)iblk << " " << err_g[i] << endl;
    }

    for (int i=0; i<nbins; i++) {
	Gave  <<  bin_size/2+i*bin_size << " " << stima_g[i] << " " << glob_av[i+igofr]/(double)iblk << " " << err_g[i] << endl;
    }    

//    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Gofr.close();
    Gave.close();
}

//Error: calcola l'errore sul blocco
double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


//ConfFinal: stampa configurazione al tempo finale t_f e configurazione al tempo t_f - dt 
void ConfFinal(void){ 

	ofstream WriteConf, WriteOld;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("Output/config.final");

	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteConf.close();

	cout << "Print old configuration to file old.final " << endl << endl;
	WriteOld.open("Output/old.final");

	for (int i=0; i<npart; ++i){
		WriteOld << _xold[i]/box << "   " <<  _yold[i]/box << "   " << _zold[i]/box << endl;
	}
	WriteOld.close();

	return;
}


//ConfXYZ: Write configuration in .xyz format
void ConfXYZ(int nconf){ 
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}

	WriteXYZ.close();
}



// Pbc: applica le condizioni al contorno
double Pbc(double r){  
	return r - box * rint(r/box);
}


// Print: stampa valori istantanei delle osservabili
void Print(int t){		
	ofstream E_pot, Press, E_kin, E_tot, Temperature;
	E_pot.open("Output/Potential_energy_inst.dat",ios::app);
	Press.open("Output/Pressure_inst.dat",ios::app);
	E_kin.open("Output/Kinetic_energy_inst.dat",ios::app);
	E_tot.open("Output/Total_energy_inst.dat",ios::app);
	Temperature.open("Output/Temperature_inst.dat",ios::app);

	E_pot << t << ";" << (walker[iv] / (double)npart + vtail) << endl;
	Press << t << ";" << rho * walker[it]/(double)npart + ( walker[ip] + ptail*(double)npart) / vol << endl;
	E_kin << t << ";" << (walker[ik] / (double)npart)  << endl;
	E_tot << t << ";" << (walker[ie] / (double)npart) << endl;
	Temperature << t << ";" << (walker[it] / (double)npart) << endl;

	E_pot.close();
	Press.close();
	E_kin.close();
	E_tot.close();
	Temperature.close();
}


// Rescale: funzione che riscala in base alla temperatura bersaglio le velocità (e di conseguenza anche le posizioni)
void Rescale (void) {

	Move();
	double sumv[3] = {0.0, 0.0, 0.0};   // riscalamento //

	for (int i=0; i<npart; ++i){
		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];
	}

	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

	double sumv2 = 0.0, fs;
	for (int i=0; i<npart; ++i){
		vx[i] =  vx[i] - sumv[0];
		vy[i] =  vy[i] - sumv[1];
		vz[i] =  vz[i] - sumv[2];
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}

	sumv2 /= (double)npart;
	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	for (int i=0; i<npart; ++i){
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;
		xold[i] = Pbc( x[i] - vx[i] * delta);
		yold[i] = Pbc( y[i] - vy[i] * delta);
		zold[i] = Pbc( z[i] - vz[i] * delta); 
	}

	return;
}


// Get_temp(): funzione per il calcolo della temperatura istantanea, serve nell'equilibrazione per avere un'idea della temperatura del sistema
double Get_temp (void) {
	double temperature=0.0;
	for (int i=0; i<npart; ++i) temperature += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	return ( (2.0 / 3.0) * temperature/(double)npart );
}


/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 4
			 ESERCIZIO 4.1

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

