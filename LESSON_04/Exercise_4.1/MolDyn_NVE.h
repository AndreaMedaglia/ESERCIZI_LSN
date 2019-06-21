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
#include <string>
using namespace std;


//parameters, observables
const int m_props=105;
double walker[m_props];
int n_props, igofr;
int iv,ik,it,ie,ip;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
double err_pot, err_kin, err_etot, err_temp, err_press;
double vtail, ptail;

//pigreco
const double pi=3.1415927;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part],_xold[m_part],_yold[m_part],_zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, n_block, n_data, n, l_block, nbins;
double delta, bin_size;
char b;

//functions
void Input(void);
void Rescale(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure();
double Force(int, int);
double Pbc(double);
double Get_temp(void);
double Error(double,double,int);
void Averages(int);
void Accumulate(void);
void Print(int t);
void Reset(int);

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
