/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 9
			 ESERCIZIO 9

*****************************************************************/
#include <vector>
using namespace std;

struct Chromosome {
	vector<int> path;
	double length;
};

bool operator<(const Chromosome &a, const Chromosome &b){
	return a.length<b.length; 
}

// functions
void input(void);
void inizialization_square(void);
void inizialization_circle(void);
void conf_0(vector<int> &v);
void print_conf(vector<int> &v);
void mutation(vector<int> &v);
void permutation(vector<int> &v);
void shift(vector<int> &v);
double chromo_length(vector<int> &v);
void check(vector<int> &v);
void inversion(vector<int> &v);
int selection(vector<Chromosome> &v);
void crossover(vector<int> &v, vector<int> &w);


// parameters
int N_cities, N_pop, N_mut, selected, N_iter, N_in;
double x_start, y_start, mut, cross, sh, theta, inv, p_crossover, p_permutation, p_shift, p_inversion, sum;
double variabile_di_prova;
vector<double> x_cities;
vector<double> y_cities;
vector<Chromosome> pop_old;
vector<Chromosome> pop_new;



/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 9
			 ESERCIZIO 9

*****************************************************************/
