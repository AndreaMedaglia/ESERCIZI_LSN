/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 10
			 ESERCIZIO 10.2

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
void print_conf(vector<int> &v, int n);
double chromo_length(vector<int> &v);
void check(vector<int> &v);
int selection(vector<Chromosome> &v);
void permutation(vector<int> &v);
void shift(vector<int> &v);


// parameters
int N_cities, selected, N_iter, N_in, N_metro, conta, N_step;
double x_start, y_start, sh, theta, sum, T, T_fin, beta, beta_fin, p, alpha, incr;
vector<double> x_cities;
vector<double> y_cities;
vector<Chromosome> pop_old;
vector<Chromosome> pop_new;



/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 10
			 ESERCIZIO 10.2

*****************************************************************/
