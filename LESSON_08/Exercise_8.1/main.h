/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 8
			 ESERCIZIO 8.1

*****************************************************************/
#include <vector>
using namespace std;

struct Variational {
	double en_best;
	vector<double> energy;
	vector<double> error;
	vector<double> isto;
	double par_sigma;
	double par_mu;
};

bool operator<(const Variational &a, const Variational &b){
	return a.en_best<b.en_best; 
}

vector<Variational> Var;

// functions
double error(double AV, double AV2, int n);
double distribution(double _x, double _sigma, double _mu);
double pot(double _x);
double kin(double _x, double _sigma, double _mu);
void input(void);

// parameters
int M, N, L, conta, nbins, prog, q;
double x_0, x, sigma, mu, _min, _max, b, g, alpha, d, x_new, x_min, x_max, bin_size, mu_min, mu_max, sigma_min, sigma_max;


/****************************************************************

			ANDREA MEDAGLIA 
			   LEZIONE 8
			 ESERCIZIO 8.1

*****************************************************************/
