#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

typedef vector<vector<int> > v2d;
double bmann_table[3] = {0,0,0};

void init_lattice(v2d &lattice, int dim);
void mc_step(v2d &lattice, double N);

bool flip_spin(v2d &lattice);
double get_energy(v2d &lattice);
double get_magnetization(v2d &lattice);
void update_bmann_table(double beta);
double bmann_weigth(int deltaE);

void output_config(v2d &lattice, int t);
void output_log(v2d &lattice, int t, double T);

int main() {

	// Lattice geometry
	int dim = 100;
	int N = dim*dim;
	
	// Start sampling after tterm mc_steps, run tmax mc_steps per Markov chain
	int tmax = 10000;
	int tterm = 3000;

	// Explore temperature range [T0,T1] using ntemp bins
	double T0 = 1.;
	double T1 = 5.;
	int ntemp = 100;
	double deltaT = (T1-T0)/ntemp;

	// Output file, append content to existing content.
	// Format: # T    <E>    <E^2>    <m>    <m^2>    N
	ofstream of;
	of.open("data.txt", ofstream::app);

	// Seed PRG.
	srand(time(NULL));
	
	// Run a MC for every temperature
	for (double T = T0; T < T1; T += deltaT)
	{
		cout << "Start MC sim T=" << T << endl;
		clock_t begin = clock();

		// Precompute positive deltaE Boltzmann weigths for given T(beta): exp(-beta*deltaE)
		double beta = 1./T;
		update_bmann_table(beta);

		// Declare and random initialize a spin lattice dim*dim. 
		v2d spin;
		init_lattice(spin, dim);

		// Energy and Magnetization accumulators
		long long AE = 0, AE2 = 0, AM = 0, AM2 = 0, AN = 0;

		// Start MC
		for (int t = 0; t < tmax; ++t)
		{
			// Update N lattice loci.
			// Actually PROPOSE N lattice loci update.
			mc_step(spin, N);

			// Run sampling
			if (t >= tterm)
			{
				int E = get_energy(spin);
				int M = get_magnetization(spin);
				AE += E;
				AE2 += E*E;
				AN++;
				AM += M;
				AM2 += M*M;
			}

			/*
			*	Some logs and check data
			*/
			if (T == T0 && (t == 0 || t == 10 || t == 100 || t == 1000 || t == 9999))
			{
				output_config(spin,t);
			}
			if (t % 100 == 0)
			{
				output_log(spin,t,T);
			}
		}

		cout << "Computation time " << float(clock() - begin)/CLOCKS_PER_SEC << " s" << endl;
		// Storing data
		of << T << "\t" << dim << "\t" << AE/((double)AN)/N << "\t" << AE2/((double)AN)/N/N << "\t" << AM/((double)AN)/N << "\t" << AM2/((double)AN)/N/N << "\t" << AN << endl; 

	}

	of.close();
	return 0;
}


void output_config(v2d &lattice, int t)
{
	ofstream cf;
	cf.open("config_"+to_string(t)+".txt");

	int dim = lattice.size();
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			cf << (lattice[i][j]>0);
		}
		cf << endl;
	}
	cf.close();
}

void output_log(v2d &lattice, int t, double T)
{
	ofstream lf;
	lf.open("logs.txt",ofstream::app);

	int dim = lattice.size();
	int sum = 0;

	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			sum += (lattice[i][j] > 0);
	
	lf << T << "\t" << dim << "\t" << (double)sum/(dim*dim) << "\t" << t << endl;
}

void init_lattice(v2d &lattice, int dim)
{
	lattice.clear();

	for (int i = 0; i < dim; ++i)
	{
		vector<int> tmp;

		for (int j = 0; j < dim; ++j)
		{
			int r = -1. + ((double)rand()/RAND_MAX > 0.5)*2.;
			tmp.push_back(r);
		}

		lattice.push_back(tmp);
		tmp.clear();
	}
}

bool flip_spin(v2d &lattice)
{
	int dim = lattice.size();
	int N = dim*dim;

	int index = (double)rand()/RAND_MAX*(N-1);
	
	int x = index/dim;
	int y = index%dim;

	int sum = lattice[(x+1)%dim][y] + lattice[(x+dim-1)%dim][y] + lattice[x][(y+1)%dim] + lattice[x][(y+dim-1)%dim];
	int deltaE = 2*lattice[x][y]*sum;

	if (deltaE <= 0)
	{
		lattice[x][y] *= -1;
		return true;
	}
	else
	{
		double r = (double)rand()/RAND_MAX;
		if (r <= bmann_weigth(deltaE))
		{
			lattice[x][y] *= -1;
			return true;
		}

		// microstep rejected
		return false;
	}
}

double get_energy(v2d &lattice)
{
	double energy = 0.;
	int dim = lattice.size();

	for (int x = 0; x < dim; ++x)
	{
		for (int y = 0; y < dim; ++y)
		{
			int sum = lattice[(x+1)%dim][y] + lattice[x][(y+1)%dim];
			energy += -lattice[x][y]*sum;
		}
	}

	return energy;
}

double get_magnetization(v2d &lattice)
{
	double m = 0.;
	int dim = lattice.size();

	for (int x = 0; x < dim; ++x)
		for (int y = 0; y < dim; ++y)
			m += lattice[x][y];

	return m;
}

void mc_step(v2d &lattice, double N)
{
	// MC step = N microsteps (spin flips)
	int microsteps = 0;
	while (microsteps < N)
	{
		if (flip_spin(lattice))
			microsteps++;
		else
			microsteps++;
	}
}

void update_bmann_table(double beta)
{
	// possible increasing energy: +8 +4
	// 0 added for consistency

	for (int i = 0; i < 3; ++i)
	{
		bmann_table[i] = exp(-4.*beta*i);
	}
}
double bmann_weigth(int deltaE)
{
	int i = deltaE/4;
	return bmann_table[i];
}
