#pragma GCC optimize "O3,omit-frame-pointer,inline"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

void update_weights();
bool metropolis_step();
bool can_move(int i, int j);
double get_prob_ratio(int i, int j);
double get_white_prob(int i, int j);
double get_black_prob(int i, int j);
void apply_move(int i, int j);
double get_mz();
double get_energy();
void print_conf(string filename, bool special = false);

vector<vector<double> > weights = {{0,0},{0,0}};
vector<vector<double> > ai_matrix = {{0,0},{0,0}};

bool DEBUG = true;

int L = 10;	// L temporal slices, Suzuki
int N = 10; // N - spin chain
double T = 1.;
int NTmax = 500;
double Jx = -10.;
double Jz = 0.;
double tau = 1. / (L*T);

int tterm = 3000;
int Mmax = 100000; 

vector<vector<int> > grid;

int main()
{
	string output = "";

	srand(time(NULL));
	
	int NT = 0;
	while (NT < NTmax)
	{
		NT++;
		T = 1.5*rand()/(double)RAND_MAX;
		L = round(10/T);
		L += L%2;
		tau = 1./(L*T);

		cout << "Init sim: " << T << endl;
		cout << "L: " << L << " tau: " << tau << endl; 	
		// FIXED tau!

		// precompute table of weights
		update_weights();

		// grid initialization L*N	
		grid.clear();
		vector<int> tmp; tmp.clear();
		int m = 0.;
		for (int i = 0; i < N; ++i)
		{
			int spin = -1 + (rand()/(double)RAND_MAX > .5)*2;
			tmp.push_back(spin);
			m += spin;
		}

		if (m == N || m == -N)
		{
			NT--;
			continue;
		}

		for (int i = 0; i < L; ++i)
			grid.push_back(tmp);

		if (DEBUG)
		{
			cout << "Init m = " << m << endl;
			print_conf("logs/0.conf");
		}
		

		double AMZ = 0., AE = 0., AEE = 0.;
		int t = 0, n = 0, M = 0;
		
		while (M < Mmax)
		{
			bool update = false;

			/*if (t % 10000 == 0 && DEBUG)
				print_conf("logs/"+to_string(t)+".conf");
			*/

			if (metropolis_step())
			{
				t++;
				update = true;
			}
			else
				n++;

			if (t < 500 && DEBUG && update && NT == 1)
				print_conf("logs/"+to_string(t)+".conf");


			if (t >= tterm && t%(N*L) == 0)
			{
				M++;
				double mz = get_mz();
				double en = get_energy();
				AMZ += mz;
				AE += en;
				AEE += en*en;
			}
		}

		if (DEBUG)
		{
			cout << "Computing terminated" << endl;
			cout << "Metropolis steps, accepted: " << t << " (" << 100*t/((double)(t+n)) << "%), rejected: " << n << " (" << 100*n/((double)(t+n)) << "%)" << endl;
			cout << "Termalization steps (Sux) " << tterm << endl;
			cout << M << " measures" << endl;
		}	
		cout << NT << "\t<E> = " << AE/((double)M) << endl;
		
		double E = AE/(double)M;
		double E2 = AEE/(double)M;
		double Cv = (E2-E*E)/(T*T);

		ofstream of;
		of.open("results/taufixed_N"+to_string(N)+".dat",ofstream::app);
		//of.open("results/L"+to_string(L)+"_N"+to_string(N)+".dat",ofstream::app);
		//of << T << "\t" << E << "\t" << (E2-E*E)/T << "\t" << m << endl;
		of << L << "\t" << tau << "\t" << T << "\t" << m << "\t" << E << "\t" << Cv << endl;
		of.close();

	}	
	
	return 0;
}

void update_weights()
{

	tau = 1./(L*T);

	weights[0][0] = -sinh(.5*Jx*tau) * exp(Jz*tau);
	weights[0][1] = cosh(.5*Jx*tau) * exp(Jz*tau);
	weights[1][0] = 0.;
	weights[1][1] = exp(-Jz*tau);

	ai_matrix[0][0] = -.5*Jx/tanh(.5*Jx*tau);
	ai_matrix[0][1] = -.5*Jx*tanh(.5*Jx*tau);
	ai_matrix[1][0] = 0.;
	ai_matrix[1][1] = 0.;
}

bool metropolis_step()
{
	int i = round(rand()/(double)RAND_MAX*(L-1));
	int j = round(rand()/(double)RAND_MAX*(N-1));

	if (!can_move(i,j))
		return false;

	double prob_ratio = get_prob_ratio(i,j);
	double r = rand()/(double)RAND_MAX;
	if (r >= prob_ratio)
	{
		// revert move applied when computing get_prob_ratio
		apply_move(i,j);
		return false;
	}
	// accept move
	return true;
}

bool can_move(int i, int j)
{

	if (grid[i][j] < 0 || grid[(i+1)%L][j] < 0) 
		return false;

	if (i%2 == j%2)
	{
		if (grid[i][(j-1+N)%N] > 0 || grid[(i+1)%L][(j-1+N)%N] > 0)
			return false;
	}
	else
	{
		if (grid[i][(j+1)%N] > 0 || grid[(i+1)%L][(j+1)%N] > 0)
			return false;
	}

	return true;
}

double get_prob_ratio(int i, int j)
{
	if (i%2 == j%2)
	{
		double old_p = get_white_prob(i, (j-1+N)%N);
		apply_move(i,j);
		double new_p = get_white_prob(i, (j-1+N)%N);

		return min(1.,new_p/old_p);
	}
	else
	{

		double old_p = get_white_prob(i, j);
		apply_move(i,j);
		double new_p = get_white_prob(i, j);

		return min(1.,new_p/old_p);
	}
}

// where (i,j) = bottom-left vertex (of white space)
double get_white_prob(int i, int j)
{
	double p1 = get_black_prob(i, (j-1+N)%N);
	double p2 = get_black_prob((i+1)%L, j);
	double p3 = get_black_prob(i,(j+1)%N);
	double p4 = get_black_prob((i-1+L)%L,j);

	return p1*p2*p3*p4;
}

// where (i,j) = bottom-left vertex (of a black space)
double get_black_prob(int i, int j)
{
	int Ph = (grid[i][j] * grid[i][(j+1)%N]) > 0;
	int Pv = (grid[i][j] * grid[(i+1)%L][j]) > 0;

	return weights[Ph][Pv];
}

// if left == true, move = top-left, else, move = top-right
void apply_move(int i, int j)
{
	if (i%2 == j%2)
	{
		grid[i][(j+N-1)%N] *= -1;
		grid[(i+1)%L][(j+N-1)%N] *= -1;
	}
	else
	{
		grid[i][(j+1)%N] *= -1;
		grid[(i+1)%L][(j+1)%N] *= -1;
	}

	grid[i][j] *= -1;
	grid[(i+1)%L][j] *= -1;
}

double get_mz()
{
	double mz = 0.;

	for (int i = 0; i < L; ++i)
	{
		int s = 0;
		for (int j = 0; j < N; ++j)
			s += grid[i][j];

		mz += s;
	}

	return mz/(double)L;
}

/*	Assuming only black boxes */
double get_energy()
{
	double AE = 0.;

	for (int i = 0; i < L; ++i)
	{	

		// loop only on black boxes
		for (int j = i%2; j < N; j += 2)
		{
			int Ph = (grid[i][j]*grid[i][(j+1)%N]) > 0;
			int Pv = (grid[i][j]*grid[(i+1)%L][j]) > 0;

			AE += ai_matrix[Ph][Pv] + Jz*grid[i][j]*grid[i][(j+1)%N];
		}
	}

	return AE/(double)L;
}

void print_conf(string filename, bool special)
{

	ofstream of;
	of.open(filename);
	if (special)
	{	
		for (int i = L-1; i >= 0; i--)
		{
			for (int j = 0; j < N; ++j)
			{
				of << (int)(grid[i][j] > 0);	
			}
			of << endl;
		}
	}
	else
	{
		for (int i = 0; i < L; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				of << (int)(grid[i][j] > 0) << " ";	
			}
			of << endl;
		}	
	}
	of.close();

}