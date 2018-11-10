#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>

int debug = 1;

using namespace std;

typedef vector<vector<int> > v2d;

double randd(); // get uniform rand double between 0 and 1

void init_lattice(v2d &lattice, int dim, double p);
void set_lattice(v2d &lattice, int a);
void print_lattice(v2d &lattice, const char *filename);

int find_real(vector<int> &clusters, int x);
void link(vector<int> &clusters, int x, int y);

double get_energy(v2d &lattice);
double get_magnetization(v2d &lattice);

v2d lattice_bonds(v2d &lattice, double T);
void reduce_lattice(v2d &lattice, vector<int> &clusters);

void mc_step(v2d &lattice, double T);



/*
*	Run one MC step:
*	1. make bonds (= clusterize) with probability 1-e^(-2/T)
*	2. run on every cluster and flip it with probability = 1/2
*/
void mc_step(v2d &lattice, double T)
{
	if (debug)
		print_lattice(lattice,"spins.txt");

	v2d label(lattice);		set_lattice(label,-1);
	vector<int> clusters;	clusters.clear();

	int dim = lattice.size();

	// Precompute probability of forming a bond between aligned sites
	double b_prob = 1. - exp(-2./T);

	// Clusterize lattice: form a bond with prob. b_prob
	int last_id = 0;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			int spin = lattice[i][j];
			int up = lattice[(i+dim-1)%dim][j];
			int left = lattice[i][(j+dim-1)%dim];
			
			int my_label = label[i][j];
			int label_up = label[(i+dim-1)%dim][j];
			int label_left = label[i][(j+dim-1)%dim];


			if (my_label == -1)
			{
				label[i][j] = last_id;
				my_label = last_id;
				clusters.push_back(last_id);
				last_id += 1;
			}

			if (spin*up > 0)
			{
				if (randd() <= b_prob)
				{
					if (label_up == -1)
					{
						label[(i+dim-1)%dim][j] = my_label; 
					}
					else
					{
						link(clusters,my_label,label_up);
					}
				}
			}

			if (spin*left > 0)
			{
				if (randd() <= b_prob)
				{
					if (label_left == -1)
					{
						label[i][(j+dim-1)%dim] = my_label; 
					}
					else
					{
						link(clusters,my_label,label_left);
					}
				}
			}
		}
	}

	reduce_lattice(label,clusters);
	
	// Reduce clusters
	int k = 0;
	for (int i = 0; i < clusters.size(); ++i)
	{
		if (clusters[i] == i)
		{
			clusters[i] = k;
			k += 1;
		}
	}

	// number of unique clusters
	int NC = k;
	// normalizing label and clusters ids from (RANDOM,RANDOM) to (0,Nc-1)	
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			label[i][j] = clusters[label[i][j]];
		}
	}

	// clusters[i] = 1 => flip spins in cluster i
	clusters.clear();
	for (int i = 0; i < NC; ++i)
		clusters.push_back(randd() <= 0.5);

	// Flip clusters
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			if (clusters[label[i][j]])
				lattice[i][j] *= -1;
	
	if (debug)
		print_lattice(lattice,"flipped.txt");
	
}

int main()
{	
	// Lattice geometry
	int dim = 100;
	int N = dim*dim;

	// probability of spin UP
	double p = 0.5; 

	// Explore temperature range [T0,T1] using ntemp bins
	double T0 = 1.;
	double T1 = 5.;
	int ntemp = 100;
	double deltaT = (T1-T0)/ntemp;

	// Start sampling after iterm mc_steps, run imax steps per MC
	int imax = 1000;
	int iterm = 100;

	// Spin lattice
	v2d spin;
	
	// Seed PRG.
	srand(time(NULL));

	for (double T = T0; T < T1; T += deltaT)
	{
		if (debug)
			cout << "Start MC sim T=" << T << endl;
		
		// Timing
		clock_t begin = clock();

		// Reset Energy and Magnetization accumulators
		long long AE = 0, AE2 = 0, AE4 = 0, AM = 0, AM2 = 0, AM4 = 0, AN = 0;

		// Reset lattice to random
		init_lattice(spin,dim,p);

		for (int i = 0; i < imax; ++i)
		{	
			mc_step(spin,T);

			if (i >= iterm)
			{
				int E = get_energy(spin);
				int M = get_magnetization(spin);
				int E2 = E*E;
				int M2 = M*M;
				AE += E;
				AE2 += E2;
				AE4 += E2*E2;
				AM += M;
				AM2 += M2;
				AM4 += M2*M2;
				AN++;
			}
		}
		
		if (debug)
			cout << "Computation time " << float(clock() - begin)/CLOCKS_PER_SEC << " s" << endl;
		
		ofstream of;
		of.open("sims.txt",ofstream::app);
		of << T << "\t" << dim << "\t" << AE/((double)AN)/N << "\t" << AE2/((double)AN)/N/N << "\t" << AM/((double)AN)/N << "\t" << AM2/((double)AN)/N/N << "\t" << AM4/((double)AN)/N/N/N/N << "\t" << AN << endl; 
		of.close();
	}	
	
	return 0;
}

/**FROM PERCOLATION**/
/*
*	Initialize a dim*dim lattice and set every entry to occupied (1) with probability p
*/
void init_lattice(v2d &lattice, int dim, double p)
{
	lattice.clear();
	vector<int> tmp;
	for (int i = 0; i < dim; ++i)
	{
		tmp.clear();

		for (int j = 0; j < dim; ++j)
		{
			tmp.push_back(-1+(rand()/(double)RAND_MAX < p)*2);
		}
		
		lattice.push_back(tmp);
	}
	//cout << "Lattice initialized" << endl;
}

/*
*	Set every entry of lattice to int a
*/
void set_lattice(v2d &lattice, int a)
{	
	int dim = lattice.size();
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			lattice[i][j] = a;
}

void print_lattice(v2d &lattice, const char *filename)
{
	ofstream of;
	of.open(filename);

	for (int i = 0; i < lattice.size(); ++i)
	{
		for (int j = 0; j < lattice.size(); ++j)
		{
			if (lattice[i][j]<0)
				of << "x";
			else
				of << lattice[i][j];
			of << "\t";
		}
		of << endl;
	}
	of.close();
}

/*
*	Substitutes dead clusters id according to clusters content.
*/
void reduce_lattice(v2d &lattice, vector<int> &clusters)
{
	int dim = lattice.size();

	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			if (lattice[i][j] > -1)
			{
				int l = lattice[i][j];

				if (clusters[l] == l)
					continue;
				
				lattice[i][j] = find_real(clusters,l);
			}
		}
	}
}


/*
*	Loop into clusters[] to find primitive ID of cluster x
*	and (2nd loop) updates every node of the chain to that ID
*	so that next lookup will be faster
*/
int find_real(vector<int> &clusters, int x)
{
	// Search for primitive of cluster labelled x
	int y = x;
	while (clusters[y] != y)
		y = clusters[y];

	// Update every node of the chain from x to primitive of x
	while (clusters[x] != x)
	{
		int z = clusters[x];
		clusters[x] = y;
		x = z;
	}
	return y;
}
/*
*	Execute a collision between cluster labelled x and cluster labelled y.
*	Look for primitive of x and y, and set every node of both chains to minimum
*	primitive
*/
void link(vector<int> &clusters, int x, int y)
{
	if (x < y)
		swap(x,y);

	clusters[find_real(clusters,x)] = find_real(clusters,y);
	//cout << "Linked " << x << " to " << y << endl;
}

/**FROM ISING**/
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

double randd()
{
	return rand()/((double)RAND_MAX);
}
