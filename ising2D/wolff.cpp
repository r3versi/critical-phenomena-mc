#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

typedef vector<vector<int> > v2d;
int debug = 1;	// debug flag

double pb = 0.;	// 1 - exp(-2J/T)  probability of bond forming (updated every time wolff_step is called)
void update_pb(double T);	

void init_lattice(v2d &lattice, int dim, double p);
void print_lattice(v2d &lattice, const char *filename);

void wolff_step(v2d &lattice, double T);
void flip_cluster(v2d &lattice, v2d &cluster);
void widen_cluster(v2d &lattice, v2d &cluster, int x, int y);
void try_bond(v2d &lattice, v2d &cluster, int x, int y);



int main()
{
	double T = 2.;	// Temperature
	double p = 0.5;	// Probability of init spin up
	int dim = 10;	// Lattice side

	int mterm = 100;	// MC termalization steps
	int mmax = 101;	// MC max steps

	srand(time(NULL));

	v2d lattice;
	init_lattice(lattice, dim, p);

	if (debug)	print_lattice(lattice, "start.txt");
	
	for (int m = 0; m < mmax; ++m)
	{
		wolff_step(lattice,T);
		
		if (m >= mterm)
		{	
			if (debug)	cout << "Sampling" << endl;
		}	
	
	}

	if (debug)	print_lattice(lattice,"end.txt");
	
	return 0;
}

/*
*	Wolff step:
*	1. update probability of bond forming
*	2. start clusterizing from a random site
*	3. flip cluster
*/
void wolff_step(v2d &lattice, double T)
{
	// geometry
	int dim = lattice.size();
	int N = dim*dim;

	// update probability of bond forming
	update_pb(T);

	// entry = 0/1, site belongs to cluster
	v2d cluster(dim, vector<int>(dim,0));

	// get random starting site
	int index = rand()/(double)RAND_MAX*(N-1);
	int x = index/dim;
	int y = index%dim;

	// recursively add sites to cluster
	widen_cluster(lattice,cluster,x,y);

	// flip cluster
	flip_cluster(lattice,cluster);

}


/*
*	Explore (x,y) neighbours and try to add them to the cluster.
*	Function called recursively by try_bond every time a new site
*	is added to the cluster.
*/
void widen_cluster(v2d &lattice, v2d &cluster, int x, int y)
{
	int dim = lattice.size();
	int s = lattice[x][y];

	// Add (x,y) to cluster
	cluster[x][y] = 1;
	
	//neighbours (periodic conditions)
	int x_up = (x+dim-1)%dim;
	int x_down = (x+1)%dim;
	int y_left = (y+dim-1)%dim;
	int y_right = (y+1)%dim; 

	// ferromagnetic interaction
	if (s*lattice[x_up][y] > 0 && cluster[x_up][y] == 0)
		try_bond(lattice,cluster,x_up,y);
	
	if (s*lattice[x_down][y] > 0 && cluster[x_down][y] == 0)
		try_bond(lattice,cluster,x_down,y);
	
	if (s*lattice[x][y_left] > 0 && cluster[x][y_left] == 0)
		try_bond(lattice,cluster,x,y_left);

	if (s*lattice[x][y_right] > 0 && cluster[x][y_right] == 0)
		try_bond(lattice,cluster,x,y_right);

}

/*
*	If bond successfully formed, call widen_cluster starting from this site.
*	Cluster is build recursively.
*/
void try_bond(v2d &lattice, v2d &cluster, int x, int y)
{	
	double r = rand()/(double)RAND_MAX;
	if (r <= pb)
		widen_cluster(lattice, cluster, x, y);
}

/*
*	Flip whole cluster
*/
void flip_cluster(v2d &lattice, v2d &cluster)
{
	int dim = lattice.size();

	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			if (cluster[i][j])
				lattice[i][j] *= -1;
}



/*
*	Update globally defined pb
*/
void update_pb(double T)
{
	pb = 1. - exp(-2./T);
}

/*
*	Initialize a dim*dim lattice and set every entry to up with probability p
*/
void init_lattice(v2d &lattice, int dim, double p)
{
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
