#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm> //max_element

using namespace std;

typedef vector<vector<int> > v2d;

void init_lattice(v2d &lattice, int dim, double p);
void set_lattice(v2d &lattice, int a);
void print_lattice(v2d &lattice, const char *filename);

int find_real(vector<int> &clusters, int x);
void link(vector<int> &clusters, int x, int y);

vector<int> clusterize_lattice(v2d &lattice);
void reduce_lattice(v2d &lattice, vector<int> &clusters);
void clusters_length(v2d &lattice, vector<int> &clusters);
int mega_cluster(v2d &lattice, vector<int> &clusters);

void sim(int dim, double p);



void sim(int dim, double p)
{	
	int debug = 0;

	v2d lattice;
	vector<int> clusters;

	// Initialize every entry of lattice randomly
	init_lattice(lattice,dim,p);
	if (debug) cout << "Init" << endl;
	// Get clusters rules (primitives and substitution chains)
	clusters = clusterize_lattice(lattice);
	if (debug) cout << "Clustered" << endl;
	// Apply clusters rules (i.e. substitute non-primitive cluster ids)
	reduce_lattice(lattice,clusters); 
	if (debug) cout << "Reduced" << endl;

	// Some log to check
	print_lattice(lattice,"lattice.txt");
	if (debug) cout << "Printed" << endl;

	// Updates clusters so that entry [x] contains number of members of cluster labeled x
	clusters_length(lattice,clusters);
	if (debug) cout << "Lengths computed" << endl;
	
	int M = mega_cluster(lattice, clusters);
	if (debug) cout << "Mega detected" << endl;

	// Time to gather simulation results
	// Length of longest cluster
	int max_s = 0;
	for (int i = 0; i < clusters.size(); ++i)
	{
		if (clusters[i] > max_s)
			max_s = clusters[i];
	}
	if (debug) cout << "MAXS detected" << endl;

	// n[s] : number of clusters having length s. NOTICE: there will be several clusters having 0 length (due to Holshen-Kopelman alg implementation)
	vector<int> count(max_s+1,0);
	for (int i = 0; i < clusters.size(); ++i)
		count[clusters[i]] += 1;
	if (debug) cout << "count[s] computed" << endl;

	// numbers of clusters having length > 0:
	int N = 0;
	// chi = SUM{ s^2 * n[s] }
	double chi_p = 0., chi = 0.;
	
	for (int s = 1; s < count.size(); ++s)
	{
		N += count[s];
		chi += s*s*count[s];
	}

	// empty lattice
	if (N == 0)
		return;

	chi_p = chi - max_s*max_s*count[max_s];
	chi = chi/N;
	chi_p = chi_p/N;

	if (debug) cout << "N,chi,chi_p computed" << endl;

	// Output sim
	ofstream of;
	of.open("sims.txt",ofstream::app);
	of << dim << "\t" << p << "\t" << max_s << "\t" << chi << "\t" << chi_p << "\t" << M << endl;
	of.close();
	
}


int main()
{
	double p = 0.5;
	int dim = 10;
	int imax = 100000;
	srand(time(NULL));

	for (int i = 0; i < imax; ++i)
	{
		p = rand()/(double)RAND_MAX;
		if (i%1000 == 0) cout << i << "/" << imax << endl;
		sim(dim,p);	
	}
	
	return 0;
}



/*
*	Initialize a dim*dim lattice and set every entry to occupied (1) with probability p
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

/*
*	Returns primitive clusters and updates lattice (-1/1) to labels
*/
vector<int> clusterize_lattice(v2d &lattice)
{
	
	v2d label(lattice);		set_lattice(label,-1);
	vector<int> clusters; 	clusters.clear();

	int dim = lattice.size();

	int last_id = 0;
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			if (lattice[i][j] > -1)
			{
				int l = -1;

				int up = label[(i+dim-1)%dim][j];
				int left = label[i][(j+dim-1)%dim];
				
				if (up >= 0 && left >= 0)
				{
					l = find_real(clusters,min(up,left));
					link(clusters,up,left);
				} 
				else if(up >= 0)
				{
					l = find_real(clusters,up);
				} 
				else if(left >= 0)
				{
					l = find_real(clusters,left);
				}
				else
				{
					l = last_id;
					clusters.push_back(last_id);
					last_id++;
					
				}

				label[i][j] = l;
			}
		}
	}

	lattice = label;

	return clusters;
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


void clusters_length(v2d &lattice, vector<int> &clusters)
{
	int dim = lattice.size();
	// set clusters to 0: this vector will be used to count length of every cluster
	for (int i = 0; i < clusters.size(); ++i)	
		clusters[i] = 0;
	
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{	
			int l = lattice[i][j];
			if (l != -1)
				clusters[l] += 1;
		}
	}
}

int mega_cluster(v2d &lattice, vector<int> &clusters)
{
	int dim = lattice.size();
	// first row vs last row
	for (int j = 0; j < dim; ++j)
	{
		for (int jj = 0; jj < dim; ++jj)
		{	
			// check if there is a collision between first row and last row OR first col and first col
			if (((lattice[0][j] != -1) && (lattice[0][j] == lattice[dim-1][jj])) ||	((lattice[j][0] != -1) && (lattice[j][0] == lattice[jj][dim-1])))
			{
				//cout << "Collision! (j,jj)" << j <<" "<< jj << endl;
				return 1;
			}
		}
	}

	// if you get here we found no collision, so there's no megacluster_bombXX432^!£*
	return 0;
}
