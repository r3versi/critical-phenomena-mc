#pragma GCC optimize "O3,omit-frame-pointer,inline"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

enum DEBUG {VERBOSE,LIGHT,NONE};

/*
*	Grid base class: periodic boundaries conditions handled.
*/
class Grid 
{
public: 
	unsigned int N, L;
	vector<vector<int> > grid;

	Grid(){}
	Grid(const Grid& grid) = default;
	Grid(int L, int N)
	{
		this->N = N;
		this->L = L;
	}

	int get(int i, int j) {return grid[(i+L)%L][(j+N)%N];}
	void set(int i, int j, int value) {grid[(i+L)%L][(j+N)%N] = value;}
	void init(int value) {grid = vector<vector<int> >(L,vector<int>(N,value));}
};

/*
*	Spin lattice.
*/
class Spins: public Grid
{
public:
	Spins(){}
	Spins(const Spins& spins) = default;
	Spins(int L, int N)
	{
		this->N = N;
		this->L = L;
	}

	/* Initialize lattice picking randomly up/down spins along time slice 0 and replicating it L times. */
	int random_init()
	{
		grid.clear();
		vector<int> tmp; tmp.clear();

		int m = 0;
		for (int i = 0; i < N; ++i)
		{
			int spin = (rand()/(double)RAND_MAX > .5)*2 -1;
			tmp.push_back(spin);
			m += spin;
		}

		for (int j = 0; j < L; ++j)
			grid.push_back(tmp);
		
		return m;
	}

	/* Flips spin at site (i,j) */
	void flip(int i, int j)
	{
		grid[(i+L)%L][(j+N)%N] *= -1;
	}

	/* Get z-magnetization along chain at time slice 0 */
	int get_sz()
	{
		int sz = 0;
		if (grid.size() == 0)
			return 0;

		for (int i = 0; i < N; ++i)
			sz += grid[0][i];

		return sz;
	}

	/* Get index of spins configuration on a black plaquette (4 spins). (i,j) bottom-left site. */
	int get_u(int i, int j)
	{
		int Ph = (get(i,j+1) * get(i,j)) > 0;
		int Pv = (get(i+1,j) * get(i,j)) > 0;
		return Ph+Pv;
	}
};

/*
*	Loops class containing Holshen-Kopelmann clusterization algorithm.
*/
class Loops: public Grid
{
public:
	vector<int> loops;

	Loops(){}
	Loops(const Loops& loops) = default;
	Loops(int L, int N)
	{
		this->N = N;
		this->L = L;
		
		loops.clear();
		grid = vector<vector<int> >(L,vector<int>(N,-1));
	}

	/* Return primitive loop ID whose x belongs to. */
	int find_real(int x)
	{
		int y = x;
		while(loops[y] != y)
			y = loops[y];
		
		while(loops[x] != x)
		{
			int z = loops[x];
			loops[x] = y;
			x = z;
		}

		return y;
	}

	/* Links loops ID x and y. */
	void link(int x, int y)
	{
		if (x < y)
			swap(x,y);

		loops[find_real(x)] = find_real(y);
	}

	/* Clusterize applying node type g to black plaquette (i,j). */
	void apply_node(int i, int j, int g)
	{
		/*	
			g=0: vertical lines
			g=1: horizontal lines
			g=2: crossed lines		
		*/
		int next_id = loops.size();

		if (g == 0)
		{
			if (get(i,j) == -1)
			{

				loops.push_back(next_id);
				set(i,j, next_id);
				set(i+1,j, next_id);
				next_id++;

				loops.push_back(next_id);
				set(i,j+1, next_id);
				set(i+1,j+1, next_id);
			
			}
			else if (get(i+1,j) == -1)
			{
				set(i+1,j, get(i,j));
				set(i+1,j+1, get(i,j+1));
			}
			else
			{
				link(get(i,j), get(i+1,j));
				link(get(i,j+1), get(i+1,j+1));
			}
		}
		else if(g == 1)
		{
			if (get(i,j) == -1)
			{
				loops.push_back(next_id);
				set(i,j, next_id);
				set(i,j+1, next_id);
				next_id++;

				loops.push_back(next_id);
				set(i+1,j, next_id);
				set(i+1,j+1, next_id);
			}
			else if (get(i+1,j) == -1)
			{
				link(get(i,j), get(i,j+1));

				loops.push_back(next_id);
				set(i+1,j, next_id);
				set(i+1,j+1, next_id);
			}
			else
			{
				link(get(i,j), get(i,j+1));
				link(get(i+1,j), get(i+1,j+1));
			}
		}
		else if (g == 2)
		{
			if (get(i,j) == -1)
			{
				loops.push_back(next_id);
				set(i,j, next_id);
				set(i+1,j+1, next_id);
				next_id++;

				loops.push_back(next_id);
				set(i,j+1, next_id);
				set(i+1,j, next_id);
			}
			else if (get(i+1,j) == -1)
			{
				set(i+1,j+1, get(i,j));
				set(i+1,j, get(i,j+1));
			}
			else
			{
				link(get(i+1,j+1), get(i,j));
				link(get(i+1,j), get(i,j+1));	
			}
		}
	}

	/* Final loop IDs substitution. */
	void reduce_loops()
	{

		for (int i = 0; i < loops.size(); ++i)
			loops[i] = find_real(i);

		for (int i = 0; i < L; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				if (get(i,j) == -1)
					continue;
				
				set(i,j, loops[get(i,j)]);
			}
		}
	}

};

class Sim
{
private:
	vector<double> conf;
	vector<double> node_weight;
	vector<vector<int>> DUG;
	vector<vector<double>> PUG;

	vector<double> Hnd_matrix;
	
	bool calc_acorr;
	DEBUG debug_level;

	double T;
	int L, N;

	double Jx, Jz, D, tau;

	int term, measures;

	Spins lattice;
	Grid nodes;
	Loops labels;

	int t, M;
	double AE, AEE;

	int tmax; // max correlation time to test
	vector<double> E_hist;
	vector<double> acorr;
	vector<int> counts;
public:

	Sim(int L, int N, double T);

	void run();
	void step();

	void update_corr(double E);

	void compute_weights();
	int get_random_g(int u);
	double get_energy();

	void set_Jx(double Jx);
	void set_Jz(double Jz);
	
	void set_termalization(int term);
	void set_nmeasures(int measures);
	
	void set_debug(DEBUG debug_level);
	void calc_autocorr(bool flag);
	void set_acorr_time(int steps);

	void debug_output(string filename);
	void save_output();
};

Sim::Sim(int L, int N, double T)
{
	if (L%2 != 0)
		cout << "ERROR: odd L." << endl;
	if (N < 0)
		cout << "ERROR: N negative." << endl;
	if (T < 0.)
		cout << "ERROR: T negative." << endl;

	this->L = L;
	this->N = N;
	this-> T = T;

	// Set defaults for termalization, measure gathering and autocorrelation.
	set_termalization(100);
	set_nmeasures(10000);
	set_acorr_time(1000);
	calc_acorr = false;

	// Physical default constants
	Jx = -10.;
	Jz = 0.;
	D = Jz/Jx;
	tau = 1./(T*L);

	// weights of spin configurations
	conf = {0,0,0};

	// energy non-diagonal part matrix elements
	Hnd_matrix = {0,0,0};

	// weights of nodes
	node_weight = {0,0,0};

	// compatibility matrix (spin/node)
	DUG = {
		{0,1,1},
		{1,1,0},
		{1,0,1}};

	// g weights normalized to 1
	PUG = {
		{0,0,0},
		{0,0,0},
		{0,0,0}};

	lattice = Spins(L,N);
	nodes = Grid(L,N);
	labels = Loops(L,N);

	lattice.random_init();
	nodes.init(-1);
}

void Sim::update_corr(double E)
{
	if (E_hist.size() < tmax)
		E_hist.push_back(E);
	else
	{
		acorr[0] += E*E;
		counts[0] ++;
		
		for (int i = 1; i < tmax; ++i)
		{
			acorr[i] += E*E_hist[tmax-i];
			counts[i] ++;
		}

		E_hist.erase(E_hist.begin());
		E_hist.push_back(E);
	}
}

void Sim::run()
{
	compute_weights();

	t = 0;
	M = 0;
	AE = 0.;
	AEE = 0.;

	while(M < measures)
	{
		if (debug_level == VERBOSE)
			debug_output(to_string(t));

		// update lattice applying global moves
		step();
		t++;


		// if lattice is termalized, get measures
		if (t >= term)
		{
			M++;
			double en = get_energy();
			AE += en;
			AEE += en*en;

			if (calc_acorr)
				update_corr(en);
		}	
	}
	
	if (debug_level != NONE)
		cout << "Simulation completed." << endl;
}

void Sim::save_output()
{
	double E = AE/(double)M;
	double E2 = AEE/(double)M;
	double Cv = (E2-E*E)/(T*T);
	string filename = "L"+to_string(L)+"_N"+to_string(N)+".dat";

	if (debug_level != NONE)
	{
		cout << "N=" << N << " L=" << L << " T=" << T << endl;
		cout << term << " term steps" << endl;
		cout << M << " measures" << endl;
		cout << "E=" << E << " E2=" << E2 << " Cv=" << Cv << endl;
		cout << "Added entry to " << filename << endl;
	}

	ofstream out;
	out.open(filename,ofstream::app);
	out << T << "\t" << E << "\t" << Cv << endl;
	out.close();

	if (calc_acorr)
	{
		cout << "Energy autocorrelation computed and saved to corr.dat" << endl;

		ofstream acorrf;
		acorrf.open("corr.dat");
		double norm = acorr[0]/counts[0]-E*E;
		for (int i = 0; i < tmax; ++i)
		{
			acorr[i] = (acorr[i]/counts[i]-E*E)/norm;
			acorrf << i << "\t" << acorr[i] << endl;
		}
		acorrf.close();
	}
}
	

double Sim::get_energy()
{
	double E = 0.;

	for (int i = 0; i < L; ++i)
	{	
		// loop only on black boxes
		for (int j = i%2; j < N; j += 2)
		{
			int Ph = (lattice.get(i,j)*lattice.get(i,j+1)) > 0;
			int Pv = (lattice.get(i,j)*lattice.get(i+1,j)) > 0;

			E += Hnd_matrix[Ph+Pv] + Jz*lattice.get(i,j)*lattice.get(i,j+1);
		}
	}

	return E/(double)L;
}

/* Given spins configuration u (on a plaquette), returns a random node g (compatible with u). */
int Sim::get_random_g(int u)
{
	double r = rand()/(double)RAND_MAX;
	
	int g = 0;
	double p = PUG[u][g];

	while (r > p)
	{
		g++;
		p += PUG[u][g];
	}

	return g;
}

void Sim::step()
{
	// reset nodes
	nodes = Grid(L,N);
	nodes.init(-1);
	for (int i = 0; i < L; ++i)
	{
		for (int j = i%2; j < N; j+=2)
		{
			// Get a random node g compatible with configuration u at black plaquette (i,j).
			int g = get_random_g(lattice.get_u(i,j));
			nodes.set(i,j,g);
		}
	}

	// reset labels
	labels = Loops(L,N);
	int last_id = 0;

	for (int i = 0; i < L; ++i)
	{
		for (int j = i%2; j < N; j+=2)
		{
			int g = nodes.get(i,j);
			labels.apply_node(i,j,g);
		}
	}

	// Clusterize loops 
	labels.reduce_loops();
	
	// Pick loops to flip
	for (int i = 0; i < labels.loops.size(); ++i)
	{
		if (labels.loops[i] != i)
			labels.loops[i] = 0;	
		else
		{
			double r = rand()/(double)RAND_MAX;
			labels.loops[i] = (r >= .5);
		}
	}

	// Apply flipping moves
	for (int i = 0; i < L; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			int l = labels.get(i,j);
			if (labels.loops[l])
				lattice.flip(i,j);
		}
	}

}

void Sim::compute_weights()
{
	conf[0] = -sinh(.5*Jx*tau)*exp(Jz*tau);
	conf[1] = cosh(.5*Jx*tau)*exp(Jz*tau);
	conf[2] = exp(-Jz*tau);

	Hnd_matrix[0] = -.5*Jx/tanh(.5*Jx*tau);
	Hnd_matrix[1] = -.5*Jx*tanh(.5*Jx*tau);
	Hnd_matrix[2] = 0.;
	
	node_weight[0] = .5*(conf[2]+conf[1]-conf[0]);
	node_weight[1] = .5*(-conf[2]+conf[1]+conf[0]);
	node_weight[2] = .5*(conf[2]-conf[1]+conf[0]);

	for (int i = 0; i < 3; ++i)
	{
		double den = 0.;
		for (int j = 0; j < 3; ++j)
			den += DUG[i][j]*node_weight[j];

		for (int j = 0; j < 3; ++j)
			PUG[i][j] = DUG[i][j] * (node_weight[j])/den;
	}

	if (debug_level != NONE)
	{	
		cout << "SPIN WEIGHTS (w0,w1,w2) = " << conf[0] << "," << conf[1] << "," << conf[2] << endl;
		cout << "NODE WEIGHTS (nw0,nw1,nw2) = " << node_weight[0] << "," << node_weight[1] << "," << node_weight[2] << endl;
		cout << "ENERGY MAT.EL. " << Hnd_matrix[0] << "," << Hnd_matrix[1] << "," << Hnd_matrix[2] << endl;
	}
}

void Sim::debug_output(string filename)
{
	ofstream nodes_f;
	nodes_f.open("debug_nodes_"+filename+".conf");
	for (int i = 0; i < L; ++i)
	{
		for (int j = i%2; j < N; j+=2)
			nodes_f << nodes.get(i,j) << " ";
		nodes_f << endl;
	}
	nodes_f.close();

	ofstream spin_f;
	spin_f.open("debug_spin_"+filename+".conf");
	for (int i = 0; i < L; ++i)
	{
		for (int j = 0; j < N; ++j)
			spin_f << (lattice.get(i,j)>0) << " ";
		spin_f << endl;
	}
	spin_f.close();
}

void Sim::set_Jx(double Jx) {this->Jx = Jx;}
void Sim::set_Jz(double Jz) {this->Jz = Jz;}
void Sim::set_termalization(int term){this->term = term;}
void Sim::set_nmeasures(int measures){this->measures = measures;}

void Sim::set_debug(DEBUG debug_level){this->debug_level = debug_level;}

void Sim::calc_autocorr(bool flag)
{
	if (flag == true)
	{
		// set flag to true
		calc_acorr = true;

		acorr = vector<double>(tmax,0.);
		counts = vector<int>(tmax,0);
	}
}

void Sim::set_acorr_time(int steps)
{
	// Reset to default
	if (steps <= 0 || steps >= measures) 
	{
		steps = measures/10;
		
		if (debug_level != NONE)
			cout << "WARNING: autocorr max timeshift overwritten to " << steps << endl;
	}

	this->tmax = steps;
	
	if (calc_acorr)
	{
		// reinit autocorrelation accumulators
		acorr = vector<double>(tmax,0.);
		counts = vector<int>(tmax,0);
	}
}

int main()
{

	clock_t begin = clock();
	int seed = time(NULL);
	cout << "Seed: " << seed << endl; 
	srand(seed);

	int L = 1000;
	int N = 10;
	double T = .01;

	Sim sim = Sim(L,N,T);
	// Physycal parameters
	sim.set_Jx(-10.);
	sim.set_Jz(0.);

	// Simulation parameters
	sim.set_termalization(100);
	sim.set_nmeasures(10000);

	// If set to true, computes normalized energy autocorrelation
	sim.calc_autocorr(false);
	sim.set_acorr_time(0); // Max timeshift to inspect

	// Debug level: VERBOSE (produces debug files), LIGHT (no debug files), NONE
	sim.set_debug(LIGHT);

	// Start simulation
	sim.run();

	// Produce output
	sim.save_output();

	double cputime = double(clock()-begin)/CLOCKS_PER_SEC*1000;
	cout << "CPU time: " << cputime << " ms" << endl;
	cout << "Bye" << endl;
	return 0;
}