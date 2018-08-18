// Metropolis algo with calculation of energy
#include <fstream>
#include <iostream>
#include <random>
#include <math.h>
using namespace std;
#define N 100
#define K 1
#define ITER 100000000
#define TempPoints 30
//#define kT 1.1
//2.27

//random seed generators
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis (0.0, 1.0);
double kT [TempPoints] = {0.1, 0.5, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.22, 2.24, 2.26, 2.27, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 4.9};
short lattice[N*N + 1];
short neighbors[4];
double mag[3];	//store moments of magnetization- <mag>, <mag^2>, <mag^4>
double energy[3]; // store moments of energy- <energy>, <energy^2>, <energy^4>

void create (short* arr)
{
	for (int i=0; i<N*N; i++)
		if (dis(gen)<0.5)
			arr[i]=-1;
		else
			arr[i]=+1;
	arr[N*N]=0;
}

short* findAdjacent (int x, short* arr)
{
	int row = x/N, col = x % N;	
	//neighbors of lattice[x]
	int l= N*row+ col-1;
	int r= N*row+ col+1;
	int u= N*(row-1)+ col;
	int d= N*(row+1)+ col;

	// set boundaries
	if (row == 0)
		u = N*N;
	if (col == 0)
		l = N*N;
	if (row == N-1)
		d = N*N;
	if (col == N-1)
		r = N*N;
	neighbors[0] = arr[l];
	neighbors[1] = arr[r];
	neighbors[2] = arr[u];
	neighbors[3] = arr[d];
	return neighbors;
} 

void update (int x, short* lattice, double kT)
{	
	//change in energy after update
	int dE; 
	// update probability
	double p;
	// lattice coordinates
	findAdjacent (x, lattice);
	//Calculate change in energy
	dE = 2*K*lattice[x]*(neighbors[0] + neighbors[1] + neighbors[2] + neighbors[3]);
	//Flip if change is negative
	if (dE<0)
		lattice[x] = -lattice[x];
	else
	{
		p = exp(-double(dE)/kT);
		//Flip with probability p if change is positive
		if (dis(gen) <= p)
			lattice[x] = -lattice[x];	
	}
}

double* magnetization (short* arr)
{
	double sum=0;
	for (int i=0; i<N*N; i++)
		sum += arr[i];		// net magnetization of system
	mag[0] = sum/(N*N);
	mag[1] = pow (sum, 2)/(N*N);
	mag[2] = pow (sum, 4)/(N*N); 
	return mag;
}

double* isingEnergy (short* arr)
{
	double sum = 0;
	for (int i=0; i<N*N; i++)
	{	
		findAdjacent (i, lattice);
		sum += -K*lattice[i]*(neighbors[0] + neighbors[1] + neighbors[2] + neighbors[3]);		// energy of system
	}
	energy[0] = sum/(N*N);
	energy[1] = pow (sum, 2)/(N*N);
	energy[2] = pow (sum, 4)/(N*N);
	return energy; 
}

int main ()
{
	ofstream opfile;
  	opfile.open ("IsingDat2.ods");
  	opfile << "N=100, iter=100000000" << endl;

	int iter=0, T=0;
	for (T=0; T<TempPoints; T++) 
	{
		create(lattice);
		opfile << "kT =" << kT[T] << endl;
		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
				opfile << lattice[N*i + j] << "\t";
			opfile << endl;
		}	// print lattice
		
		opfile << "iteration \t" << "<energy> \t" << "<energy^2> \t" << "<energy^4> \t" << "<mag> \t" << "<mag^2> \t" << "<mag^4>" << endl;
		for (iter=0; iter < ITER; iter++)
		{
			update(iter%(N*N), lattice, kT[T]);
			if (iter%(N*N)==0)
			{
				magnetization(lattice);
				isingEnergy(lattice);
				opfile << iter <<"\t"<< energy[0] <<"\t"<< energy[1] <<"\t"<< energy[2] <<"\t"<< mag[0] <<"\t"<< mag[1] <<"\t"<< mag[2] << endl;			
			}
			// print iteration number, energy and magnetization
		}
		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
				opfile << lattice[N*i + j] << "\t";
			opfile << endl;
		}	// print lattice
	}	
	opfile.close();
	return 0;
}