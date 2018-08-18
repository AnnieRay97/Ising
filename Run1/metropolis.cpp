#include <fstream>
#include <iostream>
#include <random>
#include <math.h>
using namespace std;
#define N 40
#define K 1
//#define kT 1.1
//2.27

//random seed generators
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis (0.0, 1.0);
double kT [13] = {0.5, 1.0, 1.2, 1.5, 1.7, 1.9, 2.0, 2.1, 2.2, 2.4, 2.5, 2.8, 3.1};
short lattice[N*N + 1];

void create (short* arr)
{
	for (int i=0; i<N*N; i++)
		if (dis(gen)<0.5)
			arr[i]=-1;
		else
			arr[i]=+1;
	arr[N*N]=0;
}

void update (int x, short* arr, double kT)
{	
	//change in energy after update
	int dE; 
	// update probability
	double p, a, b;
	// lattice coordinates
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
	//Calculate change in energy
	dE = 2*K*arr[x]*(arr[l] + arr[r] + arr[u] + arr[d]);
	//Flip if change is negative
	if (dE<0)
		arr[x] = -arr[x];
	else
	{
		p = exp(-double(dE)/kT);
		//Flip with probability p if change is positive
		if (dis(gen) <= p)
			arr[x] = -arr[x];	
	}
}

double magnetization (short* arr)
{
	double mag = 0;
	for (int i=0; i<N*N; i++)
		mag += arr[i];
	return mag/(N*N);
}

int main ()
{
	ofstream opfile;
  	opfile.open ("IsingDat.ods");
  	opfile << "N=40, iter=10000000" << endl;

	int iter=0, T=0;
	for (T=0; T<13; T++) 
	{
		create(lattice);
		opfile << "kT =" << kT[T] << endl;
		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
				opfile << lattice[N*i + j] << "\t";
			opfile << endl;
		}	// print lattice
		
		for (iter=0; iter < 10000000; iter++)
		{
			/*for (x=0; x<N*N; x++)
			{
			*/	update(iter%(N*N), lattice, kT[T]);
				if (iter%(N*N)==0)
					opfile << magnetization(lattice) << "\t" << iter << endl;			
			//}
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