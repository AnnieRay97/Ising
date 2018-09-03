// Metropolis algo with calculation of energy
#include <fstream>
#include <iostream>
#include <random>
#include <math.h>
#include <complex>

using namespace std;

typedef complex<double> dcomp;

#define N 10
#define K 1
#define TempPoints 30
#define MAXWIN 25
//#define kT 1.1
//2.27

//random seed generators
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis (0.0, 1.0);
double kT [TempPoints] = {0.1, 0.5, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.22, 2.24, 2.26, 2.27, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 4.9};
short lattice[N*N];
short neighbors[4];
double mag[3];	//store moments of magnetization- <mag>, <mag^2>, <mag^4>
double energy[3]; // store moments of energy- <energy>, <energy^2>, <energy^4>
double chiCorel;
double WinMag[3];	//store window average of the moments of magnetization
double WinEn[3];	//store window average of the moments of energy
double energyVar;
double WinChiCorel;
void create (short* arr)
{
	for (int i=0; i<N*N; i++)
		if (dis(gen)<0.5)
			arr[i]=-1;
		else
			arr[i]=+1;
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
		u = N*(N-1)+col;
	if (col == 0)
		l = N*(row)+N-1;
	if (row == N-1)
		d = col;
	if (col == N-1)
		r = N*(row);
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
		sum += arr[i];		
	mag[0] = sum/(N*N);	// magnetization of system
	mag[1] = pow (mag[0], 2); //square of magnetization of system
	mag[2] = pow (mag[0], 4); // fourth power of magnetization of system
	return mag;
}

double Chi (short* arr)
{
	dcomp i, chi = 0, b;
	double exponent;
	i = -1;
	i = sqrt(i);
	double pi = 2*asin(1);
	for (int iter=0; iter<N*N; iter++)
	{
		exponent = 2*(pi/N)*((iter)%N);
		//exponent = 0;	// k=0 gives X(0) = N*N*<m^2>
		b = exp(exponent*i);
		b *= arr[iter];
		b /= N*N;
		chi +=  b;
	}
	return N*N*pow(abs(chi),2);
	return 0;
}

double* isingEnergy (short* arr)
{
	double sum = 0;
	for (int i=0; i<N*N; i++)
	{	
		findAdjacent (i, lattice);
		sum += -K*arr[i]*(neighbors[0] + neighbors[1] + neighbors[2] + neighbors[3]);		// energy of system
	}
	energy[0] = sum/(N*N);
	energy[1] = pow (energy[0], 2);
	energy[2] = pow (energy[0], 4);
	return energy; 
}

int main ()
{
	ofstream opfile;
  	opfile.open ("IsingDat3_lat10_run1.csv");
  	opfile << "N=10" << endl;

	int iter=0, T=0, window=0;
	for (T=0; T<TempPoints; T++) 
	{
	 	create(lattice);
		opfile << "kT =" << kT[T] << endl;
		opfile << "window , <energy> , EnergyVar , <mag> , <energy^2> , <mag^2> , <energy^4> , <mag^4> , <chi(k)>" << endl;
		for (window=0;  ; window++)
		{
			double c1WinEn, c2WinEn, c1Winsq, c2Winsq, check1, check2;
			int windowsize = pow(2,window);
			for (int k=0;k<3;k++)
			{
				WinEn[k]=0; 
				WinMag[k]=0;
			}	// Initialize window averages to 0
			WinChiCorel = 0;
			for (int i=0; i<windowsize; i++)
			{
				for (iter=0; iter < N*N; iter++)
					update(iter, lattice, kT[T]);	//1 time step
				
				magnetization(lattice);
				isingEnergy(lattice);
				chiCorel = Chi(lattice);
				for (int k=0; k<3; k++)
				{
					WinEn[k] += energy[k]/windowsize; 
					WinMag[k] += mag[k]/windowsize;
				}	//storing Window average of magnetization and energy
				WinChiCorel += chiCorel/windowsize; //storing Window average of chi
			}
			energyVar = WinEn[1]-pow(WinEn[0],2);	//Variance(energy) of this window
			opfile << window << "," << WinEn[0] << "," << energyVar << "," << WinMag[0] << "," << WinEn[1] << "," << WinMag[1] << "," << WinEn[2] << "," << WinMag[2] << "," << WinChiCorel << endl;	// Window avg of chi to be used for correlation length

			if (window == 0)
			{
				c1WinEn = WinEn[0]; 
				continue;
			}
			if (window == 1)
			{
				c2WinEn = WinEn[0]; 
				continue;
			}		// Store values of Avg Energy for last two windows 
			// if (c1WinEn==c2WinEn & c2WinEn==WinEn[0] )
			// 	break;	// If the values of Mean Energy match for last 3 windows stop!
			
			c1Winsq = pow(c1WinEn,2);
			c2Winsq = pow(c2WinEn,2);
			check1 = abs(pow(WinEn[0],2)-c2Winsq);
			check2 = abs(c2Winsq-c1Winsq);
			cout << "c1Win ," << c1WinEn << "c1WinSq ," << c1Winsq << "c2Win ," << c2WinEn << "c1WinSq ," << c1Winsq<< endl; 
			cout << "Current-Last ," << check1 << endl;
			cout << "Last-SecondLast ," << check2 << endl;
			cout << energyVar << endl;
			if ((check1 <= energyVar) and (check2<=energyVar))
				break;	// If the values of Mean Energy match for last 3 windows stop! 
			c1WinEn = c2WinEn;
			c2WinEn = WinEn[0];			
		}
	}	
	opfile.close();
			
	return 0;
}