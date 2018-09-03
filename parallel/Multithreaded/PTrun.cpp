#include "omp.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <random>
#include <math.h>
#include <complex>
#include <string>
#include <array>
using namespace std;

#define K 1
#define TempPoints 10

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis (0.0, 1.0);
typedef complex<double> dcomp;
double kT[TempPoints] = {0.5, 1.7, 2.15, 2.2, 2.25, 2.3, 2.35, 2.6, 3.6, 4.6};
//double kT [TempPoints] = {0.1, 0.5, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.22, 2.24, 2.26, 2.27, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 4.9};

void create (short* arr, int N)
{
	for (int i=0; i<N*N; i++)
		if (dis(gen)<0.5)
			arr[i]=-1;
		else
			arr[i]=+1;
}

void swapLattices (short* a1, short* a2)
{
	short* temp;
	temp = a1;
	a1 = a2;
	a2 = temp;
	return;
}

void findAdjacent (int x, short* arr, short* neighbors, int N)
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
	return;
} 

void update (int x, short* lattice, double kT, short* neighbors, int N)
{	
	//change in energy after update
	int dE; 
	// update probability
	double p;
	// lattice coordinates
	findAdjacent (x, lattice, neighbors, N);
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

double* magnetization (short* arr, double* mag, int N)
{
	double sum=0;
	for (int i=0; i<N*N; i++)
		sum += arr[i];		
	mag[0] = abs(sum/(N*N));	// absolute magnetization of system
	mag[1] = pow (mag[0], 2); 	//square of magnetization of system
	mag[2] = pow (mag[0], 4); 	//fourth power of magnetization of system
	return mag;
}

double* isingEnergy (short* arr, short* neighbors, double* energy, int N)
{
	double sum = 0;
	for (int i=0; i<N*N; i++)
	{	
		findAdjacent (i, arr, neighbors, N);
		sum += -K*arr[i]*(neighbors[0] + neighbors[1] + neighbors[2] + neighbors[3]);		// energy of system
	}
	energy[0] = sum/(N*N);
	energy[1] = pow (energy[0], 2);
	energy[2] = pow (energy[0], 4);
	return energy; 
}

double Chi (short* arr, int N)
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

int main(int argc, char** argv)
{
	omp_set_num_threads(5);
	int runs = atoi(argv[1]);
	#pragma omp parallel
	{
		int ID = omp_get_thread_num();
		int N = 20*(ID+1) ;	// Lattice size
		string fname[5];
		ofstream opfile;
		int flag;
		double delE, delBeta, prob;
		short **arr = new short*[TempPoints];
		short **adjacents = new short*[TempPoints];
		double** magLat = new double*[TempPoints];
		double** energyLat = new double*[TempPoints];
		double** WinEn = new double* [TempPoints];
		double** WinMag = new double* [TempPoints]; 
		double chi[TempPoints], WinChi[TempPoints];
		double c1Win[TempPoints], c2Win[TempPoints], energyVar[TempPoints];
		double windowsize;
		
		for (int iter=0; iter<TempPoints; iter++)
		{	
			arr[iter] = new short [N*N];
			adjacents[iter] = new short[4];
			magLat[iter] = new double[3];
			energyLat[iter] = new double[3]; 
			WinEn[iter] = new double[3];
			WinMag[iter] = new double[3];
		}

		for (int r=0; r<runs; r++)
		{
			fname[ID] = "PT_";
			fname[ID] += to_string(N);
			fname [ID]+= "_";
			fname [ID]+= to_string(r);
			fname [ID]+= ".ods";
			opfile.open(fname[ID]);
			opfile << "N = " << N << endl; 

			flag=1;

			for (int iter=0; iter<TempPoints; iter++)
				create(arr[iter], N);	//initialize and create lattices
			
			for (int window = 0; ; window++)
			{
				windowsize = pow(2, window); 
				for (int T=0; T<TempPoints; T++)
				{
					for (int i=0; i<3; i++)
					{
						WinEn[T][i] = 0;
						WinMag[T][i] = 0;
					}
					WinChi[T] = 0;
				}
				opfile << "Window \t Temp \t <energy> \t EnergyVar \t <energy^2> \t <energy^4> \t <mag> \t <mag^2> \t <mag^4> \t <chi(k)> \t spec heat" << endl;
				for (int iter=0; iter<windowsize; iter++)
				{	
					for (int T=0; T<TempPoints; T++)
					{
						for (int i=0; i<N*N; i++)
							update(i, arr[T], kT[T], adjacents[T], N);	//Lattice sweep

						isingEnergy(arr[T], adjacents[T], energyLat[T], N);	//Energy 
						magnetization(arr[T], magLat[T], N);	//Mag calculation
						chi[T] = Chi(arr[T], N);

						for (int i=0; i<3; i++)
						{
							WinEn[T][i] += energyLat[T][i]/windowsize;
							WinMag[T][i] += magLat[T][i]/windowsize;		// storing window averages for energy and mag at a Temp Point
						}
						WinChi[T] += chi[T]/windowsize;			
					}		
					for (int k=0; k<TempPoints-1; k++)
					{
						delE = energyLat[k+1][0]-energyLat[k][0];
						delBeta = (double)1/kT[k+1] - (double)1/kT[k];
						prob = exp(delBeta*delE);
						if (prob > 1)
							swapLattices(arr[k], arr[k+1]);
						else
							if (dis(gen)<prob)
								swapLattices(arr[k], arr[k+1]);
					}					// 1 Time step	
				}			// 1 window = (2^window) time steps
				
				for (int T=0; T<TempPoints; T++)
				{
					energyVar[T] = WinEn[T][1] - pow(WinEn[T][0], 2); 
					opfile << window << "\t" << kT[T] << "\t" << WinEn[T][0] << "\t" << energyVar[T] << "\t" << WinEn[T][1] << "\t" << WinEn[T][2] << "\t" << WinMag[T][0] << "\t" << WinMag[T][1] << "\t" << WinMag[T][2] << "\t" << WinChi[T] << "\t" << energyVar[T]/(pow(kT[T],2)) << endl;
				}			// print window averages and variance in energy for all Temp
				
				if (window==0)
				{
					for (int T=0; T<TempPoints; T++)
						c1Win[T] = WinEn[T][0];
					continue;
				}
				
				flag = 0;
				for (int T=0; T<TempPoints; T++)
				{
					double check1, c1Winsq;
					c1Winsq = pow(c1Win[T], 2);
					check1 = abs(c1Winsq - pow(WinEn[T][0], 2));

					if (check1<=energyVar[T])
						flag = 0;
					else
					{
						flag = 1;
						break;
					}
				}
				if (flag==0)
					break;	
				for (int T=0; T<TempPoints; T++)
					c1Win[T] = WinEn[T][0];
			}

			opfile.close();

		}

		for (int iter=0; iter<TempPoints; iter++)
		{	
			delete(arr[iter]);
			delete(adjacents[iter]);
			delete(magLat[iter]);
			delete(energyLat[iter]);
			delete(WinEn[iter]);
			delete(WinMag[iter]);
		}

		delete(arr);
		delete(adjacents);
		delete(magLat);
		delete(energyLat);
		delete(WinEn);
		delete(WinMag);			
	}
	return 0;
}