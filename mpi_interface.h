#ifndef _MPI_INTERFACE_H_
#define _MPI_INTERFACE_H_

#include <math.h>
#include <mpi.h>
#include "matrix.h"
#include <vector>

using namespace std;

class GridTopology
{
	
	public:
		int upper; //upper neighbor id
		int lower; //lower neighbor id
		int left; //left neighbor id
		int right; //right neighbor id
		int myid; //process id
		int nproc; //number of processes
		int reorder;
		
		int local_disp; //displacement of local grid in global grid
		vector<int> disps;  //vector of local displacements; only relevant at root process
		vector<int> counts; //vector of ones used in the collective gather; only relevant at root
	
		int x_dim; //number of cells in x dimension of current process
		int y_dim; //number of cells in y dimension of current process
	
		MPI_Comm comm2d;
		vector<int> coords; // process coordinates within the grid topology
		vector<int> dims; //number of processes in each dimension
		vector<int> periods; //periodicity in specified dimension
	
		matrix hx_vec_local; //x coordinates of each process
		matrix hy_vec_local; //y coordinates of each process
		
		matrix p;
	
		double gamma = 0.0; //sum used for tolerance calculation in conjugate gradient
		double gamma_tot = 0.0;
		double rnew = 0.0;
		double rnew_tot = 0.0;
		double rold = 0.0;
		double rold_tot = 0.0;
		double alpha_tot = 0.0;
	
		GridTopology(int order);
		void CreateTopology(int Nx, int Ny, double hy, double hx, double Ly, double Lx);
		void GetDisplacements(int Ny, int Nx);
		void GetFullArray(matrix & FullArray, matrix & SubArray, int Ny, int Nx);

};

void Neighbors(matrix & A, GridTopology & Grid);

void ConjGrad(matrix & x, matrix & f, int M, int N, double hx, double hy, double dt, GridTopology & Grid);

#endif
