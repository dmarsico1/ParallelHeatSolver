#ifndef _GRID_H_
#define _GRID_H_

#include <map>
#include "matrix.h"

using namespace std;

class Grid
{

	public:
		
		int * Nx; //Grid cells in x
		int * Ny; //Grid cells in y
		
		double * Lx; //Length of domain in x
		double * Ly; //Length of domain in y
		
		double hx; //Grid size in x
		double hy; //Grid size in y
		
		double * dt; //Step size
		double * Nt; //Stopping time
		double * write_freq; //frequency with which data is written
		
		matrix hx_vec;
		matrix hy_vec;
		
		std::map<std::string,int> GridValuesInt; //Dictionary of integer grid parameters
		std::map<std::string,double> GridValuesDouble; //Dictionary of double grid parameters
		
		int dim_freq;
		matrix write_freq_vec;
		
		Grid();
		
		//Initialize default values for grid
		void InitGridValues();
		
		//Read input and change default grid values if needed
		void ReadInput();
			
		//Define coordinate vectors
		void InitCoordinates();
		
		void InitFrequencyVector();
};

#endif

