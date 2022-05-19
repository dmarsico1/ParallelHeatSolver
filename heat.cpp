#include "grid.h"
#include "mpi_interface.h"
#include "nc_data.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <mpi.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


int main(int argc, char** argv){
	
	//Initialize MPI
	MPI_Init(&argc,&argv);
	
	//Declare global grid
	Grid GlobalGrid;
	
	//Initialize global grid values
	GlobalGrid.InitGridValues();
	
	//Read non-default values if necessary
	GlobalGrid.ReadInput();
	
	//Initalize global grid coordinates
	GlobalGrid.InitCoordinates();
	
	//Vector of times at which to write the data
	GlobalGrid.InitFrequencyVector();
	
	//Declare MPI subgrids
	GridTopology SubGrid(1);
	
	//Define MPI cartesian topology
	SubGrid.CreateTopology(*GlobalGrid.Nx,*GlobalGrid.Ny, GlobalGrid.hy, GlobalGrid.hx, 
						   *GlobalGrid.Ly, *GlobalGrid.Lx);
	
	//Get locations in memory for each process
	SubGrid.GetDisplacements(*GlobalGrid.Ny, *GlobalGrid.Nx);
	
	//Declare solution and forcing arrays
	matrix u;
	matrix f;
	
	//Allocate arrays
	u.allocate(SubGrid.y_dim,SubGrid.x_dim);
	f.allocate(SubGrid.y_dim,SubGrid.x_dim);
	
	u = 0.0;
	
	//Set initial conditions and forcing
	InitialConditions(u, SubGrid.y_dim, SubGrid.x_dim, *GlobalGrid.Lx, SubGrid.hx_vec_local, SubGrid.hy_vec_local);
	
	forcing(f, SubGrid.y_dim, SubGrid.x_dim, *GlobalGrid.Lx, SubGrid.hx_vec_local, SubGrid.hy_vec_local);
	//u=f;
	//Declare netcdf data
	string file_name = "heat.nc";
	NcDataFile data;
	
	//Define netcdf data on root process
	if(SubGrid.myid==0){
		
		data.Init(file_name);
		data.InitNcData(*GlobalGrid.Ny, *GlobalGrid.Nx);
		data.WriteCoordData(GlobalGrid.hy_vec,GlobalGrid.hy_vec);

	}
	
	//Define global array used for the write operation
	matrix full_array;
	if(SubGrid.myid==0){
	
		full_array.allocate(*GlobalGrid.Ny,*GlobalGrid.Nx);
	
	}
	
	int k = 0; //record dimension index
	
	double time = 0.0;
	
	double write_tol = 0.00001;
	
	//Main time-stepper
    while(time < *GlobalGrid.Nt){
		
		matrix rhs = u + f*(*GlobalGrid.dt);
		
		ConjGrad(u,rhs,SubGrid.y_dim,SubGrid.x_dim,GlobalGrid.hx,GlobalGrid.hy,*GlobalGrid.dt, SubGrid);
		
		if((time < GlobalGrid.write_freq_vec(k) + write_tol) 
			&& (time > GlobalGrid.write_freq_vec(k) - write_tol)){
			
			//Gather all data onto root process
			SubGrid.GetFullArray(full_array, u, *GlobalGrid.Ny, *GlobalGrid.Nx);
			
			//Write data
			if(SubGrid.myid == 0){
			
				data.WriteNcData(full_array,k);
				
			}
			
			k += 1;
			
		}
		
		time += *GlobalGrid.dt;
						
    }
	
	MPI_Finalize();
	
	return 0;
	
}

 
