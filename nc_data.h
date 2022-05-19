#ifndef _NC_DATA_H_
#define _NC_DATA_H_

#include "grid.h"
#include "matrix.h"
#include <netcdf>
#include <vector>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

class NcDataFile{

	public:
		
		NcFile dataFile;
		
		vector<NcDim> dimVector;
		
		string xDimName = "Nx";
		string yDimName = "Ny";
		string recName = "time";
		string xCoordinateName = "xt";
		string yCoordinateName = "yt";
		
		NcDim xDim;
    	NcDim yDim;
    	NcDim recDim;
    
    	NcVar xtVar;
    	NcVar ytVar;
    	NcVar uVar;
    	
		//Dimensions of variable netcdf file
		vector<size_t> startp;
		vector<size_t> countp;
		
		NcDataFile(){};
		NcDataFile(string file_name);
		
		void Init(string file_name);
		
		void InitNcData(int Ny, int Nx);
		
		void WriteCoordData(matrix & xt, matrix & yt);
		
		void WriteNcData(matrix & A, int k);

};


#endif

