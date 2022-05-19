#include "nc_data.h"

NcDataFile::NcDataFile(string file_name){

	dataFile.open(file_name,NcFile::replace);

}

void NcDataFile::Init(string file_name){

	dataFile.open(file_name,NcFile::replace);

}

void NcDataFile::InitNcData(int Ny, int Nx){
	
	xDim = dataFile.addDim(xDimName, Nx);
    yDim = dataFile.addDim(yDimName, Ny);
    recDim = dataFile.addDim(recName);
    
    dimVector.push_back(recDim);
    dimVector.push_back(yDim);
    dimVector.push_back(xDim);
    
    xtVar = dataFile.addVar(xCoordinateName, ncDouble, xDim);
    ytVar = dataFile.addVar(yCoordinateName, ncDouble, yDim);
    uVar = dataFile.addVar("u", ncDouble, dimVector);
    
    dimVector.push_back(recDim);
    dimVector.push_back(yDim);
    dimVector.push_back(xDim);
    
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(1);
    countp.push_back(Ny);
    countp.push_back(Nx);
    
    
}

void NcDataFile::WriteCoordData(matrix & xt, matrix & yt){

	xtVar.putVar(&xt(0));
	ytVar.putVar(&yt(0));


}

void NcDataFile::WriteNcData(matrix & A, int k){

	startp[0] = k;
	
	uVar.putVar(startp,countp,&A(0,0));


}



