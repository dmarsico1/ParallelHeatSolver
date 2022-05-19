#include <string>
#include <iostream>
#include <fstream>
#include "grid.h"

using namespace std;

bool IsInt(string input){
	
	for (int i = 0; i < input.size(); i++){
		
		if(input[i]=='.'){
			
			return false;
			
		}
		
	}
	
	return true;
	
}

Grid::Grid(){}

void Grid::InitGridValues(){

	GridValuesInt["Nx"] = 50;
	GridValuesInt["Ny"] = 50;
	
	GridValuesDouble["Lx"] = 1.0;
	GridValuesDouble["Ly"] = 1.0;
	GridValuesDouble["Nt"] = 20.0;
	GridValuesDouble["dt"] = 0.2;
	GridValuesDouble["write_freq"] = 1.0;
	
	Nx = &GridValuesInt["Nx"];
	Ny = &GridValuesInt["Ny"];
	
	Lx = &GridValuesDouble["Lx"];
	Ly = &GridValuesDouble["Ly"];
	Nt = &GridValuesDouble["Nt"];
	dt = &GridValuesDouble["dt"];
	write_freq = &GridValuesDouble["write_freq"];

}

void Grid::ReadInput(){
	
	//Read in input file
	std::ifstream input_stream("input.txt");
	
	if(input_stream.is_open()){
		
		std::string line;
		
		while (getline(input_stream, line)) {
			
			int VarSize = 0;
			string VariableName;
			string VariableValue;
			
			//Get length of variable name
			for (int i = 0; i < line.size(); i++){
				
				if(line[i] == ':'){
					break;
				}
				else{
					VarSize++;
				}
			}
				
			VariableName = line.substr(0,VarSize);
			VariableValue = line.substr(VarSize+1,line.size()-1);
			
			//Determine if input value is an integer, and write to integer input dictionary if it is
			if(IsInt(VariableValue)){
				
				GridValuesInt[VariableName]=stoi(VariableValue);
				
			}
			
			else{
				
				GridValuesDouble[VariableName]=stod(VariableValue);
				
			}
			
		}
		
	}

}

void Grid::InitCoordinates(){
	
	this->hx_vec.allocate(*this->Nx);
	this->hy_vec.allocate(*this->Ny);
	
	this->hx = 2.0*(*this->Lx)/static_cast<double>(*this->Nx);
	this->hy = 2.0*(*this->Ly)/static_cast<double>(*this->Ny);
	
	for (int i = 0; i < *this->Nx; i++){
		
		this->hx_vec(i) = -*this->Lx + (i+0.5)*this->hx;
		
	}
	
	for (int i = 0; i < *this->Ny; i++){
		
		this->hy_vec(i) = -*this->Ly + (i+0.5)*this->hy;
		
	}
	
	
}

void Grid::InitFrequencyVector(){
	
	//dimension of the frequency vector (minus 1)
	dim_freq=static_cast<int>(floor(*Nt/(*write_freq)));
	
	//allocate frequency vector
	write_freq_vec.allocate(dim_freq+1);
	
	for (int i = 0; i <= dim_freq; i++){
		write_freq_vec(i) = (i*(*write_freq));
	}

}





