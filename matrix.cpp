#include "matrix.h"
#include <iostream>


matrix::matrix(){nrows = 0; ncols = 0;mat=NULL;}

matrix::matrix(int dim){
	
	nrows = dim;
	ncols = 0;
	mat = new double[dim];
	
}

matrix::matrix(int dim, double a){
	
	nrows = dim;
	ncols = 0;
	mat = new double[dim];
	
	for (int i = 0; i < dim; i++){
		
		mat[i] = a;
		
	}
	
}

matrix::matrix(int rows, int cols){
	
	nrows = rows;
	ncols = cols;
	mat = new double[rows*cols];
	
}

matrix::~matrix() {

	if(mat == NULL){
		delete mat;
	}
	else{
 		delete [] mat;
 	}
 	
}

int matrix::GetRows(){
	
	return this->nrows;
	
}

int matrix::GetCols(){
	
	return this->ncols;
	
}

void matrix::allocate(int rows){
	
	this->nrows = rows;
	this->ncols = 0;
	mat = new double[rows];
	
}

void matrix::allocate(int rows, int cols){
	
	this->nrows = rows;
	this->ncols = cols;
	mat = new double[rows*cols];
	
}


matrix::matrix(int rows, int cols, double a){
	
	nrows = rows;
	ncols = cols;
	mat = new double[rows*cols];
	
	for (int i = 0; i < nrows; i++){
		
		for (int j = 0; j < ncols; j++)	{
			
			mat[i*ncols+j] = a;
		
		}
		
	}
	
}

double & matrix::operator() (int irow){
	
	return this->mat[irow];
	
}

const double & matrix::operator() (int irow) const {
	
	return this->mat[irow];
	
}

double & matrix::operator() (int irow, int icol){
	
	return this->mat[irow*this->ncols + icol];
	
}

const double & matrix::operator() (int irow, int icol) const {
	
	return this->mat[irow*this->ncols + icol];
	
}

//matrix & matrix::operator= (const matrix &A){
	
	//for (int i = 0; i < this->nrows; i++){
		//for (int j = 0; j < this->ncols; j++){
				
			//this->mat[i*this->ncols + j] = A.mat[i*this->ncols + j];
			
		//}
	//}
	
	//return *this;
	
//}

matrix & matrix::operator= (const matrix & A){
	
	if(this->nrows != A.nrows){
		
		this->nrows = A.nrows;
		
	}
	
	if(this->ncols != A.ncols){
		
		this->ncols = A.ncols;
		
	}
	
	
	for (int i = 0; i < this->nrows; i++){
		for (int j = 0; j < this->ncols; j++){
				
			this->mat[j*this->nrows + i] = A.mat[j*this->nrows + i];
			
		}
	}
	
	return *this;
	
}

matrix & matrix::operator= (double a){
	
	for (int i = 0; i < this->nrows; i++){
		for (int j = 0; j < this->ncols; j++){
				
			this->mat[i*this->ncols + j] =a;
			
		}
	}
	
	return *this;
	
}

matrix & matrix::operator+= (const matrix &A){
	
	for (int i = 0; i < this->nrows; i++){
		for (int j = 0; j < this->ncols; j++){
			
			this->mat[i*this->ncols+j] = this->mat[i*this->ncols + j] + A(i,j);
			
		}
	}
	
	return *this;
	
	
}

matrix & matrix::operator-= (const matrix &A){
	
	for (int i = 0; i < this->nrows; i++){
		for (int j = 0; j < this->ncols; j++){
			
			this->mat[i*this->ncols+j] = this->mat[i*this->ncols + j] - A.mat[i*this->ncols + j] ;
			
		}
	}
	
	return *this;
	
	
}

matrix matrix::operator+ (const matrix &A){
	
	matrix result(this->nrows,this->ncols);
	
	for (int i = 0; i < this->nrows; i++){
		for (int j = 0; j < this->ncols; j++){
			
			result(i,j) = this->mat[i*this->ncols + j] + A.mat[i*this->ncols + j] ;
			
		}
	}
	
	return result;
	
}

matrix matrix::operator- (const matrix &A){
	
	matrix result(this->nrows,this->ncols);
	
	for (int i = 0; i < this->nrows; i++){
		for (int j = 0; j < this->ncols; j++){
			
			result(i,j) = this->mat[i*this->ncols + j] - A.mat[i*this->ncols + j];
			
		}
	}
	
	return result;
	
}

matrix matrix::operator* (double a){
		
	matrix result(this->nrows,this->ncols);
		
	for (int i = 0; i < nrows; i++){
		
		for (int j = 0; j < ncols; j++){
			
			//this->mat[i*this->nrows + j] = a*(this->mat[i*this->nrows + j]);
			result(i,j) = a*(this->mat[i*this->ncols + j]);
		
		}
	}
	
	//return *this;
	return result;
	
}

matrix matrix::operator* (const matrix &A){
	
	matrix mat_prod(this->nrows,this->ncols);
	
	for (int i = 0; i < this->nrows; i++){
		for (int j = 0; j < this->ncols; j++){
			
			mat_prod(i,j) = (this->mat[i*this->ncols + j])*A.mat[i*this->ncols + j];
			
		}
	}
	
	return mat_prod;
	
}

//void matrix::Add(matrix & A, matrix & B, int startI, int startJ, int countI, int countJ){
	
	//for (int i = 0; i < countI; i++){
		//for (int j = 0; j < countJ; j++){
			
			//this->mat[i*(this->ncols)+j+startI*this->ncols+startJ] = A(i+startI,j+startJ) + B(i+startI,j+startJ);
			
		//}
	//}
//}

void matrix::Add(matrix & A, matrix & B, vector<int> IStart, vector<int> JStart, int countI, int countJ){
	
	// IStartEnd is a 3x2 matrix such that row k gives the start index
	// for each matrix.
	
	int startI = IStart[0];
	int startIA = IStart[1];
	int startIB = IStart[2];
	
	int startJ = JStart[0];
	int startJA = JStart[1];
	int startJB = JStart[2];

	for (int i = 0; i < countI; i++){
		for (int j = 0; j < countJ; j++){
			
			this->mat[i*(this->ncols)+j+startI*this->ncols+startJ] = A(i+startIA,j+startJA) + B(i+startIB,j+startJB);
			
		}
	}
	
	
}

void matrix::Subtract(matrix & A, matrix & B, vector<int> IStart, vector<int> JStart, int countI, int countJ){
	
	// IStartEnd is a 3x2 matrix such that row k gives the start index
	// for each matrix.
	
	int startI = IStart[0];
	int startIA = IStart[1];
	int startIB = IStart[2];
	
	int startJ = JStart[0];
	int startJA = JStart[1];
	int startJB = JStart[2];

	for (int i = 0; i < countI; i++){
		for (int j = 0; j < countJ; j++){
			
			this->mat[i*(this->ncols)+j+startI*this->ncols+startJ] = A(i+startIA,j+startJA) - B(i+startIB,j+startJB);
			
		}
	}
	
	
}


matrix matrix::GetSubArray(int startI, int startJ, int countI, int countJ){
	
	matrix result(countI,countJ);
	
	for (int i = 0; i < countI; i++){
		for (int j = 0; j < countJ; j++){
			
			result(i,j) = this->mat[i*(this->ncols)+j+startI*this->ncols+startJ];
			
		}
	}
	
	return result;
	
}

void matrix::SetRow(matrix &A, int startI, int startJ, int countJ){
	
	for (int j = 0; j < countJ; j++){
			
		this->mat[startI*(this->ncols)+startJ+j] = A(0,j);//A(startI,startJ+j);
		
	}
	
}

void matrix::SetCol(matrix &A, int startI, int startJ, int countI){
	
	for (int i = 0; i < countI; i++){
			
		this->mat[startI*(this->ncols)+startJ+i*this->ncols] = A(i,0);//A(i+startI,startJ);
		
	}
	
}

void matrix::SetArray(matrix & A, vector<int> IStart, vector<int> JStart, int countI, int countJ){
	
	int startI = IStart[0];
	int startIA = IStart[1];
	
	int startJ = JStart[0];
	int startJA = JStart[1];

	for (int i = 0; i < countI; i++){
		for (int j = 0; j < countJ; j++){
			
			this->mat[i*(this->ncols)+j+startI*this->ncols+startJ] = A(i+startIA,j+startJA);
			
		}
	}
	
	
	
	
}

//Sum elements of array
double SumMat(matrix x, int M, int N){
	
	double sum = 0.0;
	
	for (int i = 0; i < M; i++){
		for (int j = 0; j < N; j++){
		
			sum += x(i,j);
			
		}
	}
	
	return sum;
	
}

//trig forcing function

void forcing(matrix & f, int M, int N, double L, matrix & hx_vec, matrix & hy_vec){
	
	for (int i = 0; i < M; i++){
		for (int j = 0; j < N; j++){
			
			f(i,j) = 0.1*sin(2.0*M_PI*hx_vec(j)/L)*cos(2.0*M_PI*hy_vec(i)/L);
			
		}
	}	
}

//trig initial conditions

void InitialConditions(matrix & u, int M, int N, double L, matrix & hx_vec, matrix & hy_vec){
	
	for (int i = 0; i < M; i++){
		for (int j = 0; j < N; j++){
			
			u(i,j) = 2.0*sin(0.5*3.0*M_PI*hx_vec(j)/L)*cos(6.0*M_PI*hy_vec(i)/L);
			
		}
	}	
}


//Matrix multiplication
void MatMult(matrix & p, const matrix & q, int M, int N, double hx, double hy, double dt){
	
	for (int i = 1; i < M+1; i++){
		for (int j = 1; j < N+1; j++){
			
			p(i-1,j-1) = q(i,j)-(dt/pow(hx,2))*(-2.0*q(i,j)+q(i,j+1)+q(i,j-1))
							-(dt/pow(hy,2))*(-2.0*q(i,j)+q(i-1,j)+q(i+1,j));
				
			}
		}
		
}


