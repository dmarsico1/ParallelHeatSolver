#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <math.h>
#include <vector>

using namespace std;


class matrix
{
	
	public:
		
		matrix();
		matrix(int dim);
		matrix(int dim, double a);
		matrix(int rows, int cols);
		matrix(int rows, int cols, double a);
		~matrix();
		
		void allocate(int rows);
		void allocate(int rows, int cols);
		
		matrix & operator= (double a);
		matrix & operator= (const matrix &A);
		
		int GetRows();
		int GetCols();
		
		double & operator() (int irow); //subscript operator (vector)
		const double & operator() (int irow) const; //subscript operator (vector)
		
		double & operator() (int irow, int icol); //subscript operator
		const double & operator() (int irow, int icol) const; //subscript operator
		
		matrix & operator+= (const matrix &A);
		matrix & operator-= (const matrix &A);
		matrix & operator*= (double a); //scalar multiplication
		matrix & operator/= (double a); //scalar division
		
		matrix operator+ (const matrix &A); //addition
		matrix operator- (const matrix &A); //subtraction
		matrix operator* (const double a); //scalar multipliation
		matrix operator* (const matrix &A); //element-wise multiplication
		matrix operator/ (const double a);
		
		void Add(matrix & A, matrix & B, vector<int> IStart, vector<int> JStart, int countI, int countJ);
		void Subtract(matrix & A, matrix & B, vector<int> IStart, vector<int> JStart, int countI, int countJ);
		matrix GetSubArray(int startI, int startJ, int countI, int countJ);
		void SetArray(matrix & A, vector<int> IStart, vector<int> JStart, int countI, int countJ);
		void SetCol(matrix &A, int startI, int startJ, int countI);
		void SetRow(matrix &A, int startI, int startJ, int countI);
		
			
	private:
		
		int nrows;
		int ncols;
		double *mat;
	
};


//Non member functions

double SumMat(matrix x, int M, int N);

void MatMult(matrix & p, const matrix & q, int M, int N, double hx,double hy, double dt);

void forcing(matrix & f, int M, int N, double L, matrix & hx_vec, matrix & hy_vec);

void InitialConditions(matrix & u, int M, int N, double L, matrix & hx_vec, matrix & hy_vec);

//double Add(matrix &x, int M, int N);

#endif
