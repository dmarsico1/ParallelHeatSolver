#include "mpi_interface.h" 

using namespace std;


GridTopology::GridTopology(int order){

	reorder = order;
	
}

void GridTopology::CreateTopology(int Nx, int Ny, double hy, double hx, double Ly, double Lx){
	
	periods.resize(2);
	coords.resize(2);
	dims.resize(2);
	
	periods[0] = 1;
	periods[1] = 1;
	
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	
	MPI_Dims_create(nproc,2,&dims[0]);
	
	MPI_Cart_create(MPI_COMM_WORLD,2,&dims[0],&periods[0],reorder,&comm2d);
	
	MPI_Cart_get(comm2d,2,&dims[0],&periods[0],&coords[0]);
	
	//Get number of cells in each process, and define local coordinates
	
	if(coords[0] == 0){
	
		y_dim = Ny - (dims[0]-1)*(Ny/dims[0]);
		
		hy_vec_local.allocate(y_dim);
		
		std::cout << hy_vec_local.GetRows() << " " << hy_vec_local.GetCols() << std::endl;
		
		for (int i = 0; i < y_dim;  i++){
		
			hy_vec_local(i) = -Ly+(static_cast<double>(i)+0.5)*hy;
		
		}
		
	}
	else{
		
		y_dim = Ny/dims[0];
		
		hy_vec_local.allocate(y_dim);
		
		for (int i = 0; i < y_dim; i++){
		
			int offset = (Ny - (dims[0]-1)*(Ny/dims[0])) + y_dim*(coords[0]-1);
			
			hy_vec_local(i) = -Ly+static_cast<double>(offset)*hy + (static_cast<double>(i)+0.5)*hy;
		
		}
		
	}
	
	if(coords[1] == 0){
	
		x_dim = Nx - (dims[1]-1)*(Nx/dims[1]);
		
		hx_vec_local.allocate(x_dim);
		
		for (int i = 0; i < x_dim; i++){
		
			hx_vec_local(i) = -Lx+(static_cast<double>(i)+0.5)*hx;
		
		}
		
	}
	else{
		
		x_dim = Nx/dims[1];
		
		hx_vec_local.allocate(x_dim);
		
		for (int i = 0; i < x_dim; i++){
		
			int offset = (Nx - (dims[1]-1)*(Nx/dims[1])) + x_dim*(coords[1]-1);
			
			hx_vec_local(i) = -Lx+static_cast<double>(offset)*hx + (static_cast<double>(i)+0.5)*hx;
		
		}
		
	}
	
	//Get current process id
	MPI_Comm_rank(comm2d,&myid);
	
	//Get neighboring process ids
	
	MPI_Cart_shift(comm2d,0,1,&lower,&upper);
	
	MPI_Cart_shift(comm2d,1,1,&left,&right);
	
	//Allocate array for conjugate gradient
	p.allocate(y_dim+2,x_dim+2);

}


void GridTopology::GetDisplacements(int Ny, int Nx){
	
	local_disp = 0;
	
	if(coords[1] > 0){
	
		local_disp = (Nx - (dims[1] - 1)*(Nx/dims[1])) + (coords[1] - 1)*(Nx/dims[1]);
	
	}
	
	if(coords[0] > 0){
	
		local_disp = local_disp + (Ny - (dims[0] - 1)*(Ny/dims[0]))*Nx + (coords[0] - 1)*(Ny/dims[0])*Nx;
	
	}
	
	if(myid == 0){
		
		disps.resize(nproc);
		counts.resize(nproc);
		counts.assign(nproc,1);
		
	}
	
	MPI_Gather(&local_disp,1,MPI_INT,&disps[0],1,MPI_INT,0,comm2d);

}


void GridTopology::GetFullArray(matrix & FullArray, matrix & SubArray, int Ny, int Nx){
	
	int begin = 0;
	int extent = sizeof(double);
	
	MPI_Datatype subarray;
	MPI_Datatype resized_subarray;
	vector<int> dims_full_array(2);
	vector<int> dims_subarray(2);
	vector<int> starts(2);
	
	starts[0] = 0;
	starts[1] = 0;
	
	dims_full_array[0] = Ny;
	dims_full_array[1] = Nx;
	
	dims_subarray[0] = y_dim;
	dims_subarray[1] = x_dim;
	
	int dims = 2;
	
	MPI_Type_create_subarray(dims,&dims_full_array[0],&dims_subarray[0],&starts[0],MPI_ORDER_C,MPI_DOUBLE,&subarray);
	MPI_Type_commit(&subarray);
	
	MPI_Type_create_resized(subarray, begin, extent, &resized_subarray);
	MPI_Type_commit(&resized_subarray);
	
	MPI_Gatherv(&SubArray(0,0),x_dim*y_dim,MPI_DOUBLE,&FullArray(0,0),&counts[0],&disps[0],
				resized_subarray,0,comm2d);
	
}


void Neighbors(matrix & A, GridTopology & Grid){
	
	int x_dim = Grid.x_dim;
	int y_dim = Grid.y_dim;
	
	MPI_Status status;
	
	//Send upper and lower ghost cells
	MPI_Sendrecv(&A(y_dim,1),x_dim,MPI_DOUBLE,Grid.upper,0,&A(0,1),x_dim,MPI_DOUBLE,Grid.lower,0,Grid.comm2d,&status);
	MPI_Sendrecv(&A(1,1),x_dim,MPI_DOUBLE,Grid.lower,0,&A(y_dim+1,1),x_dim,MPI_DOUBLE,Grid.upper,0,Grid.comm2d,&status);	

	//Left and right ghost cells
	MPI_Datatype Column;
	MPI_Type_vector(Grid.y_dim,1,Grid.x_dim+2,MPI_DOUBLE,&Column);
	MPI_Type_commit(&Column);
	
	MPI_Sendrecv(&A(1,1),1,Column,Grid.left,1,&A(1,Grid.x_dim+1),1,Column,Grid.right,1,Grid.comm2d,&status);
	
	MPI_Type_vector(Grid.y_dim,1,Grid.x_dim+2,MPI_DOUBLE,&Column);
	MPI_Type_commit(&Column);
	
	MPI_Sendrecv(&A(1,Grid.x_dim),1,Column,Grid.right,1,&A(1,0),1,Column,Grid.left,1,Grid.comm2d,&status);
	
}

void ConjGrad(matrix & x, matrix & f, int M, int N, double hx, double hy, double dt, GridTopology & Grid){
	
	matrix r(M,N);
	////matrix p(M+2,N+2);
	matrix q(M,N);
	
	x = 0.0;
	
	r = f;

	////p = r;
	Grid.p = 0.0;
	
	vector<int> startsI(2);
	vector<int> startsJ(2);
	
	startsI[0] = 1; startsI[1] = 0;
	startsJ[0] = 1; startsJ[1] = 0;
	
	Grid.p.SetArray(r, startsI, startsJ, M, N);
	
	Grid.rold = SumMat(r*r,M,N);
	
	MPI_Allreduce(&Grid.rold,&Grid.rold_tot,1,MPI_DOUBLE,MPI_SUM,Grid.comm2d);
	
	std::cout << "rold is " << Grid.rold << std::endl;
	
	double tol = 0.001;
	
	int k=1;
	
	for (k = 0; k < pow(M,2); k++){
		
		Neighbors(Grid.p,Grid);
		
		//The commented block below was for zero Dirichlet boundary conditions.  
		//We're using periodic BCs in both directions
		
		//matrix Top = (Grid.p.GetSubArray(1,1,1,N))*-1.0;
		//matrix Bottom = (Grid.p.GetSubArray(M,1,1,N))*-1.0;
		//matrix Left = (Grid.p.GetSubArray(1,1,M,1))*-1.0;
		//matrix Right = (Grid.p.GetSubArray(1,N,M,1))*-1.0;
	
		//Grid.p.SetRow(Top,0,1,N);
		//Grid.p.SetRow(Bottom,M+1,1,N);
		//Grid.p.SetCol(Left,1,0,M);
		//Grid.p.SetCol(Right,1,N+1,M);
		
		MatMult(q,Grid.p,M,N,hx,hy,dt);
		
		matrix p_sub = Grid.p.GetSubArray(1,1,M,N);
		
		Grid.gamma = SumMat(q*p_sub,M,N);
		
		MPI_Allreduce(&Grid.gamma,&Grid.gamma_tot,1,MPI_DOUBLE,MPI_SUM,Grid.comm2d);
		
		Grid.alpha_tot = Grid.rold_tot/Grid.gamma_tot;
		
		matrix p_alpha = Grid.p*Grid.alpha_tot;
		
		vector<int> startI(3);
		vector<int> startJ(3);
		
		startI[0] = 0; startI[1] = 0; startI[2] = 1; 
		startJ[0] = 0; startJ[1] = 0; startJ[2] = 1;
		
		x.Add(x,p_alpha,startI,startJ,M,N);
		
		r -= q*Grid.alpha_tot;
		
		//r.Subtract(r,q,startI,startJ,M,N);
		
		Grid.rnew = SumMat(r*r,M,N);
		
		MPI_Allreduce(&Grid.rnew,&Grid.rnew_tot,1,MPI_DOUBLE,MPI_SUM,Grid.comm2d);
		
		std::cout << "rnew_tot " << Grid.rnew_tot << " " << k << std::endl;
		
		if (sqrt(Grid.rnew_tot) < tol){
			break;
		}
		
		matrix p_r = Grid.p*(Grid.rnew_tot/Grid.rold_tot);
		
		startI[0] = 1; startI[1] = 0; startI[2] = 1; 
		startJ[0] = 1; startJ[1] = 0; startJ[2] = 1;
		
		Grid.p.Add(r,p_r,startI,startJ,M,N);
		
		Grid.rold_tot=Grid.rnew_tot;
		
	}
	
} 


void GetGridDims(int Ny, int & y_dim, int nproc, int myid){

	if(myid == 0){

		y_dim = Ny - (nproc-1)*(Ny/nproc);

	}
	else{

		y_dim = Ny/nproc;
		
	}

}

void GetNeighbors(int & lower, int & upper, int nproc, int myid){

	if(nproc == 1){
	
		lower = 0;
		upper = 0;
	
	}
	else{
	
		if(myid == 0){
		
			lower = nproc - 1;
			upper = myid+1;
		
		}
		else if(myid == nproc-1){
		
			lower = myid-1;
			upper = 0;
		
		}
		else{
		
			lower = myid - 1;
			upper = myid+1;
		
		}
	
	}



}
