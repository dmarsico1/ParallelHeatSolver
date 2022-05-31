# ParallelHeatSolver

Parallel C++ code to solve the forced heat equation in two dimensions with periodic boundary conditions on a cell-centered rectangular grid.  The MPI cartesian topology routines are used to decompose the domain into rectangles, and a parallel conjugate gradient method is used to solve the pressure equation.

To compile: make heat

To run: mpirun -n (# of processors) heat

The file "input.txt" contains an example of how to change the default values of the parameters.  In this case, it's the number of cells in each dimension, as well as the time step.
