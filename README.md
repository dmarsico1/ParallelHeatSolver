# ParallelHeatSolver

Parallel C++ code to solve the forced heat equation with periodic boundary conditions.  The MPI cartesian topology routines are used to decompose the domain into rectangles.

To compile: make heat

To run: mpirun -n (# of processors) heat

The file "input.txt" contains an example of how to change the default values of the parameters.  In this case, it's the number of cells in each dimension, as well as the time step.
