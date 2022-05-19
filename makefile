.PHONY: clean

heat: matrix.o mpi_interface.o grid.o nc_data.o heat.o
	mpic++ -o heat matrix.o mpi_interface.o grid.o nc_data.o heat.o -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdf_c++4

matrix.o: matrix.cpp
	mpic++ -c matrix.cpp

mpi_interface.o: mpi_interface.cpp
	mpic++ -c mpi_interface.cpp

grid.o: grid.cpp
	mpic++ -c grid.cpp	

nc_data.o: nc_data.cpp
	mpic++ -c nc_data.cpp

heat.o: heat.cpp
	mpic++ -c heat.cpp

clean:
		rm -f matrix.o grid.o heat.o mpi_interface.o nc_data.o heat
