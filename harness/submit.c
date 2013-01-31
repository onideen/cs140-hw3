/*
Assignment 3 
Team Member 1 : Arne Bjune
Team Member 2 : Vegar Engen
*/

#include "nBody.h"

void readnbody(double** s, double** v, double* m, int n) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	// This is an example of reading the body parameters from the input file. 
	if (myrank == 0) {
		for (i = 0; i < n; i++) {
			double x, y, z, vx, vy, vz, m;

			int result = scanf(INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &m);
			if (result != 7) {
				fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
				exit(0);
			}
		}
	}
}

void gennbody(double** s, double** v, double* m, int n) {
	int i, j;
	double dist, theta;
	for (i = 0; i < n; i++) {
		m[i] = 1e30 * rand();
		dist = 0.5e13 * rand();

		s[i][0] = dist;

		for (j = 0; j < 3; j++) {
			v[i][j] = 0;
		}
	}

	printf("Generate nBody initial condition here.\n");
}

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	// This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
	// in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
	// sure you dont add extra debugging statements in stderr.

	if (myrank == 0) {
		for (i = 0; i < n / nprocs; i++) {
			fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
		}
	}
}

