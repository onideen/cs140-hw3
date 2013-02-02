/*
Assignment 3 
Team Member 1 : Arne Bjune
Team Member 2 : Vegar Engen
*/

#include "nBody.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
double norm(double * x);


void readnbody(double** s, double** v, double* m, int n) {
	int myrank, nprocs, nbody;
	int i, j, cpu;
	double* tmp;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
			
	nbody = n/nprocs;	
	tmp = (double *)malloc(sizeof(double)*7*nbody);


	// Node 0 reads the file and distributes the data to the right proc 
	if (myrank == 0) {
		for (cpu = 0; cpu < nprocs; cpu++) {
			
			//for each CPU we want to send the right data. !!! START OPP IGJEN HER
			for (j = 0; j < nbody; j++) {

				double x, y, z, vx, vy, vz, ma;

				int result = scanf(INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &ma);
				if (result != 7) {
					fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
					exit(0);
				} else {
					 
					tmp[j*7] = x;
					tmp[j*7+1] = y;
					tmp[j*7+2] = z;
					tmp[j*7+3] = vx;
					tmp[j*7+4] = vy;
					tmp[j*7+5] = vz;
					tmp[j*7+6] = ma;					
				}
			}
			MPI_Send(&tmp[0], nbody*7, MPI_DOUBLE, cpu, 0, MPI_COMM_WORLD);
			
		}
	
	}

	MPI_Recv(&tmp[0], nbody*7, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	for (i = 0; i < nbody; i++) {
		s[i][0] = tmp[i*7];
		s[i][1] = tmp[i*7+1];
		s[i][2] = tmp[i*7+2];

		v[i][0] = tmp[i*7+3];
		v[i][1] = tmp[i*7+4];
		v[i][2] = tmp[i*7+5];

		m[i] = tmp[i*7+6];
		printf("CPU %d: ", myrank);
		printf(OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
	}
	free(tmp);
	 
}

void gennbody(double** s, double** v, double* m, int n) {
	int i, j;
	double dist, theta;
	int myrank, nprocs, nbody;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	nbody = n/nprocs;

	srand(time(NULL)*myrank);
	for (i = 0; i < nbody; i++) {
		m[i] = 1e30 * (float)rand()/RAND_MAX;
		dist = 0.5e13 * (float)rand()/RAND_MAX;
		theta = 2*M_PI*(float)rand()/RAND_MAX;

		s[i][0] = dist*cos(theta);
		s[i][1] = dist*sin(theta);
		s[i][2] = 1e11*((float)rand()/RAND_MAX-.5);

		for (j = 0; j < 3; j++) {
			v[i][j] = 0;
		}	
	}

}

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	int myrank;
	int nprocs;
	int i,j,k, size, l;
	double* distance;
	double** acceleration;
	double r, G,f;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	size = n / nprocs;
	distance = (double *)malloc(sizeof(double) * 3);	
	acceleration = (double **)malloc(sizeof(double *) * size);	

	for (i = 0; i < size; i++) {
		acceleration[i] = (double*)malloc(sizeof(double) * 3);

		for(j = 0; j < 3; j++) {
			acceleration[i][j] = 0;
		}
	}	

	G = 6.674e-11;
	

	// s = posisjon
	// v = hastighet
	// m = masse
	// n = antall planeter
	// iter = antall steg
	// timestep = tid mellom steg
	// distance = avstandsvektor mellom 2 planeter


	for(i = 0; i < iter; i++){				//for loop over iterasjoner
		for(j = 0; j < size;j++){			//for loop for en spesefikk planet
			//printf("Masser, m[0] %1.4e, m[1] %1.4e, m[2] %1.4e, m[3] %1.4e \n",m[0],m[1],m[2],m[3]);	

			for(k = 0; k < 3; k++){
				acceleration[j][k] = 0;	
			}	
			
			for(k = 0; k< size;k++){		//alle planeter for en spesifikk
				for(l = 0; l < 3;l++){
					distance[l] = s[j][l] - s[k][l];
					
				}
				r = norm(distance);
				//printf("Avstand mellom %i og %i er %1.4e \n",j,k,r);
				if(r==0) continue;
				//f = 2.0;
				//printf("F er: %1.4e , G er: %1.4e , m[%i]: %1.4e, m[%i]: %1.4e \n", f, G,j,m[j],k,m[k]);
				f = (G * m[j] * m[k] ) / (r*r);
				
				for(l = 0; l < 3;l++){
					distance[l] = ((distance[l] * f )/ r) / m[j];
					acceleration[j][l] = acceleration[j][l] - distance[l];
				}

			}
			for(k=0;k<3;k++){
				//printf("Setter hastighet til %1.4e \n", timestep * acceleration[j][k]);
				v[j][k] = v[j][k] + timestep * acceleration[j][k];
				s[j][k] = s[j][k] + timestep * v[j][k];
			}		
		}
		
	}


	
	// This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
	// in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
	// sure you dont add extra debugging statements in stderr.

	if (myrank == 0) {
		for (i = 0; i < n / nprocs; i++) {
			fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
		}
	}
}

double norm(double * x){
	int i;
	double sum = 0;

	for(i = 0; i<3;i++){
		sum = sum + x[i]*x[i];
	}
	sum = sqrt(sum);

	return sum;
}

