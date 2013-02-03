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
void resetMatrix(double** matrix);
void printInOrder(int rank, int nprocs, int nbodies, double** s, double** v, double* m);


void readnbody(double** s, double** v, double* m, int n) {
	int myrank, nprocs, nbodies;
	int i, j, cpu;
	double* tmp;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
			
	nbodies = n/nprocs;	
	tmp = (double *)malloc(sizeof(double)*7*nbodies);


	// Node 0 reads the file and distributes the data to the right proc 
	if (myrank == 0) {
		for (cpu = 0; cpu < nprocs; cpu++) {
			
			//for each CPU we want to send the right data
			for (j = 0; j < nbodies; j++) {

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
			MPI_Send(&tmp[0], nbodies*7, MPI_DOUBLE, cpu, 0, MPI_COMM_WORLD);	
		}
	}

	MPI_Recv(&tmp[0], nbodies*7, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	for (i = 0; i < nbodies; i++) {
		s[i][0] = tmp[i*7];
		s[i][1] = tmp[i*7+1];
		s[i][2] = tmp[i*7+2];

		v[i][0] = tmp[i*7+3];
		v[i][1] = tmp[i*7+4];
		v[i][2] = tmp[i*7+5];

		m[i] = tmp[i*7+6];
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
	int myrank, nprocs;
	int i,j,k, size, l, p;
	double* distance;
	double* currentplanets;
	double** acceleration;
	double r, G,f;
	MPI_Status status;	
	

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	size = n / nprocs;
	distance = (double *)malloc(sizeof(double) * 3);	
	acceleration = (double **)malloc(sizeof(double *) * size);	
	// array med planeter x,y,z,masse	
	currentplanets = (double *)malloc(sizeof(double) * 4 * size);	
	
	

	for (i = 0; i < size; i++) {  
		acceleration[i] = (double*)malloc(sizeof(double) * 3);
	}	

	G = 6.674e-11;

	for(i = 0; i < iter; i++){				//for loop over iterasjoner
		
		resetMatrix(acceleration);
		
		for (k = 0; k < size; k++)
			for (j = 0; j < 3; j++)
				printf("A[%d][%d]: %1.4e\n", k, j, acceleration[k][j]);
		

		for (k = 0; k < size; k++) {  
			for(j = 0; j < 4; j++){
				currentplanets[j+k*4] = (j == 3) ? m[k] : s[k][j];
			}
			
		}		
		for(p = 0; p < nprocs ; p++){


			for(j = 0; j < size;j++){			//for loop for en spesefikk planet
				for(k = 0; k< size;k++){		//alle planeter for en spesifikk
					for(l = 0; l < 3;l++){
						distance[l] = s[j][l] - currentplanets[k*4 +l];	
					}
					r = norm(distance);
					if(r==0) continue;
					f = (G * m[j] * currentplanets[k*4+3] ) / (r*r);	
					for(l = 0; l < 3;l++){
						distance[l] = ((distance[l] * f )/ r) / m[j];
						acceleration[j][l] = acceleration[j][l] - distance[l];
					}

				}
					
			}
			
			if (p == (nprocs-1)) continue;

			if(myrank % 2 == 0){
				//MPI Send first, then recieve
				MPI_Send(&currentplanets[0], size*4, MPI_DOUBLE, (myrank+1)%nprocs, 0, MPI_COMM_WORLD);
				MPI_Recv(&currentplanets[0], size*4, MPI_DOUBLE, (myrank-1)%nprocs, 0, MPI_COMM_WORLD, &status);
			}else{					
				double* tmp;
				tmp  = (double *)malloc(sizeof(double)*size*4);

				//MPI Recieve first, then send
				MPI_Recv(&tmp[0], size*4, MPI_DOUBLE, (myrank -1)%nprocs, 0, MPI_COMM_WORLD, &status);
				MPI_Send(&currentplanets[0], size*4, MPI_DOUBLE, (myrank +1)%nprocs, 0, MPI_COMM_WORLD);			
				
				for (j = 0; j < size*4; j++) {
					currentplanets[j] = tmp[j];
				}

				free(tmp);
			}
		}
		for(j=0; j < size;j++){
			for(k=0;k<3;k++){
				v[j][k] = v[j][k] + timestep * acceleration[j][k];
				s[j][k] = s[j][k] + timestep * v[j][k];
			}
		}			
	}

	
	free(distance);
	free(currentplanets);
	free(acceleration);

	printInOrder(myrank, nprocs, size, s, v, m);
	
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

void resetMatrix(double** matrix) {
	int i, j, len, len0;
	len = sizeof(matrix)/sizeof(matrix[0]);
	len0 = sizeof(matrix[0])/sizeof(matrix[0][0]);
	printf("Len: %d\nLen0: %d\n", len, len0);
	for (i = 0; i < 2; i++)
		for (j = 0; j < 3; j++)
			matrix[i][j] = 0;
}
void printInOrder(int rank, int nprocs, int nbodies, double** s, double** v, double* m) {	
	int p, i;	
	for (p = 0; p < nprocs; p++) {
		if (rank == p) {
			for (i = 0; i < nbodies; i++) {
				fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
