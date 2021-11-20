#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include "timer.h"

void matrixMultiply(int *A, int *B, int *C, int n){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				C[i*n + j] += A[i*n + k] * B[k*n + j];
			}
		}
	}
}

// to compute MPI rank in the communicator based on grid coordinates  
int computeRank(int i, int j, int k, int l){
	return k*l*l + i*l + j;
}

void initialiseBuffers(int **buffA, int b){
	(*buffA) = (int *)malloc (b * sizeof(int));			
}

// initialize matrices. Random values for k==0 and 0 for k!=0
void initialiseMatrices(int **buffA, int **buffB, int **buffC, int n, int rank, int k, int c){
	
	initialiseBuffers(buffA, n*n);
	initialiseBuffers(buffB, n*n);
	initialiseBuffers(buffC, n*n);

	srand(rank/c+1);
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(k==0){
				(*buffA)[i*n + j] = rand()%50;
				(*buffB)[i*n + j] = rand()%50;
			}
			else{
				(*buffA)[i*n + j] = 0;
				(*buffB)[i*n + j] = 0;
			}
			(*buffC)[i*n + j] = 0;	
		}
	}	
}

// printing matrix
void printMatrix(int *A, int n){
	printf("[");
	for(int i=0; i<n; i++){
		printf("[");
		for(int j=0; j<n-1; j++){
			printf("%d, ", A[i*n + j]);
		}
		printf("%d]", A[i*n + n-1]);
		if(i!=n-1) printf(",\n");
	}
	printf("]\n");
}

/* 
	l is number of processors along length of the front face of grid
	matSize is matrix size for each processor to process
	blockSize is the no of elements in matrix for each processor
*/ 
void calculateGridDimensions(int *l, int *matSize, int *blockSize, int p, int c, int n){
	*l = sqrt(p/c);
	*matSize = n/(*l);
	*blockSize = (*matSize)*(*matSize);
}

// Verification using Cannon's blocked algo 
int cannon(int n, int p, int c, int i, int j, int k, int rank, int *C, MPI_Comm *icomm, MPI_Comm *jcomm){

	if(k!=0) return 1;
	
	int l, matSize, blockSize;
	calculateGridDimensions(&l, &matSize, &blockSize, p, c, n);

	// buffers to hold A, B and C matrices. These are stored as 1-D arrays
	int *buffA, *buffB, *buffC;

	initialiseMatrices(&buffA, &buffB, &buffC, matSize, rank, k, c);

	MPI_Status status[2];
	
	// Cannon's Initial shift
	int r = (l + j + i)%l;
	int s = (l + j - i)%l;	
	int r1 = (l + i + j)%l;
	int s1 = (l + i - j)%l;
	MPI_Sendrecv_replace(buffA, blockSize, MPI_INT, s, 0, r, 0, *icomm, &status[0]);
	MPI_Sendrecv_replace(buffB, blockSize, MPI_INT, s1, 1, r1, 1, *jcomm, &status[1]); 
	matrixMultiply(buffA, buffB, buffC, matSize);
	
	
	// SUMMA block iterations
	s = (l + j + 1)%l;
	s1 = (l + i + 1)%l;
	r = (l + j - 1)%l;	
	r1 = (l + i - 1)%l;

	int nrounds = l;
	for(int t=1; t < nrounds; t++){
		MPI_Sendrecv_replace(buffA, blockSize, MPI_INT, s, 0, r, 0, *icomm, &status[0]);		
		MPI_Sendrecv_replace(buffB, blockSize, MPI_INT, s1, 1, r1, 1, *jcomm, &status[1]);
		matrixMultiply(buffA,buffB,buffC,matSize);	
	}

	// check if the calculations match the 2.5D calculations
	for(int ii=0;ii<matSize;ii++){
		for(int jj=0;jj<matSize;jj++){
			if(C[ii*matSize + jj]!=buffC[ii*matSize + jj]){
				return 0;
			}
		}
	}
	return 1;
}

int main(int argc, char** argv) {
	
	MPI_Init(&argc, &argv);

	int rank, cartRank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// parse command line args
	int n = atoi(argv[1]), p = atoi(argv[2]), c = atoi(argv[3]), verbose = 0, diff=0;
	if(argc > 4) verbose = atoi(argv[4]);
	if(argc > 5) diff = atoi(argv[5]);

	int l, matSize, blockSize;
	calculateGridDimensions(&l, &matSize, &blockSize, p, c, n);

	int dims[3] = {l, l, c}, periodicity[3] = {0,0,0}, coords[3];

	MPI_Comm cartComm, kcomm, icomm, jcomm;	

	// create MPI cartesian grid and find coords & rank of current MPI process
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periodicity, true, &cartComm);	
	MPI_Comm_rank(cartComm, &cartRank);
	MPI_Cart_coords(cartComm, cartRank, 3, coords);
	int i = coords[0], j = coords[1], k = coords[2];

	// construct communicators along i,j and k directions
	MPI_Comm_split(cartComm, computeRank(i,j,0,l), computeRank(i,j,k,l), &kcomm);
	MPI_Comm_split(cartComm, computeRank(i,0,k,l), computeRank(i,j,k,l), &icomm);
	MPI_Comm_split(cartComm, computeRank(0,j,k,l), computeRank(i,j,k,l), &jcomm);

	// buffers to hold A, B and C matrices. These are stored as 1-D arrays
	int *buffA, *buffB, *buffC;
	
	initialiseMatrices(&buffA, &buffB, &buffC, matSize, rank, k, c);

	double totalTime, waitTime, t;
	totalTime = MPI_Wtime();

	// printing matrices
	if(k==0 and verbose){
		printf("rank: %d \n", rank);
		printMatrix(buffA, matSize);
		printf("---------------------------- \n");
		printMatrix(buffB, matSize);
		printf("---------------------------- \n");
	}

	MPI_Status status[2];

	// first broadcast along k direction
	MPI_Bcast(buffA, blockSize, MPI_INT, 0, kcomm);
	MPI_Bcast(buffB, blockSize, MPI_INT, 0, kcomm);

	// Cannon's Initial shift
	int r = (l + j + i - k*(l/c))%l;
	int s = (l + j - i + k*(l/c))%l;
	int r1 = (l + i + j - k*(l/c))%l;
	int s1 = (l + i - j + k*(l/c))%l;

	t = MPI_Wtime();
	MPI_Sendrecv_replace(buffA, blockSize, MPI_INT, s, 0, r, 0, icomm, &status[0]);
	MPI_Sendrecv_replace(buffB, blockSize, MPI_INT, s1, 1, r1, 1, jcomm, &status[1]);
	t = MPI_Wtime() - t;
	waitTime += t;

	matrixMultiply(buffA,buffB,buffC,matSize);
	
	// 1/c of SUMMA block iterations
	s = (l + j + 1)%l;
	s1 = (l + i + 1)%l;
	r = (l + j - 1)%l;	
	r1 = (l + i - 1)%l;

	// this is to ensure correct calculation for k=c-1 when when l%c!=0 
	int nrounds = l/c;
	if(k==c-1) nrounds = l-(l/c)*(c-1);
	for(int t=1; t < nrounds; t++){
		t = MPI_Wtime();
		MPI_Sendrecv_replace(buffA, blockSize, MPI_INT, s, 0, r, 0, icomm, &status[0]);		
		MPI_Sendrecv_replace(buffB, blockSize, MPI_INT, s1, 1, r1, 1, jcomm, &status[1]);
		t = MPI_Wtime() - t;
		waitTime += t;
		matrixMultiply(buffA,buffB,buffC,matSize);	
	}		
	
	// In place reduce to use the same buffer along k direction
	if(k==0) MPI_Reduce(MPI_IN_PLACE, buffC, blockSize, MPI_INT, MPI_SUM, 0, kcomm);
	else MPI_Reduce(buffC, buffC, blockSize, MPI_INT, MPI_SUM, 0, kcomm);
	
	// printing final matrix
	if(k==0 and verbose){
		printMatrix(buffC, matSize);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	totalTime = MPI_Wtime() - totalTime;
	double maxWaitTime;
	//MPI_Reduce(&waitTime, &maxWaitTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if(rank==0) {
		printf("execution time: %f, max wait time: %f, n: %d, p: %d, c: %d \n", totalTime, maxWaitTime, n,p,c);		
	}

	// if verification is to be done, call Cannon's algo, perform verification for each processor and collect result
	if(diff){
		int res = cannon(n, p, c, i, j, k, rank, buffC, &icomm, &jcomm);
		int finalRes = 0;
		MPI_Reduce(&res, &finalRes, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
		if(rank==0) printf("Diff is %d \n", finalRes);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
