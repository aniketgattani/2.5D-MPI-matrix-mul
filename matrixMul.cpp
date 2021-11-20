#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include "timer.h"

void matrix_multiply(int *A, int *B, int *C, int n){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				C[i*n + j] += A[i*n + k] * B[k*n + j];
			}
		}
	}
}

int computeRank(int i, int j, int k, int l){
	return k*l*l + i*l + j;
}

void initialise_buffers(int ***buffA, int b){
	*buffA = (int **)malloc (2 * sizeof(int*));
	(*buffA)[0] = (int *)malloc (b * sizeof(int));
	(*buffA)[1] = (int *)malloc (b * sizeof(int));			
}
void initialise_matrices(int ***buffA, int ***buffB, int ***buffC, int n, int rank, int k, int c){
	
	initialise_buffers(buffA, n*n);
	initialise_buffers(buffB, n*n);
	initialise_buffers(buffC, n*n);

	srand(rank/c+1);
	
	for(int i=0; i<n; i++){

		for(int j=0; j<n; j++){
			if(k==0){
				(*buffA[0])[i*n + j] = rand()%50;
				(*buffB[0])[i*n + j] = rand()%50;
			}
			else{
				(*buffA[0])[i*n + j] = 0;
				(*buffB[0])[i*n + j] = 0;
			}
			(*buffC)[0][i*n + j] = 0;	
		}
	}	
}

void print_matrix(int *A, int n){
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

int cannon(int n, int p, int c, int i, int j, int k, int rank, int *C, MPI_Comm *icomm, MPI_Comm *jcomm){
	if(k!=0) return 1;
	int l = sqrt(p/c);
	int size = n/l;
	int b = size*size;

	int **buffA, **buffB, **buffC;
	
	MPI_Status status[4];
	MPI_Request reqs[4];

	initialise_matrices(&buffA, &buffB, &buffC, size, rank, k, c);
	
	int r = (l + j + i)%l;
	int s = (l + j - i)%l;
	
	int r1 = (l + i + j)%l;
	int s1 = (l + i - j)%l;


	MPI_Isend(buffA[0], b, MPI_INT, s, 0, *icomm, &reqs[0]);
	MPI_Irecv(buffA[1], b, MPI_INT, r, 0, *icomm, &reqs[1]);

	MPI_Isend(buffB[0], b, MPI_INT, s1, 0, *jcomm, &reqs[2]);
	MPI_Irecv(buffB[1], b, MPI_INT, r1, 0, *jcomm, &reqs[3]);

	MPI_Waitall(4, reqs, status);
 
	matrix_multiply(buffA[1],buffB[1],buffC[0],size);
	
	s = (l + j + 1)%l;
	s1 = (l + i + 1)%l;
	r = (l + j - 1)%l;	
	r1 = (l + i - 1)%l;

	int nrounds = l;
	
	for(int t=1; t < nrounds; t++){
		MPI_Isend(buffB[t%2], b, MPI_INT, s1, 0, *jcomm, &reqs[0]);
		MPI_Irecv(buffB[1-t%2], b, MPI_INT, r1, 0, *jcomm, &reqs[1]);
	
		MPI_Isend(buffA[t%2], b, MPI_INT, s, 0, *icomm, &reqs[2]);
		MPI_Irecv(buffA[1-t%2], b, MPI_INT, r, 0, *icomm, &reqs[3]);
		
		MPI_Waitall(4, reqs, status);
		
		matrix_multiply(buffA[1-t%2],buffB[1-t%2],buffC[0],size);	
	}		
	for(int ii=0;ii<size;ii++){
		for(int jj=0;jj<size;jj++){
			if(C[ii*size + jj]!=buffC[0][ii*size + jj]){
				return 0;
			}
		}
	}
	return 1;

}
/* 
	l is number of processors along length of the front face of grid
	matSize is matrix size for each processor to process
	blockSize is the no of elements in matrix for each processor
*/ 
void calculateGridDimensions(int *l, int *matSize, int *blockSize, int p, int c, int n){
	*l = sqrt(p/c);
	*matSize = n/l;
	*blockSize = matSize*matSize;
}


int main(int argc, char** argv) {
	
	MPI_Init(&argc, &argv);

	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// parse command line args
	int n = atoi(argv[1]), p = atoi(argv[2]), c = atoi(argv[3]), verbose = 0, diff=0;
	if(argc > 4) verbose = atoi(argv[4]);
	if(argc > 5) diff = atoi(argv[5]);


	int l, matSize, blockSize;
	calculateGridDimensions(&l, &matSize, &blockSize, p, c, n);
	
	int dims[3] = {l, l, c}, periodicity[3] = {0,0,0}, coords[3];

	MPI_Comm cartComm, kcomm, icom, jcomm;	
	
	int cartRank;

	// create MPI cartesian grid and find coords & rank of current MPI process
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periodicity, true, &cartComm);	
	MPI_Comm_rank(cartComm, &cartRank);
	MPI_Cart_coords(cartComm, cartRank, 3, coords);

	int i = coords[0], j = coords[1], k = coords[2];

	// construct communicators along i,j and k directions
	MPI_Comm_split(cartComm, computeRank(i,j,0,l), computeRank(i,j,k,l), &kcomm);
	MPI_Comm_split(cartComm, computeRank(i,0,k,l), computeRank(i,j,k,l), &icomm);
	MPI_Comm_split(cartComm, computeRank(0,j,k,l), computeRank(i,j,k,l), &jcomm);

	// buffers to hold A, B and C matrices. These are stored as one 
	int **buffA, **buffB, **buffC;
	
	initialise_matrices(&buffA, &buffB, &buffC, size, rank, k, c);

	if(rank == 0) timer_start();
	if(k==0 and verbose){
	
		printf("rank: %d \n", rank);
		print_matrix(buffA[0], size);
		printf("---------------------------- \n");
		print_matrix(buffB[0], size);
		printf("---------------------------- \n");
		
	}
	MPI_Bcast(buffA[0], b, MPI_INT, 0, kcomm);
	MPI_Bcast(buffB[0], b, MPI_INT, 0, kcomm);

	int r = (l + j + i - k*(l/c))%l;
	int s = (l + j - i + k*(l/c))%l;
	
	int r1 = (l + i + j - k*(l/c))%l;
	int s1 = (l + i - j + k*(l/c))%l;

	MPI_Status status[4];
	MPI_Request reqs[4];

	MPI_Isend(buffA[0], b, MPI_INT, s, 0, icomm, &reqs[0]);
	MPI_Irecv(buffA[1], b, MPI_INT, r, 0, icomm, &reqs[1]);

	MPI_Isend(buffB[0], b, MPI_INT, s1, 0, jcomm, &reqs[2]);
	MPI_Irecv(buffB[1], b, MPI_INT, r1, 0, jcomm, &reqs[3]);

	MPI_Waitall(4, reqs, status);
 
	matrix_multiply(buffA[1],buffB[1],buffC[0],size);
	
	s = (l + j + 1)%l;
	s1 = (l + i + 1)%l;
	r = (l + j - 1)%l;	
	r1 = (l + i - 1)%l;

	int nrounds = l/c;
	if(k==c-1) nrounds = l-(l/c)*(c-1);
	for(int t=1; t < nrounds; t++){
		MPI_Isend(buffB[t%2], b, MPI_INT, s1, 0, jcomm, &reqs[0]);
		MPI_Irecv(buffB[1-t%2], b, MPI_INT, r1, 0, jcomm, &reqs[1]);
	
		MPI_Isend(buffA[t%2], b, MPI_INT, s, 0, icomm, &reqs[2]);
		MPI_Irecv(buffA[1-t%2], b, MPI_INT, r, 0, icomm, &reqs[3]);
		
		MPI_Waitall(4, reqs, status);
		
		matrix_multiply(buffA[1-t%2],buffB[1-t%2],buffC[0],size);	
	}		

	MPI_Reduce(buffC[0], buffC[1], b, MPI_INT, MPI_SUM, 0, kcomm);
	
	if(k==0 and verbose){
		print_matrix(buffC[1], n/l);
	}

	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0) {
		double s = timer_elapsed();
		printf("time is %f, n: %d, p: %d, c: %d \n", s,n,p,c);
	}

	if(diff){
		int res = cannon(n, p, c, i, j, k, rank, buffC[1], &icomm, &jcomm);
		int finalRes = 0;
		MPI_Reduce(&res, &finalRes, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
		if(rank==0) printf("Diff is %d \n", finalRes);
	}
/*
	int nranks;
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	printf("no of ranks: %d \n", nranks);
*/	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}