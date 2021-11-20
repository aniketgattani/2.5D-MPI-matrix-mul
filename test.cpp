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

void initialise_buffers(int **buffA, int b){
	(*buffA) = (int *)malloc (b * sizeof(int));			
}
void initialise_matrices(int **buffA, int **buffB, int **buffC, int n, int rank, int k, int c){
	
	initialise_buffers(buffA, n*n);
	initialise_buffers(buffB, n*n);
	initialise_buffers(buffC, n*n);

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
/*	if(k!=0) return 1;
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
*/	return 1;

}

int main(int argc, char** argv) {
	
	MPI_Init(&argc, &argv);

	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
	int c = atoi(argv[3]);
	int verbose = 0;
	int diff=0;
	if(argc > 4) verbose = atoi(argv[4]);
	if(argc > 5) diff = atoi(argv[5]);
	int l = sqrt(p/c);
	int size = n/l;
	int b = size*size;
	int dims[3] = {l, l, c};
	const int periodicity[3] = {0,0,0}; 
	int coords[3];

	MPI_Comm cartComm;
	MPI_Comm kcomm;
	MPI_Comm icomm;
	MPI_Comm jcomm;	
	
	int cartRank;

	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periodicity, true, &cartComm);	
	MPI_Comm_rank(cartComm, &cartRank);
	MPI_Cart_coords(cartComm, cartRank, 3, coords);

	int i = coords[0];
	int j = coords[1];
	int k = coords[2];

	MPI_Comm_split(cartComm, computeRank(i,j,0,l), computeRank(i,j,k,l), &kcomm);
	MPI_Comm_split(cartComm, computeRank(i,0,k,l), computeRank(i,j,k,l), &icomm);
	MPI_Comm_split(cartComm, computeRank(0,j,k,l), computeRank(i,j,k,l), &jcomm);

	
	int *buffA, *buffB, *buffC;
	
	initialise_matrices(&buffA, &buffB, &buffC, size, rank, k, c);

	if(rank == 0) timer_start();
	if(k==0 and verbose){
	
		printf("rank: %d \n", rank);
		print_matrix(buffA, size);
		printf("---------------------------- \n");
		print_matrix(buffB, size);
		printf("---------------------------- \n");
		
	}
	MPI_Bcast(buffA, b, MPI_INT, 0, kcomm);
	MPI_Bcast(buffB, b, MPI_INT, 0, kcomm);

	int r = (l + j + i - k*(l/c))%l;
	int s = (l + j - i + k*(l/c))%l;
	
	int r1 = (l + i + j - k*(l/c))%l;
	int s1 = (l + i - j + k*(l/c))%l;

	MPI_Status status[2];
	MPI_Request reqs[2];

	MPI_Sendrecv_replace(buffA, b, MPI_INT, s, 0, r, 0, icomm, &status[0]);
	//MPI_Irecv(buffA[1], b, MPI_INT, r, 0, icomm, &reqs[1]);

	MPI_Sendrecv_replace(buffB, b, MPI_INT, s1, 1, r1, 1, jcomm, &status[1]);
	//MPI_Irecv(buffB[1], b, MPI_INT, r1, 0, jcomm, &reqs[3]);

	//MPI_Waitall(2, reqs, status);
 
	matrix_multiply(buffA,buffB,buffC,size);
	
	s = (l + j + 1)%l;
	s1 = (l + i + 1)%l;
	r = (l + j - 1)%l;	
	r1 = (l + i - 1)%l;

	int nrounds = l/c;
	if(k==c-1) nrounds = l-(l/c)*(c-1);
	for(int t=1; t < nrounds; t++){
		MPI_Sendrecv_replace(buffB, b, MPI_INT, s1, 1, r1, 1, jcomm, &status[0]);
		//MPI_Irecv(buffB[1-t%2], b, MPI_INT, r1, 0, jcomm, &reqs[1]);
	
		MPI_Sendrecv_replace(buffA, b, MPI_INT, s, 0, r, 0, icomm, &status[1]);
		//MPI_Irecv(buffA[1-t%2], b, MPI_INT, r, 0, icomm, &reqs[3]);
		
	//	MPI_Waitall(2, reqs, status);
		
		matrix_multiply(buffA,buffB,buffC,size);	
	}		
	
	if(k==0) 
	MPI_Reduce(MPI_IN_PLACE, buffC, b, MPI_INT, MPI_SUM, 0, kcomm);
	else 
	MPI_Reduce(buffC, buffC, b, MPI_INT, MPI_SUM, 0, kcomm);
	
	if(k==0 and verbose){
		print_matrix(buffC, size);
	}

	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0) {
		double s = timer_elapsed();
		printf("time is %f, n: %d, p: %d, c: %d \n", s,n,p,c);
	}

	if(diff){
		int res = cannon(n, p, c, i, j, k, rank, buffC, &icomm, &jcomm);
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
