#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include "timer.h"
struct matrix{
	int* mat;
	int size;
};

/*void copy_matrix(int* buff, matrix &A){
	int n = A.size;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			buff[i*n+j] = A.mat[i][j];
		}
	}
}

void fill_matrix(int* buff, matrix &A){
	int n = A.size;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			A.mat[i][j] = buff[i*n+j];
		}
	}
}*/

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
void initialise_matrices(int ***buffA, int ***buffB, int ***buffC, int n, int rank, int k){
	
	initialise_buffers(buffA, n*n);
	initialise_buffers(buffB, n*n);
	initialise_buffers(buffC, n*n);
	//printf("%d %d %d %d \n", rank, n, sizeof(*buffA), sizeof(*buffA[0]));		
	srand(rank+1);
	
	for(int i=0; i<n; i++){

		for(int j=0; j<n; j++){
			if(k==0){
				(*buffA[0])[i*n + j] = rand()%100;
				(*buffB[0])[i*n + j] = rand()%100;
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

int main(int argc, char** argv) {

	
	MPI_Init(&argc, &argv);
	
	int rank;
	int ks, is, js;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
 	int c = atoi(argv[3]);
	int verbose = 0;
	if(argc > 4) verbose = atoi(argv[4]);

	int l = sqrt(p/c);
	int k = rank/(l*l);
	int i = (rank - (l*l)*k)/l;
	int j = rank%l;
	int size = n/l;
	int b = size*size;
	
	matrix A, B, C;

	MPI_Comm kcomm;
    	MPI_Comm icomm;
    	MPI_Comm jcomm;
    
    	MPI_Comm_split(MPI_COMM_WORLD, computeRank(i,j,0,l), rank, &kcomm);
	MPI_Comm_split(MPI_COMM_WORLD, computeRank(i,0,k,l), rank, &icomm);
	MPI_Comm_split(MPI_COMM_WORLD, computeRank(0,j,k,l), rank, &jcomm);

	MPI_Comm_rank(kcomm, &ks);
	MPI_Comm_rank(icomm, &is);
	MPI_Comm_rank(jcomm, &js);
	//printf("rank is js ks %d %d %d %d\n", rank, is, js, ks);
	
	int **buffA, **buffB, **buffC;
	
	initialise_matrices(&buffA, &buffB, &buffC, size, rank, k);

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



		int r = (l + j + i - (k*l)/c)%l;
		int s = (l + j - i + (k*l)/c)%l;
		
		int r1 = (l + i + j - (k*l)/c)%l;
		int s1 = (l + i - j + (k*l)/c)%l;
		


		MPI_Status status[4];
		MPI_Request reqs[4];

		MPI_Isend(buffA[0], b, MPI_INT, s, 0, icomm, &reqs[0]);
		MPI_Irecv(buffA[1], b, MPI_INT, r, 0, icomm, &reqs[1]);
	
		MPI_Isend(buffB[0], b, MPI_INT, s1, 0, jcomm, &reqs[2]);
		MPI_Irecv(buffB[1], b, MPI_INT, r1, 0, jcomm, &reqs[3]);

		MPI_Waitall(4, reqs, status);
 
		//printf("rank %d sendA %d recvA %d sendB %d recvB %d\n", rank, s, r, s1, r);	
		//fill_matrix(recvA, A);		
		//fill_matrix(recvB, B);
		
		matrix_multiply(buffA[1],buffB[1],buffC[0],size);
		
		s = (l + j + 1)%l;
		s1 = (l + i + 1)%l;
		r = (l + j - 1)%l;	
		r1 = (l + i - 1)%l;

		for(int t=1; t < l/c; t++){
			
			//copy_matrix(sendB, B);
			//copy_matrix(sendA, A);
			
	//		printf("rank %d sendA %d recvA %d sendB %d recvB %d \n", rank, s, r, s1, r1);	
			MPI_Isend(buffB[t%2], b, MPI_INT, s1, 0, jcomm, &reqs[0]);
			MPI_Irecv(buffB[1-t%2], b, MPI_INT, r1, 0, jcomm, &reqs[1]);
		

			MPI_Isend(buffA[t%2], b, MPI_INT, s, 0, icomm, &reqs[2]);
			MPI_Irecv(buffA[1-t%2], b, MPI_INT, r, 0, icomm, &reqs[3]);
			MPI_Waitall(4, reqs, status);
			
			//fill_matrix(recvA, A);
			//fill_matrix(recvB, B);
		
			matrix_multiply(buffA[1-t%2],buffB[1-t%2],buffC[0],size);
			
		}		
	

	MPI_Reduce(buffC[0], buffC[1], b, MPI_INT, MPI_SUM, 0, kcomm);
	if(k==0 and verbose){
		print_matrix(buffC[1], n/l);
		
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0) {
		double s = timer_elapsed();
		printf("time is %f \n", s);
	}
 	MPI_Finalize();
}
