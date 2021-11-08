#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
struct matrix{
	int** mat;
	int size;
};

void copy_matrix(int* buff, matrix &A){
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
}

void matrix_multiply(matrix &A, matrix &B, matrix &C){
	int n=A.size;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				C.mat[i][j] += A.mat[i][k] * B.mat[k][j];
			}
		}
	}
}

int computeRank(int i, int j, int k, int l){
	return k*l*l + i*l + j;
}

void initialise_matrices(matrix &A, matrix &B, matrix &C, int n, int rank, int k){
	A.size = n;
	B.size = n;
	C.size = n;

	A.mat = (int **) malloc(n * sizeof(int*));	
	B.mat = (int **) malloc(n * sizeof(int*));
	C.mat = (int **) malloc(n * sizeof(int*));
	

	srand(rank+1);
	
	for(int i=0; i<n; i++){
		A.mat[i] = (int *) malloc(n * sizeof(int));	
		B.mat[i] = (int *) malloc(n * sizeof(int));
		C.mat[i] = (int *) malloc(n * sizeof(int));

		for(int j=0; j<n; j++){
			if(k==0){
				A.mat[i][j] = rand()%100;
				B.mat[i][j] = rand()%100;
			}
			else{
				A.mat[i][j] = 0;
				B.mat[i][j] = 0;
			}
			C.mat[i][j] = 0;
			//printf(" i is %d j is %d a[i][j] %d b[i][j] %d\n", i,j, A.mat[i][j], B.mat[i][j]);	
		}
	}	
}

void print_matrix(matrix &A){
	int n  = A.size;
	printf("[");
	for(int i=0; i<n; i++){
		printf("[");
		for(int j=0; j<n-1; j++){
			printf("%d, ", A.mat[i][j]);
		}
		printf("%d]", A.mat[i][n-1]);
		if(i!=n-1) printf(",\n");
	}
	printf("]\n");
}

int main(int argc, char** argv) {

	
	MPI_Init(&argc, &argv);
	
	int rank;
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int n = atoi(argv[1]);
	int p = atoi(argv[2]);
 	int c = 2;
	
	int l = sqrt(p/c);
	int k = rank/(l*l);
	int i = (rank - (l*l)*k)/l;
	int j = rank%l;
	int b = (n*n)/(l*l);
	
	matrix A, B, C;

	MPI_Comm kcomm;
    	MPI_Comm icomm;
    	MPI_Comm jcomm;
    
    	MPI_Comm_split(MPI_COMM_WORLD, computeRank(i,j,0,l), rank, &kcomm);
	MPI_Comm_split(MPI_COMM_WORLD, computeRank(0,j,k,l), rank, &icomm);
	MPI_Comm_split(MPI_COMM_WORLD, computeRank(i,0,k,l), rank, &jcomm);
	
	int *sendA, *recvA, *sendB, *recvB, *sendC, *recvC;
	sendA = (int *)malloc (b * sizeof(int)); 
	recvA = (int *)malloc (b * sizeof(int)); 

	sendB = (int *)malloc (b * sizeof(int)); 
	recvB = (int *)malloc (b * sizeof(int)); 

	sendC = (int *)malloc (b * sizeof(int)); 
	recvC = (int *)malloc (b * sizeof(int)); 

	initialise_matrices(A,B,C,n/l,rank,k);
	
	if(k==0){
		copy_matrix(sendA, A);
		copy_matrix(sendB, B);
		print_matrix(A);
		print_matrix(B);
	}
	printf("rank: %d computed %d \n", rank, computeRank(i,j,0,l));
	MPI_Bcast(sendA, b, MPI_INT, computeRank(i,j,0,l), kcomm);
	MPI_Bcast(sendB, b, MPI_INT, computeRank(i,j,0,l), kcomm);

	if(k==0){
		MPI_Reduce(sendC, recvC, b, MPI_INT, MPI_SUM, rank, kcomm);
	}
	else{
		fill_matrix(sendA, A);
		fill_matrix(sendB, B);
		int r = (l + j + i - (k*l)/c)%l;
		int s = (l + j - i + (k*l)/c)%l;
		print_matrix(A);
		print_matrix(B);
		
		MPI_Status statusA;
		MPI_Status statusB;

		copy_matrix(sendA, A);

		MPI_Send(sendA, b, MPI_INT, computeRank(i,s,k,l), 0, icomm);
		MPI_Recv(recvA, b, MPI_INT, computeRank(i,r,k,l), 0, icomm, &statusA);

		fill_matrix(recvA, A);
		
		int r1 = (l + i + j - (k*l)/c)%l;
		int s1 = (l + i - j + (k*l)/c)%l;
		
		copy_matrix(sendB, B);
		
		MPI_Send(sendB, b, MPI_INT, computeRank(s1,j,k,l), 0, jcomm);
		MPI_Recv(recvB, b, MPI_INT, computeRank(r1,j,k,l), 0, jcomm, &statusB);

		fill_matrix(recvB, B);
		
		matrix_multiply(A,B,C);
		
		s = (l + j + 1)%l;
		s1 = (l + i + 1)%l;

		for(int t=1; t < l/c; t++){
			
			copy_matrix(sendB, B);
			copy_matrix(sendA, A);
			
			MPI_Send(sendA, b, MPI_INT, computeRank(i,s,k,l), 0, icomm);
			MPI_Send(sendB, b, MPI_INT, computeRank(s1,j,k,l), 0, jcomm);

			r = (l + r - 1)%l;

			MPI_Recv(recvA, b, MPI_INT, computeRank(i,r,k,l), 0, icomm, &statusA);
			MPI_Recv(recvB, b, MPI_INT, computeRank(r,j,k,l), 0, jcomm, &statusB);

			matrix_multiply(A,B,C);
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
 	MPI_Finalize();
}
