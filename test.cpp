#include <math.h>
#include <mpi.h>

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
				C[i][j] += A.mat[i][k] * B.mat[k][j];
			}
		}
	}
}

int computeRank(int i, int j, int k, int l){
	return k*l*l + i*l + j;
}

int main(int argc, char** argv) {

	
	MPI_Init(NULL, NULL);
	
	int n = atoi(argv[0]);
	int p = atoi(argv[1]);
 	int c;
	
	int l = sqrt(p/c);
	int k = rank/(l*l);
	int i = (rank - (l*l)*k)/l;
	int j = rank%l;
	int b = (n*n)/(l*l);

	int rank;
	MPI_Comm kcomm;
    MPI_Comm icomm;
    MPI_Comm jcomm;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, k, rank, &kcomm);
	MPI_Comm_split(MPI_COMM_WORLD, i, rank, &icomm);
	MPI_Comm_split(MPI_COMM_WORLD, j, rank, &jcomm);
	
	int *sendA, *recvA, *sendB, *recvB, *sendC, *recvC;
	sendA = malloc (b * sizeof(int)); 
	recvA = malloc (b * sizeof(int)); 

	sendB = malloc (b * sizeof(int)); 
	recvB = malloc (b * sizeof(int)); 

	sendC = malloc (b * sizeof(int)); 
	recvC = malloc (b * sizeof(int)); 

	if(k==0){
		copy_matrix(sendA, A);
		copy_matrix(sendB, B);
		copy_matrix(sendC, C);

		MPI_Bcast(sendA, b, MPI_INT, rank, kcomm);

		MPI_Reduce(sendC, recvC, b, MPI_INT, MPI_SUM, rank, kcomm);
	}
	else{
		MPI_RECV(A,B);
		
		int r = (l + j + i − (k*l)/c)%l;
		int s = (l + j - i + (k*l)/c)%l;
		
		MPI_Status statusA;
		MPI_Status statusB;

		copy_matrix(sendA, A);

		MPI_Send(sendA, b, MPI_INT, computeRank(i,s,k), 0, icomm);
		MPI_Recv(recvA, b, MPI_INT, computeRank(i,r,k), 0, icomm, statusA);

		fill_matrix(recvA, A);
		
		int r1 = (l + i + j - (k*l)/c)%l;
		int s1 = (l + i - j + (k*l)/c)%l;
		

		copy_matrix(sendB, B);
		
		MPI_Send(sendB, b, MPI_INT, computeRank(s1,j,k), 0, jcomm);
		MPI_Recv(recvB, b, MPI_INT, computeRank(r1,j,k), 0, jcomm, statusB);

		fill_matrix(recvB, B);
		
		matrix_multiply(A,B,C);
		
		s = (l + j + 1)%l;
		s1 = (l + i + 1)%l;

		for(int t=1; t < l/c; t++){
			
			copy_matrix(sendB, B);
			copy_matrix(sendA, A);
			
			MPI_Send(sendA, b, MPI_INT, computeRank(i,s,k), 0, icomm);
			MPI_Send(sendB, b, MPI_INT, computeRank(s1,j,k), 0, jcomm);

			r = (l + r - 1)%l;

			MPI_Recv(recvA, b, MPI_INT, computeRank(i,r,k), 0, icomm, statusA);
			MPI_Recv(recvB, b, MPI_INT, computeRank(r,j,k), 0, jcomm, statusB);

			matrix_multiply(A,B,C);
		}		
	}

	MPI_Barrier(MPI_COMM_WORLD);
 	MPI_Finalize();
}