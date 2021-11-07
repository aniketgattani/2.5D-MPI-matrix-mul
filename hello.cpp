//
// An MPI hello world program
//

#include <stdio.h>
#include <mpi.h>

main(int argc, char **argv)
{
   int my_rank, num_ranks;
   
   MPI_Init(&argc,&argv);

   MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
     
   printf("Hello World from rank %d of %d\n", my_rank, num_ranks);
            
   MPI_Finalize();
}
