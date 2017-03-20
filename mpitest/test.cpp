#include <mpi.h>
#include <cstdio>

int main(int argc, char** argv)
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // get name of processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // output message
    std::cout << world_rank << ": processor " << processor_name << ", total: " << world_size << " processors" << std::endl;

    // Finalize the MPI environment.
    MPI_Finalize();
}
