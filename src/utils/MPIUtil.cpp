#include "MPIUtil.h"
#include "mpi.h"

bool cMPIUtil::Init()
{
    return MPI_SUCCESS == MPI_Init(NULL, NULL);
}
int cMPIUtil::GetWorldSize()
{
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    return world_size;
}
int cMPIUtil::GetWorldRank()
{
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return world_rank;
}
void cMPIUtil::End()
{
    MPI_Finalize();
}