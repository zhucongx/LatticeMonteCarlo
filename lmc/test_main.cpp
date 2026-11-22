#include <gtest/gtest.h>
#include <mpi.h>

// The simulation classes ask MPI for their rank, so MPI has to be running while the tests execute.
int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  const int result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
