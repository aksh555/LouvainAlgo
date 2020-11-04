#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <omp.h>
#include <mpi.h>
#include "dspl.hpp"
using namespace std;

static string inputFileName;
static int me, nprocs;
static int ranksPerNode = 1;
static GraphWeight threshold = 1.0E-2;

static void parseCommandLine(const int argc, char *const argv[]);

int main(int argc, char *argv[])
{
  double t0, t1, t2, t3, ti = 0.0;
  int max_threads;

  max_threads = omp_get_max_threads();

  if (max_threads > 1)
  {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
    {
      cerr << "MPI library does not support MPI_THREAD_FUNNELED." << endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
    }
  }
  else
  {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  parseCommandLine(argc, argv);
  createCommunityMPIType();
  double td0, td1, td, tdt;
  MPI_Barrier(MPI_COMM_WORLD);
  td0 = MPI_Wtime();

  Graph *g = nullptr;

  // read input graph
  BinaryEdgeList rm;
  g = rm.read(me, nprocs, ranksPerNode, inputFileName);

  assert(g != nullptr);
#ifdef PRINT_DIST_STATS
  g->print_dist_stats();
#endif

  MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINTF
  assert(g);
#endif
  td1 = MPI_Wtime();
  td = td1 - td0;

  MPI_Reduce(&td, &tdt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  GraphWeight currMod = -1.0;
  GraphWeight prevMod = -1.0;
  double total = 0.0;

  vector<GraphElem> ssizes, rsizes, svdata, rvdata;
  size_t ssz = 0, rsz = 0;
  int iters = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  t1 = MPI_Wtime();

  currMod = distLouvainMethod(me, nprocs, *g, ssz, rsz, ssizes, rsizes,svdata, rvdata, currMod, threshold, iters);

  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  total = t0 - t1;
  double tot_time = 0.0;
  MPI_Reduce(&total, &tot_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (me == 0)
  {
    double avgt = (tot_time / nprocs);
    cout<< "\n************************************* \n\n";
    cout << "\nExecution Summary\n\n";
    cout << "\tAverage total time: " << avgt << "s\n" << "\tNo. of Processes used: " << nprocs << "\n";
    cout << "\tFinal Modularity: " << currMod << "\n";
    cout << "\tTotal no. of phases: " << iters << "\n";
  }

  MPI_Barrier(MPI_COMM_WORLD);
  delete g;
  destroyCommunityMPIType();
  MPI_Finalize();
  return 0;
}

void parseCommandLine(const int argc, char *const argv[])
{
  int ret;
  while ((ret = getopt(argc, argv, "f:br:t:n:wlp:")) != -1)
  {
    switch (ret)
    {
    case 'f':
      inputFileName.assign(optarg);
      break;
    case 't':
      threshold = atof(optarg);
      break;
    default:
      assert(0 && "Should not reach here!!");
      break;
    }
  }

  if (me == 0 && (argc == 1))
  {
    cerr << "Must specify some options." << endl;
    MPI_Abort(MPI_COMM_WORLD, -99);
  }

  if (me == 0 && inputFileName.empty())
  {
    cerr << "Must specify a binary file name with -f or provide parameters for generating a graph." << endl;
    MPI_Abort(MPI_COMM_WORLD, -99);
  }
}
