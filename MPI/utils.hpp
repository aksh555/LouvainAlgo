#pragma once
#ifndef UTILS_HPP
#define UTILS_HPP
#define PI (3.14159)
#define MAX_PRINT_NEDGE (100000)
#define MLCG (2147483647) // 2^31 - 1
#define ALCG (16807)      // 7^5
#define BLCG (0)
#define SR_UP_TAG 100
#define SR_DOWN_TAG 101
#define SR_SIZES_UP_TAG 102
#define SR_SIZES_DOWN_TAG 103
#define SR_X_UP_TAG 104
#define SR_X_DOWN_TAG 105
#define SR_Y_UP_TAG 106
#define SR_Y_DOWN_TAG 107
#define SR_LCG_TAG 108

#include <random>
#include <utility>
#include <cstring>
using namespace std;

using GraphElem = int64_t;
using GraphWeight = double;
const MPI_Datatype MPI_GRAPH_TYPE = MPI_INT64_T;
const MPI_Datatype MPI_WEIGHT_TYPE = MPI_DOUBLE;

extern unsigned seed;

// Is nprocs a power-of-2?
int is_pwr2(int nprocs)
{
    return ((nprocs != 0) && !(nprocs & (nprocs - 1)));
}

// return unint32_t seed
GraphElem reseeder(unsigned initseed)
{
    seed_seq seq({initseed});
    vector<uint32_t> seeds(1);
    seq.generate(seeds.begin(), seeds.end());

    return (GraphElem)seeds[0];
}

// Local random number generator
template <typename T, typename G = default_random_engine>
T genRandom(T lo, T hi)
{
    thread_local static G gen(seed);
    using Dist = typename conditional<
        is_integral<T>::value, uniform_int_distribution<T>, uniform_real_distribution<T>>::type;

    thread_local static Dist utd{};
    return utd(gen, typename Dist::param_type{lo, hi});
}

// Parallel Linear Congruential Generator
// x[i] = (a*x[i-1] + b)%M
class LCG
{
public:
    LCG(unsigned seed, GraphWeight *drand,
        GraphElem n, MPI_Comm comm = MPI_COMM_WORLD) : seed_(seed), drand_(drand), n_(n)
    {
        comm_ = comm;
        MPI_Comm_size(comm_, &nprocs_);
        MPI_Comm_rank(comm_, &rank_);

        // allocate long random numbers
        rnums_.resize(n_);

        // init x0
        if (rank_ == 0)
            x0_ = reseeder(seed_);

        // step #1: bcast x0 from root
        MPI_Bcast(&x0_, 1, MPI_GRAPH_TYPE, 0, comm_);

        // step #2: parallel prefix to generate first random value per process
        parallel_prefix_op();
    }

    ~LCG() { rnums_.clear(); }

    // matrix-matrix multiplication for 2x2 matrices
    void matmat_2x2(GraphElem c[], GraphElem a[], GraphElem b[])
    {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                GraphElem sum = 0;
                for (int k = 0; k < 2; k++)
                {
                    sum += a[i * 2 + k] * b[k * 2 + j];
                }
                c[i * 2 + j] = sum;
            }
        }
    }

    // x *= y
    void matop_2x2(GraphElem x[], GraphElem y[])
    {
        GraphElem tmp[4];
        matmat_2x2(tmp, x, y);
        memcpy(x, tmp, sizeof(GraphElem[4]));
    }

    // find kth power of a 2x2 matrix
    void mat_power(GraphElem mat[], GraphElem k)
    {
        GraphElem tmp[4];
        memcpy(tmp, mat, sizeof(GraphElem[4]));

        // mat-mat multiply k times
        for (GraphElem p = 0; p < k - 1; p++)
            matop_2x2(mat, tmp);
    }

    // parallel prefix for matrix-matrix operation
    // `x0 is the very first random number in the series
    // `ab is a 2-length array which stores a and b
    // `n_ is (n/p)
    // `rnums is n_ length array which stores the random nums for a process
    void parallel_prefix_op()
    {
        GraphElem global_op[4];
        global_op[0] = ALCG;
        global_op[1] = 0;
        global_op[2] = BLCG;
        global_op[3] = 1;

        mat_power(global_op, n_);              // M^(n/p)
        GraphElem prefix_op[4] = {1, 0, 0, 1}; // I in row-major

        GraphElem global_op_recv[4];

        int steps = (int)(log2((double)nprocs_));

        for (int s = 0; s < steps; s++)
        {

            int mate = rank_ ^ (1 << s); // toggle the sth LSB to find my neighbor

            // send/recv global to/from mate
            MPI_Sendrecv(global_op, 4, MPI_GRAPH_TYPE, mate, SR_LCG_TAG,
                         global_op_recv, 4, MPI_GRAPH_TYPE, mate, SR_LCG_TAG,
                         comm_, MPI_STATUS_IGNORE);

            matop_2x2(global_op, global_op_recv);

            if (mate < rank_)
                matop_2x2(prefix_op, global_op_recv);

            MPI_Barrier(comm_);
        }

        // populate the first random number entry for each process
        // (x0*a + b)%P
        if (rank_ == 0)
            rnums_[0] = x0_;
        else
            rnums_[0] = (x0_ * prefix_op[0] + prefix_op[2]) % MLCG;
    }

    // generate random number based on the first
    // random number on a process
    // TODO check the 'quick'n dirty generators to
    // see if we can avoid the mod
    void generate()
    {

        for (GraphElem i = 1; i < n_; i++)
        {
            rnums_[i] = (rnums_[i - 1] * ALCG + BLCG) % MLCG;
        }

        GraphWeight mult = 1.0 / (GraphWeight)(1.0 + (GraphWeight)(MLCG - 1));

        for (GraphElem i = 0; i < n_; i++)
            drand_[i] = (GraphWeight)((GraphWeight)fabs(rnums_[i]) * mult); // 0-1
    }

    // copy from drand_[idx_start] to new_drand,

    void rescale(GraphWeight *new_drand, GraphElem idx_start, GraphWeight const &lo)
    {
        GraphWeight range = (1.0 / (GraphWeight)nprocs_);

        for (GraphElem i = idx_start, j = 0; i < n_; i++, j++)
            new_drand[j] = lo + (GraphWeight)(range * drand_[i]); // lo-hi
    }

private:
    MPI_Comm comm_;
    int nprocs_, rank_;
    unsigned seed_;
    GraphElem n_, x0_;
    GraphWeight *drand_;
    vector<GraphElem> rnums_;
};

// locks
#ifdef USE_OPENMP_LOCK
#else
#ifdef USE_SPINLOCK
#include <atomic>
atomic_flag lkd_ = ATOMIC_FLAG_INIT;
#else
#include <mutex>
mutex mtx_;
#endif
void lock()
{
#ifdef USE_SPINLOCK
    while (lkd_.test_and_set(memory_order_acquire))
    {
        ;
    }
#else
    mtx_.lock();
#endif
}
void unlock()
{
    mtx_.unlock();
#endif
}
#endif // UTILS
