#pragma once
#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <bits/stdc++.h>
#include <mpi.h>
#include "utils.hpp"
using namespace std;

unsigned seed;

struct Edge
{
    GraphElem tail_;
    GraphWeight weight_;

    Edge() : tail_(-1), weight_(0.0) {}
};

struct EdgeTuple
{
    GraphElem ij_[2];
    GraphWeight w_;

    EdgeTuple(GraphElem i, GraphElem j, GraphWeight w) : ij_{i, j}, w_(w)
    {
    }
    EdgeTuple(GraphElem i, GraphElem j) : ij_{i, j}, w_(1.0)
    {
    }
    EdgeTuple() : ij_{-1, -1}, w_(0.0)
    {
    }
};

// per process graph instance
class Graph
{
public:
    Graph() : lnv_(-1), lne_(-1), nv_(-1),
              ne_(-1), comm_(MPI_COMM_WORLD)
    {
        MPI_Comm_size(comm_, &size_);
        MPI_Comm_rank(comm_, &rank_);
    }

    Graph(GraphElem lnv, GraphElem lne,
          GraphElem nv, GraphElem ne,
          MPI_Comm comm = MPI_COMM_WORLD) : lnv_(lnv), lne_(lne),
                                            nv_(nv), ne_(ne),
                                            comm_(comm)
    {
        MPI_Comm_size(comm_, &size_);
        MPI_Comm_rank(comm_, &rank_);

        edge_indices_.resize(lnv_ + 1, 0);
        edge_list_.resize(lne_); // this is usually populated later

        parts_.resize(size_ + 1);
        parts_[0] = 0;

        for (GraphElem i = 1; i < size_ + 1; i++)
            parts_[i] = ((nv_ * i) / size_);
    }

    ~Graph()
    {
        edge_list_.clear();
        edge_indices_.clear();
        parts_.clear();
    }

    // update vertex partition information
    void repart(vector<GraphElem> const &parts)
    {
        memcpy(parts_.data(), parts.data(), sizeof(GraphElem) * (size_ + 1));
    }

    void set_edge_index(GraphElem const vertex, GraphElem const e0)
    {
#if defined(DEBUG_BUILD)
        assert((vertex >= 0) && (vertex <= lnv_));
        assert((e0 >= 0) && (e0 <= lne_));
        edge_indices_.at(vertex) = e0;
#else
        edge_indices_[vertex] = e0;
#endif
    }

    void edge_range(GraphElem const vertex, GraphElem &e0,
                    GraphElem &e1) const
    {
        e0 = edge_indices_[vertex];
        e1 = edge_indices_[vertex + 1];
    }

    // collective
    void set_nedges(GraphElem lne)
    {
        lne_ = lne;
        edge_list_.resize(lne_);

        // compute total number of edges
        ne_ = 0;
        MPI_Allreduce(&lne_, &ne_, 1, MPI_GRAPH_TYPE, MPI_SUM, comm_);
    }

    GraphElem get_base(const int rank) const
    {
        return parts_[rank];
    }

    GraphElem get_bound(const int rank) const
    {
        return parts_[rank + 1];
    }

    GraphElem get_range(const int rank) const
    {
        return (parts_[rank + 1] - parts_[rank] + 1);
    }

    int get_owner(const GraphElem vertex) const
    {
        const vector<GraphElem>::const_iterator iter =
            upper_bound(parts_.begin(), parts_.end(), vertex);

        return (iter - parts_.begin() - 1);
    }

    GraphElem get_lnv() const { return lnv_; }
    GraphElem get_lne() const { return lne_; }
    GraphElem get_nv() const { return nv_; }
    GraphElem get_ne() const { return ne_; }
    MPI_Comm get_comm() const { return comm_; }

    Edge const &get_edge(GraphElem const index) const
    {
        return edge_list_[index];
    }

    Edge &set_edge(GraphElem const index)
    {
        return edge_list_[index];
    }

    // local <--> global index translation
    GraphElem local_to_global(GraphElem idx)
    {
        return (idx + get_base(rank_));
    }

    GraphElem global_to_local(GraphElem idx)
    {
        return (idx - get_base(rank_));
    }

    GraphElem local_to_global(GraphElem idx, int rank)
    {
        return (idx + get_base(rank));
    }

    GraphElem global_to_local(GraphElem idx, int rank)
    {
        return (idx - get_base(rank));
    }

    void print(bool print_weight = true) const
    {
        if (lne_ < MAX_PRINT_NEDGE)
        {
            for (int p = 0; p < size_; p++)
            {
                MPI_Barrier(comm_);
                if (p == rank_)
                {
                    cout << "###############" << endl;
                    cout << "Process #" << p << ": " << endl;
                    cout << "###############" << endl;
                    GraphElem base = get_base(p);
                    for (GraphElem i = 0; i < lnv_; i++)
                    {
                        GraphElem e0, e1;
                        edge_range(i, e0, e1);
                        if (print_weight)
                        {
                            for (GraphElem e = e0; e < e1; e++)
                            {
                                Edge const &edge = get_edge(e);
                                cout << i + base << " " << edge.tail_ << " " << edge.weight_ << endl;
                            }
                        }
                        else
                        {
                            for (GraphElem e = e0; e < e1; e++)
                            {
                                Edge const &edge = get_edge(e);
                                cout << i + base << " " << edge.tail_ << endl;
                            }
                        }
                    }
                    MPI_Barrier(comm_);
                }
            }
        }
        else
        {
            if (rank_ == 0)
                cout << "Graph size per process is {" << lnv_ << ", " << lne_ << "}, which will overwhelm STDOUT." << endl;
        }
    }

    // print statistics about edge distribution
    void print_dist_stats()
    {
        long sumdeg = 0, maxdeg = 0;
        long lne = (long)lne_;

        MPI_Reduce(&lne, &sumdeg, 1, MPI_LONG, MPI_SUM, 0, comm_);
        MPI_Reduce(&lne, &maxdeg, 1, MPI_LONG, MPI_MAX, 0, comm_);

        long my_sq = lne * lne;
        long sum_sq = 0;
        MPI_Reduce(&my_sq, &sum_sq, 1, MPI_LONG, MPI_SUM, 0, comm_);

        double average = (double)sumdeg / size_;
        double avg_sq = (double)sum_sq / size_;
        double var = avg_sq - (average * average);
        double stddev = sqrt(var);

        MPI_Barrier(comm_);

        if (rank_ == 0)
        {
            cout << endl;
            cout << "-------------------------------------------------------" << endl;
            cout << "Graph edge distribution characteristics" << endl;
            cout << "-------------------------------------------------------" << endl;
            cout << "Number of vertices: " << nv_ << endl;
            cout << "Number of edges: " << ne_ << endl;
            cout << "Maximum number of edges: " << maxdeg << endl;
            cout << "Average number of edges: " << average << endl;
            cout << "Expected value of X^2: " << avg_sq << endl;
            cout << "Variance: " << var << endl;
            cout << "Standard deviation: " << stddev << endl;
            cout << "-------------------------------------------------------" << endl;
        }
    }
    vector<GraphElem> edge_indices_;
    vector<Edge> edge_list_;

private:
    GraphElem lnv_, lne_, nv_, ne_;
    vector<GraphElem> parts_;
    MPI_Comm comm_;
    int rank_, size_;
};

class BinaryEdgeList
{
public:
    BinaryEdgeList() : M_(-1), N_(-1),
                       M_local_(-1), N_local_(-1),
                       comm_(MPI_COMM_WORLD)
    {
    }
    BinaryEdgeList(MPI_Comm comm) : M_(-1), N_(-1),
                                    M_local_(-1), N_local_(-1),
                                    comm_(comm)
    {
    }

    Graph *read(int me, int nprocs, int ranks_per_node, string file)
    {
        int file_open_error;
        MPI_File fh;
        MPI_Status status;

        MPI_Info info;
        MPI_Info_create(&info);
        int naggr = (ranks_per_node > 1) ? (nprocs / ranks_per_node) : ranks_per_node;
        if (naggr >= nprocs)
            naggr = 1;
        stringstream tmp_str;
        tmp_str << naggr;
        string str = tmp_str.str();
        MPI_Info_set(info, "cb_nodes", str.c_str());

        file_open_error = MPI_File_open(comm_, file.c_str(), MPI_MODE_RDONLY, info, &fh);
        MPI_Info_free(&info);

        if (file_open_error != MPI_SUCCESS)
        {
            cout << " Error opening file! " << endl;
            MPI_Abort(comm_, -99);
        }
        MPI_File_read_all(fh, &M_, sizeof(GraphElem), MPI_BYTE, &status);
        MPI_File_read_all(fh, &N_, sizeof(GraphElem), MPI_BYTE, &status);
        M_local_ = ((M_ * (me + 1)) / nprocs) - ((M_ * me) / nprocs);

        // create local graph
        Graph *g = new Graph(M_local_, 0, M_, N_);

        // Let N = array length and P = number of processors.
        // From j = 0 to P-1,
        // Starting point of array on processor j = floor(N * j / P)
        // Length of array on processor j = floor(N * (j + 1) / P) - floor(N * j / P)

        uint64_t tot_bytes = (M_local_ + 1) * sizeof(GraphElem);
        MPI_Offset offset = 2 * sizeof(GraphElem) + ((M_ * me) / nprocs) * sizeof(GraphElem);

        if (tot_bytes < INT_MAX)
            MPI_File_read_at(fh, offset, &g->edge_indices_[0], tot_bytes, MPI_BYTE, &status);
        else
        {
            int chunk_bytes = INT_MAX;
            uint8_t *curr_pointer = (uint8_t *)&g->edge_indices_[0];
            uint64_t transf_bytes = 0;

            while (transf_bytes < tot_bytes)
            {
                MPI_File_read_at(fh, offset, curr_pointer, chunk_bytes, MPI_BYTE, &status);
                transf_bytes += chunk_bytes;
                offset += chunk_bytes;
                curr_pointer += chunk_bytes;

                if ((tot_bytes - transf_bytes) < INT_MAX)
                    chunk_bytes = tot_bytes - transf_bytes;
            }
        }

        N_local_ = g->edge_indices_[M_local_] - g->edge_indices_[0];
        g->set_nedges(N_local_);

        tot_bytes = N_local_ * (sizeof(Edge));
        offset = 2 * sizeof(GraphElem) + (M_ + 1) * sizeof(GraphElem) + g->edge_indices_[0] * (sizeof(Edge));

        if (tot_bytes < INT_MAX)
            MPI_File_read_at(fh, offset, &g->edge_list_[0], tot_bytes, MPI_BYTE, &status);
        else
        {
            int chunk_bytes = INT_MAX;
            uint8_t *curr_pointer = (uint8_t *)&g->edge_list_[0];
            uint64_t transf_bytes = 0;

            while (transf_bytes < tot_bytes)
            {
                MPI_File_read_at(fh, offset, curr_pointer, chunk_bytes, MPI_BYTE, &status);
                transf_bytes += chunk_bytes;
                offset += chunk_bytes;
                curr_pointer += chunk_bytes;

                if ((tot_bytes - transf_bytes) < INT_MAX)
                    chunk_bytes = (tot_bytes - transf_bytes);
            }
        }

        MPI_File_close(&fh);

        for (GraphElem i = 1; i < M_local_ + 1; i++)
            g->edge_indices_[i] -= g->edge_indices_[0];
        g->edge_indices_[0] = 0;

        return g;
    }

private:
    GraphElem M_;
    GraphElem N_;
    GraphElem M_local_;
    GraphElem N_local_;
    MPI_Comm comm_;
};

#endif