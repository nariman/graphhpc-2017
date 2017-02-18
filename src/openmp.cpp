/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <cstring>  // std::memset
#include <iostream> // std::cout

#include <omp.h>

#include "main.h"

using namespace std;


/**
 * Calculates betweenness centrality by the algorithm with OpenMP stack.
 */
void run(vertex_size_t  vertices, 
         edge_size_t    edges, 
         vertex_size_t* ends, 
         edge_size_t*   indices, 
         double*        result) { // don't touch the result array!

    // Constants

    unsigned threads = omp_get_max_threads();


    // Allocations

    vertex_size_t* order   = new vertex_size_t[threads * vertices];
    edge_size_t*   p_last  = new   edge_size_t[threads * vertices];
    edge_size_t*   p_prev  = new   edge_size_t[threads * (edges + 1)];
    vertex_size_t* p_val   = new vertex_size_t[threads * (edges + 1)];
    edge_size_t*   sigma   = new   edge_size_t[threads * vertices];
    vertex_size_t* depth   = new vertex_size_t[threads * vertices];
    double*        delta   = new        double[threads * vertices];
    vertex_size_t* queue   = new vertex_size_t[threads * vertices];

    double*        res     = new        double[threads * vertices];

    memset(res, 0, sizeof(double) * threads * vertices);


    // Main logic

    #pragma omp parallel firstprivate(order, p_last, p_prev, p_val, sigma, depth, delta, queue, res)
    {
        unsigned tid = omp_get_thread_num();

        order   += vertices * tid;
        p_last  += vertices * tid;
        p_prev  += edges    * tid + 1;
        p_val   += edges    * tid + 1;
        sigma   += vertices * tid;
        depth   += vertices * tid;
        delta   += vertices * tid;
        queue   += vertices * tid;

        res     += vertices * tid;

        vertex_size_t* order_pos;
        edge_size_t    p_pos;
        vertex_size_t* queue_front;
        vertex_size_t* queue_back;

        #pragma omp for
        for (vertex_size_t s = 0; s < vertices; ++s) {
            memset(p_last, 0, sizeof(edge_size_t)   * vertices);
            memset(sigma,  0, sizeof(edge_size_t)   * vertices);
            memset(depth,  0, sizeof(vertex_size_t) * vertices);
            memset(delta,  0, sizeof(double)        * vertices);

            order_pos    = order;
            p_pos        = 1;
            queue_front  = queue;
            queue_back   = queue;

            sigma[s]      = 1;
            depth[s]      = 1;
            *queue_back++ = s;

            while (queue_front != queue_back) {
                vertex_size_t v = *queue_front++;
                *order_pos++ = v;

                for (vertex_size_t *t = &ends[indices[v]], 
                                   *r = &ends[indices[v + 1]];
                     t != r; 
                     ++t) {
                    if (depth[*t] == 0) {
                        depth[*t] = depth[v] + 1;
                        *queue_back++ = *t;
                    }

                    if (depth[*t] > depth[v]) {
                        sigma[*t] += sigma[v];

                        p_prev[p_pos] = p_last[*t];
                        p_val[p_pos] = v;
                        p_last[*t] = p_pos++;
                    }
                }
            }

            while (--order_pos != order)  {
                edge_size_t i = p_last[*order_pos];
                double d = (1 + delta[*order_pos]) / (double) sigma[*order_pos];

                while (i != 0) {
                    delta[p_val[i]] += sigma[p_val[i]] * d;
                    i = p_prev[i];
                }

                res[*order_pos] += delta[*order_pos] / 2;
            }
        }
    }


    // Result calculation

    #pragma omp parallel firstprivate(res)
    {
        unsigned tid = omp_get_thread_num();
        res += vertices * tid;

        for (vertex_size_t v = 0; v < vertices; v++) {
            #pragma omp atomic
            result[v] += res[v];
        }
    }


    // Deal with it!

    delete[] res;

    delete[] queue;
    delete[] delta;
    delete[] depth;
    delete[] sigma;
    delete[] p_val;
    delete[] p_prev;
    delete[] p_last;
    delete[] order;
}
