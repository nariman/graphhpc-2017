/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <cstring>  // std::memset
#include <iostream>  // std::cout

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
         double*        result) {
    unsigned threads = omp_get_max_threads();

    vertex_size_t* order   = new vertex_size_t[threads * vertices];
    edge_size_t*   p_last  = new   edge_size_t[threads * vertices];
    edge_size_t*   p_prev  = new   edge_size_t[threads * edges + 1];
    vertex_size_t* p_val   = new vertex_size_t[threads * edges + 1];
    edge_size_t*   sigma   = new   edge_size_t[threads * vertices];
    vertex_size_t* depth   = new vertex_size_t[threads * vertices];
    bool*          visited = new          bool[threads * vertices];
    double*        delta   = new        double[threads * vertices];
    vertex_size_t* queue   = new vertex_size_t[threads * vertices];

    p_prev++;
    p_val++;

    #pragma omp parallel firstprivate(order, p_last, p_prev, p_val, sigma, depth, visited, delta, queue)
    {
        unsigned tid = omp_get_thread_num();

        order   += vertices * tid;
        p_last  += vertices * tid;
        p_prev  += edges    * tid;
        p_val   += edges    * tid;
        sigma   += vertices * tid;
        depth   += vertices * tid;
        visited += vertices * tid;
        delta   += vertices * tid;
        queue   += vertices * tid;

        #pragma omp for
        for (vertex_size_t s = 0; s < vertices; s++) {
            memset(p_last, 0, sizeof(edge_size_t) * vertices);
            memset(sigma,  0, sizeof(edge_size_t) * vertices);

            fill(visited, visited + vertices, false);
            fill(delta,   delta   + vertices, 0.0);

            vertex_size_t order_pos   = 0;
            vertex_size_t p_pos       = 0;
            vertex_size_t queue_front = 0;
            vertex_size_t queue_back  = 0;

            sigma[s]      = 1;
            depth[s]      = 0;
            visited[s]    = true;
            queue[queue_back++] = s;   

            while (queue_front != queue_back) {
                vertex_size_t v = queue[queue_front++];
                order[order_pos++] = v;

                for (edge_size_t i = indices[v]; i < indices[v + 1]; i++) {
                    vertex_size_t t = ends[i];

                    if (!visited[t]) {
                        depth[t] = depth[v] + 1;
                        queue[queue_back++] = t;
                        visited[t] = true;
                    }

                    if (depth[t] > depth[v] && depth[t] == depth[v] + 1) {
                        sigma[t] += sigma[v];
                        
                        p_prev[p_pos] = p_last[t];
                        p_val[p_pos] = v;
                        p_last[t] = p_pos;
                        p_pos++;
                    }
                }
            }

            while (order_pos --> 1)  {
                vertex_size_t v = order[order_pos];
                edge_size_t i = p_last[v];
                double d = (1 + delta[v]) / (double) sigma[v];

                while (i != 0) {
                    vertex_size_t t = p_val[i];
                    delta[t] += sigma[t] * d;
                    i = p_prev[i];
                }

                #pragma omp atomic
                result[v] += delta[v] / 2;
            }
        }
    }

    delete[] queue;
    delete[] delta;
    delete[] visited;
    delete[] depth;
    delete[] sigma;
    delete[] --p_val;
    delete[] --p_prev;
    delete[] p_last;
    delete[] order;
}