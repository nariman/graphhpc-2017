/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <iostream>  // std::cout

#include <cuda.h>
#include <thrust/fill.h>

#include "main.h"

#ifndef THREADS_PER_ITERATION
#define THREADS_PER_ITERATION 64
#endif

#ifndef BLOCKS_PER_ITERATION
#define BLOCKS_PER_ITERATION 2
#endif

#define THREADS_PER_BLOCK THREADS_PER_ITERATION / BLOCKS_PER_ITERATION

using namespace std;


__global__ void kernel(vertex_size_t* d_vertices,
                       edge_size_t*   d_edges,
                       vertex_size_t* d_ends,
                       edge_size_t*   d_indices, 
                       double*        d_result,
                       vertex_size_t* d_order,
                       edge_size_t*   d_p_last,
                       edge_size_t*   d_p_prev,
                       vertex_size_t* d_p_val,
                       edge_size_t*   d_sigma,
                       vertex_size_t* d_depth,
                       bool*          d_visited,
                       double*        d_delta,
                       vertex_size_t* d_queue) {
    int ID = threadIdx.x + blockIdx.x * blockDim.x;

    vertex_size_t vertices = *d_vertices;
    edge_size_t   edges    = *d_edges;

    if (ID >= vertices) {
        return;
    }

    d_result  += vertices * ID;

    d_order   += vertices * ID;
    d_p_last  += vertices * ID;
    d_p_prev  += edges    * ID + 1;
    d_p_val   += edges    * ID + 1;
    d_sigma   += vertices * ID;
    d_depth   += vertices * ID;
    d_visited += vertices * ID;
    d_delta   += vertices * ID;
    d_queue   += vertices * ID;

    for (vertex_size_t s = ID; s < vertices; s += THREADS_PER_ITERATION) {
        memset(d_p_last, 0, sizeof(edge_size_t) * vertices);
        memset(d_sigma,  0, sizeof(edge_size_t) * vertices);

        thrust::fill(d_visited, d_visited + vertices, false);
        thrust::fill(d_delta,   d_delta   + vertices, 0.0);

        vertex_size_t order_pos   = 0;
        vertex_size_t p_pos       = 0;
        vertex_size_t queue_front = 0;
        vertex_size_t queue_back  = 0;

        d_sigma[s]            = 1;
        d_depth[s]            = 0;
        d_visited[s]          = true;
        d_queue[queue_back++] = s;

        while (queue_front != queue_back) {
            vertex_size_t v = d_queue[queue_front++];
            d_order[order_pos++] = v;

            for (edge_size_t i = d_indices[v]; i < d_indices[v + 1]; i++) {
                vertex_size_t t = d_ends[i];

                if (!d_visited[t]) {
                    d_depth[t] = d_depth[v] + 1;
                    d_queue[queue_back++] = t;
                    d_visited[t] = true;
                }

                if (d_depth[t] > d_depth[v] && d_depth[t] == d_depth[v] + 1) {
                    d_sigma[t] += d_sigma[v];

                    d_p_prev[p_pos] = d_p_last[t];
                    d_p_val[p_pos] = v;
                    d_p_last[t] = p_pos;
                    p_pos++;
                }
            }
        }

        while (order_pos --> 1)  {
            vertex_size_t v = d_order[order_pos];
            edge_size_t i = d_p_last[v];
            double d = (1 + d_delta[v]) / (double) d_sigma[v];

            while (i != 0) {
                vertex_size_t t = d_p_val[i];
                d_delta[t] += d_sigma[t] * d;
                i = d_p_prev[i];
            }

            d_result[v] += d_delta[v] / 2;
        }
    }
}

/**
 * Calculates betweenness centrality by the algorithm with CUDA stack.
 */
void run(vertex_size_t  h_vertices, 
         edge_size_t    h_edges, 
         vertex_size_t* h_ends, 
         edge_size_t*   h_indices, 
         double*        result) {

    if (THREADS_PER_ITERATION % BLOCKS_PER_ITERATION != 0) {
        // you're idiot;
        exit(1);
    }

    // Allocatins for the graph data

    vertex_size_t* d_vertices;
    edge_size_t*   d_edges;
    vertex_size_t* d_ends;
    edge_size_t*   d_indices;

    cudaMalloc((void **) &d_vertices, sizeof(vertex_size_t));
    cudaMalloc((void **) &d_edges,    sizeof(edge_size_t));
    cudaMalloc((void **) &d_ends,     sizeof(vertex_size_t) * h_edges);
    cudaMalloc((void **) &d_indices,  sizeof(edge_size_t) * (h_vertices + 1));

    cudaMemcpy(d_vertices, &h_vertices, sizeof(vertex_size_t),                  cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges,    &h_edges,    sizeof(edge_size_t),                    cudaMemcpyHostToDevice);
    cudaMemcpy(d_ends,      h_ends,     sizeof(vertex_size_t) * h_edges,        cudaMemcpyHostToDevice);
    cudaMemcpy(d_indices,   h_indices,  sizeof(edge_size_t) * (h_vertices + 1), cudaMemcpyHostToDevice);


    // Allocating an array for results

    double* h_result = (double*) malloc(sizeof(double) * 
                                        h_vertices * 
                                        THREADS_PER_ITERATION);
    double* d_result;

    cudaMalloc((void **) &d_result, sizeof(double) * 
                                    h_vertices * 
                                    THREADS_PER_ITERATION);
    cudaMemset(d_result, 0, sizeof(double) * 
                            h_vertices * 
                            THREADS_PER_ITERATION);


    // Allocating an temporary arrays for calculations

    vertex_size_t* d_order;
    edge_size_t*   d_p_last;
    edge_size_t*   d_p_prev;
    vertex_size_t* d_p_val;
    edge_size_t*   d_sigma;
    vertex_size_t* d_depth;
    bool*          d_visited;
    double*        d_delta;
    vertex_size_t* d_queue;

    cudaMalloc((void **) &d_order,   sizeof(vertex_size_t) * h_vertices * THREADS_PER_ITERATION);
    cudaMalloc((void **) &d_p_last,  sizeof(edge_size_t)   * h_vertices * THREADS_PER_ITERATION);
    cudaMalloc((void **) &d_p_prev,  sizeof(edge_size_t)   * (h_edges   * THREADS_PER_ITERATION + 1));
    cudaMalloc((void **) &d_p_val,   sizeof(vertex_size_t) * (h_edges   * THREADS_PER_ITERATION + 1));
    cudaMalloc((void **) &d_sigma,   sizeof(edge_size_t)   * h_vertices * THREADS_PER_ITERATION);
    cudaMalloc((void **) &d_depth,   sizeof(vertex_size_t) * h_vertices * THREADS_PER_ITERATION);
    cudaMalloc((void **) &d_visited, sizeof(bool)          * h_vertices * THREADS_PER_ITERATION);
    cudaMalloc((void **) &d_delta,   sizeof(double)        * h_vertices * THREADS_PER_ITERATION);
    cudaMalloc((void **) &d_queue,   sizeof(vertex_size_t) * h_vertices * THREADS_PER_ITERATION);

    // Calculating

    kernel<<<BLOCKS_PER_ITERATION, THREADS_PER_BLOCK>>>(d_vertices, 
                                                        d_edges,
                                                        d_ends, 
                                                        d_indices,
                                                        d_result,
                                                        d_order,
                                                        d_p_last,
                                                        d_p_prev,
                                                        d_p_val,
                                                        d_sigma,
                                                        d_depth,
                                                        d_visited,
                                                        d_delta,
                                                        d_queue);
    cudaDeviceSynchronize();

    // Result

    cudaMemcpy(h_result, 
               d_result, 
               sizeof(double) * h_vertices * THREADS_PER_ITERATION,
               cudaMemcpyDeviceToHost);

    for (vertex_size_t i = 0; i < h_vertices * THREADS_PER_ITERATION; i++) {
        result[i % h_vertices] += h_result[i];
    }

    // Cleaning

    cudaFree(d_order);
    cudaFree(d_p_last);
    cudaFree(d_p_prev);
    cudaFree(d_p_val);
    cudaFree(d_sigma);
    cudaFree(d_depth);
    cudaFree(d_visited);
    cudaFree(d_delta);
    cudaFree(d_queue);

    cudaFree(d_result);
    cudaFree(d_indices);
    cudaFree(d_ends);
    cudaFree(d_edges);
    cudaFree(d_vertices);

    free(h_result);
}
