/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <stdio.h> // std::printf

#include <cuda.h>

#include "main.h"

#ifndef BLOCKS_PER_ITERATION
#define BLOCKS_PER_ITERATION 64
#endif

#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 512
#endif

using namespace std;


__device__ __inline__ double atomic_add_double(double* address, 
                                               double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*) address; 
    unsigned long long int old = *address_as_ull, assumed;
    
    do {
        assumed = old;
        old = atomicCAS(address_as_ull,
                        assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);
    
    return __longlong_as_double(old);
}

__global__ void kernel_prepare(vertex_size_t* d_vertices,
                               double*        d_result,
                               double*        d_result_block) {
    vertex_size_t vertices = *d_vertices;

    double* result_block = d_result_block + vertices * blockIdx.x;

    for (vertex_size_t i = threadIdx.x; i < vertices; i += blockDim.x) {
        result_block[i] = 0;
    }

    if (blockIdx.x == 0) {
        for (vertex_size_t i = threadIdx.x; i < vertices; i += blockDim.x) {
            d_result[i] = 0;
        }
    }

    __threadfence_block();
}

__global__ void kernel_calculate_threads_optimum(vertex_size_t* d_vertices,
                                                 edge_size_t*   d_edges,
                                                 vertex_size_t* d_ends,
                                                 edge_size_t*   d_indices) {
    // nothing...
}

__global__ void kernel_summarize(vertex_size_t* d_vertices,
                                 double*        d_result,
                                 double*        d_result_block) {
    vertex_size_t vertices = *d_vertices;

    double* result_block = d_result_block + vertices * blockIdx.x;

    for (vertex_size_t v = threadIdx.x; v < vertices; v += blockDim.x) {
        atomic_add_double(d_result + v, result_block[v]);
    }

    __threadfence_block();
}

__global__ void kernel_debug_graph_data(vertex_size_t* d_vertices,
                                        edge_size_t*   d_edges,
                                        vertex_size_t* d_ends,
                                        edge_size_t*   d_indices) {
    vertex_size_t vertices = *d_vertices;
    edge_size_t   edges    = *d_edges;

    printf("\n\nGPU GRAPH DATA\n");

    printf("Vertices: %d\n", vertices);
    printf("Edges:    %lld\n", edges);

    printf("Ends:     ");
    for (edge_size_t i = 0; i < edges; i++) {
        printf("%d ", d_ends[i]);
    }

    printf("\n");

    printf("Indices:   ");
    for (vertex_size_t i = 0; i < vertices + 1; i++) {
        printf("%lld ", d_indices[i]);
    }

    printf("\n\n");
}

__global__ void kernel_debug_result_data(vertex_size_t* d_vertices,
                                         double*        d_result,
                                         double*        d_result_block) {
    vertex_size_t vertices = *d_vertices;

    printf("\n\nGPU RESULT DATA\n");

    printf("Global:   ");
    for (vertex_size_t v = 0; v < vertices; v++) {
        printf("%f ", d_result[v]);
    }

    for (unsigned block = 0; block < BLOCKS_PER_ITERATION; block++) {
        printf("Block #%d: ", block);

        double* result_block = d_result_block + vertices * block;
        for (vertex_size_t v = 0; v < vertices; v++) {
            printf("%f ", result_block[v]);
        }

        printf("\n");
    }

    printf("\n\n");
}

__global__ void kernel_calculate_bc(vertex_size_t* d_vertices,
                                    edge_size_t*   d_edges,
                                    vertex_size_t* d_ends,
                                    edge_size_t*   d_indices,
                                    double*        d_result_block,
                                    vertex_size_t* d_order,
                                    edge_size_t*   d_p_last,
                                    edge_size_t*   d_p_prev,
                                    vertex_size_t* d_p_val,
                                    edge_size_t*   d_sigma,
                                    vertex_size_t* d_depth,
                                    int*           d_visited,
                                    double*        d_delta,
                                    vertex_size_t* d_queue) {
    vertex_size_t  vertices     = *d_vertices;
    edge_size_t    edges        = *d_edges;

    double*        result_block = d_result_block + vertices * blockIdx.x;

    vertex_size_t* order        = d_order        + vertices * blockIdx.x * 2; // magic
    edge_size_t*   p_last       = d_p_last       + vertices * blockIdx.x;
    edge_size_t*   p_prev       = d_p_prev       + edges    * blockIdx.x + 1;
    vertex_size_t* p_val        = d_p_val        + edges    * blockIdx.x + 1;
    edge_size_t*   sigma        = d_sigma        + vertices * blockIdx.x;
    vertex_size_t* depth        = d_depth        + vertices * blockIdx.x;
    int*           visited      = d_visited      + vertices * blockIdx.x;
    double*        delta        = d_delta        + vertices * blockIdx.x;
    vertex_size_t* queue        = d_queue        + vertices * blockIdx.x;

    __syncthreads();

    for (vertex_size_t s = blockIdx.x; s < vertices; s += gridDim.x) {
        __shared__ vertex_size_t* order_pos;
        __shared__ edge_size_t    p_pos;
        __shared__ vertex_size_t* queue_front;
        __shared__ vertex_size_t* queue_back;
        __shared__ vertex_size_t  queue_back_cumulative;
        __shared__ vertex_size_t  level; // depth, as a constant

        __shared__ vertex_size_t  need;

        __syncthreads();

        for (vertex_size_t i = threadIdx.x; i < vertices; i += blockDim.x) {
            p_last[i]  = 0;
            sigma[i]   = 0;
            depth[i]   = vertices;
            visited[i] = 0;
            delta[i]   = 0;
        }

        __syncthreads();

        if (threadIdx.x == 0) {
            order_pos     = order;
            p_pos         = 0;
            queue_front   = queue;
            queue_back    = queue;
            level         = 0;

            sigma[s]      = 1;
            depth[s]      = 0;
            visited[s]    = 1;
            *queue_back++ = s;
        }

        __syncthreads();

        while (queue_front != queue_back) {
            if (threadIdx.x == 0) {
                queue_back_cumulative = 0;

                need                  = queue_back - queue_front;
            }

            __syncthreads();

            for (vertex_size_t i = threadIdx.x; i < need; i += blockDim.x) {
                vertex_size_t v  = *(queue_front + i);
                *(order_pos + i) = v;

                for (edge_size_t j = d_indices[v]; j < d_indices[v + 1]; j++) {
                    vertex_size_t t = d_ends[j];

                    if (0 == atomicCAS(visited + t, 0, 1)) {
                        atomicExch(depth + t, level + 1);
                        *(queue_back + atomicAdd(&queue_back_cumulative, 
                                                 (vertex_size_t) 1)) = t;
                    }

                    if (depth[t] > level) {
                        atomicAdd(sigma + t, sigma[v]);
                        edge_size_t p_pos_cache = atomicAdd(&p_pos, 1);
                        p_prev[p_pos_cache] = atomicExch(p_last + t, 
                                                         p_pos_cache);
                        p_val[p_pos_cache] = v;
                    }
                }
            }

            __syncthreads();

            if (threadIdx.x == 0) {
                order_pos    += need;
                queue_front  += need;
                queue_back   += queue_back_cumulative;
                level        ++;
                *order_pos++  = need;
            }

            __threadfence_block();
            __syncthreads();
        }

        __syncthreads();
        
        if (threadIdx.x == 0) {
            // Set to zero the first level of depth (where is only `s` vertex)
            order[1] = 0;

            // Set `order_pos` to a valid cell for the reverse iteration
            order_pos--;
        }

        #ifdef DEBUG_SOLUTION
        if (blockIdx.x == 0 && threadIdx.x == 0) {
            printf("\n\n");

            printf("P_LAST:\n");
            for (vertex_size_t i = 0; i < vertices; i++) {
                printf("%d: %lld | ", i, p_last[i]);
            }

            printf("\n\n");

            printf("P_PREV:\n");
            for (edge_size_t i = 0; i < p_pos; i++) {
                printf("%d: %lld | ", i, p_prev[i]);
            }

            printf("\n\n");

            printf("P_VAL:\n");
            for (edge_size_t i = 0; i < p_pos; i++) {
                printf("%d: %lld | ", i, p_val[i]);
            }

            printf("\n\n");

            printf("SIGMA:\n");
            for (edge_size_t i = 0; i < vertices; i++) {
                printf("%d: %lld | ", i, sigma[i]);
            }

            printf("\n\n");

            printf("DEPTH:\n");
            for (edge_size_t i = 0; i < vertices; i++) {
                printf("%d: %lld | ", i, depth[i]);
            }

            printf("\n\n");
        }
        #endif // DEBUG_SOLUTION

        __syncthreads();

        while (true) {
            if (threadIdx.x == 0) {
                need = *order_pos--;
            }

            __threadfence_block();
            __syncthreads();

            if (need == 0) {
                break;
            }

            for (vertex_size_t i = threadIdx.x; i < need; i += blockDim.x) {
                vertex_size_t v = *(order_pos - i);
                edge_size_t j = p_last[v];
                double d = (1 + delta[v]) / (double) sigma[v];

                while (j != 0) {
                    vertex_size_t t = p_val[j];
                    atomic_add_double(delta + t, sigma[t] * d);
                    j = p_prev[j];
                }

                result_block[v] += delta[v] / 2;
            }

            __syncthreads();

            if (threadIdx.x == 0) {
                order_pos -= need;
            }

            __syncthreads();
        }

        __threadfence_block();
        __syncthreads();
    }
}

/**
 * Calculates betweenness centrality by the algorithm with CUDA stack.
 */
void run(vertex_size_t  h_vertices, 
         edge_size_t    h_edges, 
         vertex_size_t* ends, 
         edge_size_t*   indices, 
         double*        h_result) {

    // Allocations for the graph data

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
    cudaMemcpy(d_ends,      ends,       sizeof(vertex_size_t) * h_edges,        cudaMemcpyHostToDevice);
    cudaMemcpy(d_indices,   indices,    sizeof(edge_size_t) * (h_vertices + 1), cudaMemcpyHostToDevice);

    #ifdef DEBUG_SOLUTION
    kernel_debug_graph_data<<<1, 1>>>(d_vertices, 
                                      d_edges, 
                                      d_ends, 
                                      d_indices);
    cudaDeviceSynchronize();
    #endif // DEBUG_SOLUTION


    // Allocating an arrays for results

    double* d_result;
    double* d_result_block;

    cudaMalloc((void **) &d_result, sizeof(double) * h_vertices);
    cudaMalloc((void **) &d_result_block, sizeof(double) * 
                                          h_vertices * 
                                          BLOCKS_PER_ITERATION);

    #ifndef DEBUG_SOLUTION
    cudaMemset(d_result, 0, sizeof(double) * h_vertices);
    cudaMemset(d_result_block, 0, sizeof(double) * 
                                h_vertices *
                                BLOCKS_PER_ITERATION);
    #endif // DEBUG_SOLUTION


    // Allocating an temporary arrays for main calculations

    vertex_size_t* d_order;
    edge_size_t*   d_p_last;
    edge_size_t*   d_p_prev;
    vertex_size_t* d_p_val;
    edge_size_t*   d_sigma;
    vertex_size_t* d_depth;
    int*           d_visited;
    double*        d_delta;
    vertex_size_t* d_queue;

    // Array with a level-ordering.
    // Each level contains some vertex indices, and number of indices,
    // as a last element in a level.
    // This structure allows us to iterate for each level in reverse 
    // ordering.
    cudaMalloc((void **) &d_order,   sizeof(vertex_size_t) * h_vertices * BLOCKS_PER_ITERATION * 2);

    cudaMalloc((void **) &d_p_last,  sizeof(edge_size_t)   * h_vertices * BLOCKS_PER_ITERATION);
    cudaMalloc((void **) &d_p_prev,  sizeof(edge_size_t)   * (h_edges   * BLOCKS_PER_ITERATION + 1));
    cudaMalloc((void **) &d_p_val,   sizeof(vertex_size_t) * (h_edges   * BLOCKS_PER_ITERATION + 1));
    cudaMalloc((void **) &d_sigma,   sizeof(edge_size_t)   * h_vertices * BLOCKS_PER_ITERATION);
    cudaMalloc((void **) &d_depth,   sizeof(vertex_size_t) * h_vertices * BLOCKS_PER_ITERATION);
    cudaMalloc((void **) &d_visited, sizeof(int)           * h_vertices * BLOCKS_PER_ITERATION);
    cudaMalloc((void **) &d_delta,   sizeof(double)        * h_vertices * BLOCKS_PER_ITERATION);
    cudaMalloc((void **) &d_queue,   sizeof(vertex_size_t) * h_vertices * BLOCKS_PER_ITERATION);


    // Calculating

    #ifdef DEBUG_SOLUTION
    kernel_prepare<<<BLOCKS_PER_ITERATION, THREADS_PER_BLOCK>>>(d_vertices,
                                                                d_result, 
                                                                d_result_block);
    cudaDeviceSynchronize();
    #endif // DEBUG_SOLUTION

    kernel_calculate_bc<<<BLOCKS_PER_ITERATION, THREADS_PER_BLOCK>>>(d_vertices, 
                                                                     d_edges,
                                                                     d_ends, 
                                                                     d_indices,
                                                                     d_result_block,
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

    kernel_summarize<<<BLOCKS_PER_ITERATION, THREADS_PER_BLOCK>>>(d_vertices,
                                                                  d_result, 
                                                                  d_result_block);
    cudaDeviceSynchronize();

    // Result

    #ifdef DEBUG_SOLUTION
    kernel_debug_result_data<<<1, 1>>>(d_vertices,
                                       d_result, 
                                       d_result_block);
    cudaDeviceSynchronize();
    #endif // DEBUG_SOLUTION

    cudaMemcpy(h_result, 
               d_result, 
               sizeof(double) * h_vertices,
               cudaMemcpyDeviceToHost);

    // Cleaning

    cudaFree(d_queue);
    cudaFree(d_delta);
    cudaFree(d_visited);
    cudaFree(d_depth);
    cudaFree(d_sigma);
    cudaFree(d_p_val);
    cudaFree(d_p_prev);
    cudaFree(d_p_last);
    cudaFree(d_order);

    cudaFree(d_result_block);
    cudaFree(d_result);

    cudaFree(d_indices);
    cudaFree(d_ends);
    cudaFree(d_edges);
    cudaFree(d_vertices);
}
