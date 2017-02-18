/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <cstring>  // std::memset
#include <iomanip>  // std::setprecision
#include <iostream> // std::cout, std::fixed

#ifdef DEBUG_SOLUTION
#include <map>      // std::map
#endif

#include <cuda.h>
#include <omp.h>

#include "main.h"

#if defined(CLOCK_MONOTONIC)
#define CLOCK CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define CLOCK CLOCK_REALTIME
#else
#error "Failed to find a timing clock."
#endif

#ifndef BLOCKS_PER_ITERATION
#define BLOCKS_PER_ITERATION 64
#endif

#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 256
#endif

#define BPI BLOCKS_PER_ITERATION
#define TPB THREADS_PER_BLOCK

#define CPU
#define GPU

using namespace std;


#ifdef DEBUG_SOLUTION
static map<string, timespec> timers;
static map<string, double> timings;
#endif

static inline void timer_start(const string &name) {
#ifdef DEBUG_SOLUTION
    clock_gettime(CLOCK, &timers[name]);
#endif
}

static inline void timer_end(const string &name) {
#ifdef DEBUG_SOLUTION
    timespec finish_time;
    clock_gettime(CLOCK, &finish_time);

    timespec start_time = timers[name];

    double time = finish_time.tv_sec - (double) start_time.tv_sec + 
        (finish_time.tv_nsec - (double) start_time.tv_nsec) * 1.0e-9;
    timings[name] += time;
#endif
}

#ifdef GPU
static __device__ __inline__ double atomic_add_double(double* address, 
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
#endif

#ifdef GPU
__global__ void kernel(vertex_size_t* d_vertices,
                       edge_size_t*   d_edges,
                       vertex_size_t* d_ends,
                       edge_size_t*   d_indices,
                       vertex_size_t* d_order,
                       edge_size_t*   d_p_last,
                       edge_size_t*   d_p_prev,
                       vertex_size_t* d_p_val,
                       edge_size_t*   d_sigma,
                       vertex_size_t* d_depth,
                       int*           d_visited,
                       double*        d_delta,
                       vertex_size_t* d_queue,
                       double*        d_result,
                       double*        d_result_block,
                       vertex_size_t  gpu_start,
                       vertex_size_t  gpu_end) {
    // Offsetting
    
    vertex_size_t  vertices     = *d_vertices;
    edge_size_t    edges        = *d_edges;

    vertex_size_t* order        = d_order        + vertices * blockIdx.x * 2; // magic
    edge_size_t*   p_last       = d_p_last       + vertices * blockIdx.x;
    edge_size_t*   p_prev       = d_p_prev       + edges    * blockIdx.x + 1;
    vertex_size_t* p_val        = d_p_val        + edges    * blockIdx.x + 1;
    edge_size_t*   sigma        = d_sigma        + vertices * blockIdx.x;
    vertex_size_t* depth        = d_depth        + vertices * blockIdx.x;
    int*           visited      = d_visited      + vertices * blockIdx.x;
    double*        delta        = d_delta        + vertices * blockIdx.x;
    vertex_size_t* queue        = d_queue        + vertices * blockIdx.x;

    double*        result_block = d_result_block + vertices * blockIdx.x;

    __syncthreads();

    // Parallel cycle

    for (vertex_size_t s = gpu_start + blockIdx.x; 
         s < gpu_end;
         s += gridDim.x) {
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

                result_block[v] += delta[v];
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

    __threadfence_block();
    __syncthreads();

    // Result reduction

    for (vertex_size_t v = threadIdx.x; v < vertices; v += blockDim.x) {
        atomic_add_double(d_result + v, result_block[v]);
    }

    __threadfence_block();
}
#endif

/**
 * Calculates betweenness centrality by the algorithm with mixed
 * OpenMP-CUDA stack.
 */
void run(vertex_size_t  vertices,
         edge_size_t    edges,
         vertex_size_t* ends,     // don't touch the ends array!
         edge_size_t*   indices,  // don't touch the indices array!
         double*        result) { // be careful!

    //
    // Pre-init
    //

    #ifdef GPU
    timer_start("GPU Pre-initialization");
    cudaSetDevice(0);
    timer_end("GPU Pre-initialization");
    #endif


    //
    // Constants
    //

    // CPU

    #ifdef CPU

    unsigned threads = omp_get_max_threads();

    vertex_size_t* order;
    edge_size_t*   p_last;
    edge_size_t*   p_prev;
    vertex_size_t* p_val;
    edge_size_t*   sigma;
    vertex_size_t* depth;
    double*        delta;
    vertex_size_t* queue;

    double*        result_thread;

    #endif

    #ifdef GPU
    double*        h_result;
    #endif

    // GPU

    #ifdef GPU

    vertex_size_t* d_vertices;
    edge_size_t*   d_edges;
    vertex_size_t* d_ends;
    edge_size_t*   d_indices;

    vertex_size_t* d_order;
    edge_size_t*   d_p_last;
    edge_size_t*   d_p_prev;
    vertex_size_t* d_p_val;
    edge_size_t*   d_sigma;
    vertex_size_t* d_depth;
    int*           d_visited;
    double*        d_delta;
    vertex_size_t* d_queue;

    double*        d_result;
    double*        d_result_block;

    #endif


    //
    // Allocations
    //

    timer_start("Allocations");

    // # linter-disable: length-limit
    #pragma omp parallel sections
    {
        // CPU
        #pragma omp section
        {
            #ifdef CPU
            order         = new vertex_size_t[threads * vertices];
            p_last        = new   edge_size_t[threads * vertices];
            p_prev        = new   edge_size_t[threads * (edges + 1)];
            p_val         = new vertex_size_t[threads * (edges + 1)];
            sigma         = new   edge_size_t[threads * vertices];
            depth         = new vertex_size_t[threads * vertices];
            delta         = new        double[threads * vertices];
            queue         = new vertex_size_t[threads * vertices];

            result_thread = new        double[threads * vertices];
            #endif

            #ifdef GPU
            h_result      = new        double[vertices];
            #endif

            #ifdef CPU
            memset(result_thread, 0, sizeof(double) * threads * vertices);
            #endif
        }

        // GPU
        #ifdef GPU
        #pragma omp section
        {
            cudaMalloc((void **) &d_vertices, sizeof(vertex_size_t));
            cudaMalloc((void **) &d_edges,    sizeof(edge_size_t));
            cudaMalloc((void **) &d_ends,     sizeof(vertex_size_t) * edges);
            cudaMalloc((void **) &d_indices,  sizeof(edge_size_t)   * (vertices + 1));

            // Array with a level-ordering.
            // Each level contains some vertex indices, and number of indices,
            // as a last element in a level.
            // This structure allows us to iterate for each level in reverse 
            // ordering.
            cudaMalloc((void **) &d_order,   sizeof(vertex_size_t) * vertices * BPI * 2);

            cudaMalloc((void **) &d_p_last,  sizeof(edge_size_t)   * vertices * BPI);
            cudaMalloc((void **) &d_p_prev,  sizeof(edge_size_t)   * (edges   * BPI + 1));
            cudaMalloc((void **) &d_p_val,   sizeof(vertex_size_t) * (edges   * BPI + 1));
            cudaMalloc((void **) &d_sigma,   sizeof(edge_size_t)   * vertices * BPI);
            cudaMalloc((void **) &d_depth,   sizeof(vertex_size_t) * vertices * BPI);
            cudaMalloc((void **) &d_visited, sizeof(int)           * vertices * BPI);
            cudaMalloc((void **) &d_delta,   sizeof(double)        * vertices * BPI);
            cudaMalloc((void **) &d_queue,   sizeof(vertex_size_t) * vertices * BPI);

            cudaMalloc((void **) &d_result,       sizeof(double) * vertices);
            cudaMalloc((void **) &d_result_block, sizeof(double) * vertices * BPI);

            cudaMemcpy(d_vertices, &vertices, sizeof(vertex_size_t),                  cudaMemcpyHostToDevice);
            cudaMemcpy(d_edges,    &edges,    sizeof(edge_size_t),                    cudaMemcpyHostToDevice);
            cudaMemcpy(d_ends,      ends,     sizeof(vertex_size_t) * edges,          cudaMemcpyHostToDevice);
            cudaMemcpy(d_indices,   indices,  sizeof(edge_size_t)   * (vertices + 1), cudaMemcpyHostToDevice);

            cudaMemset(d_result,       0, sizeof(double) * vertices);
            cudaMemset(d_result_block, 0, sizeof(double) * vertices * BPI);
        }
        #endif
    }

    timer_end("Allocations");


    //
    // Partitioning
    //

    // Splits for CPU 2/3 and GPU 1/3 of all vertices

    #if defined(CPU) && defined(GPU)
    vertex_size_t cpu_start = 0;
    vertex_size_t cpu_end   = vertices / 3 * 2;
    vertex_size_t gpu_start = vertices / 3 * 2;
    vertex_size_t gpu_end   = vertices;
    #elif defined(CPU)
    vertex_size_t cpu_start = 0
    vertex_size_t cpu_end   = vertices;
    #elif defined(GPU)
    vertex_size_t gpu_start = 0;
    vertex_size_t gpu_end   = vertices;
    #else
    #error "What a strange code, that does nothing? :("
    #endif


    //
    // Main logic
    //

    #ifdef GPU
    timer_start("GPU main logic");
    kernel<<<BLOCKS_PER_ITERATION, THREADS_PER_BLOCK>>>(d_vertices, 
                                                        d_edges,
                                                        d_ends, 
                                                        d_indices,
                                                        d_order,
                                                        d_p_last,
                                                        d_p_prev,
                                                        d_p_val,
                                                        d_sigma,
                                                        d_depth,
                                                        d_visited,
                                                        d_delta,
                                                        d_queue,
                                                        d_result,
                                                        d_result_block,
                                                        gpu_start,
                                                        gpu_end);
    #endif

    // GPU is working async, that means, we can calculate some data on CPU

    #ifdef CPU
    timer_start("CPU main logic");

    #pragma omp parallel firstprivate(order, p_last, p_prev, p_val, sigma, depth, delta, queue, result_thread)
    {
        // Constants

        unsigned tid = omp_get_thread_num();

        // Offsetting

        order         += vertices * tid;
        p_last        += vertices * tid;
        p_prev        += edges    * tid + 1;
        p_val         += edges    * tid + 1;
        sigma         += vertices * tid;
        depth         += vertices * tid;
        delta         += vertices * tid;
        queue         += vertices * tid;

        result_thread += vertices * tid;

        // Variables

        vertex_size_t* order_pos;
        edge_size_t    p_pos;
        vertex_size_t* queue_front;
        vertex_size_t* queue_back;

        // Parallel cycle

        #pragma omp for
        for (vertex_size_t s = cpu_start; s < cpu_end; ++s) {
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

                result_thread[*order_pos] += delta[*order_pos];
            }
        }

        // Result reduction

        for (vertex_size_t v = 0; v < vertices; v++) {
            #pragma omp atomic
            result[v] += result_thread[v];
        }
    }

    timer_end("CPU main logic");
    #endif

    // Waiting for GPU
    #ifdef GPU
    cudaDeviceSynchronize();
    timer_end("GPU main logic");
    #endif


    //
    // Result calculation
    //

    #ifdef GPU
    timer_start("GPU result reduction");

    cudaMemcpy(h_result, 
               d_result, 
               sizeof(double) * vertices,
               cudaMemcpyDeviceToHost);

    #pragma omp parallel for
    for (vertex_size_t v = 0; v < vertices; v++) {
        result[v] += h_result[v];
    }

    timer_end("GPU result reduction");
    #endif

    #pragma omp parallel for
    for (vertex_size_t v = 0; v < vertices; v++) {
        result[v] /= 2;
    }


    //
    // Cleaning
    // Deal with it!
    //

    timer_start("Cleaning");

    #pragma omp parallel sections
    {
        // CPU
        #pragma omp section
        {
            #ifdef GPU
            delete[] h_result;
            #endif

            #ifdef CPU
            delete[] result_thread;

            delete[] queue;
            delete[] delta;
            delete[] depth;
            delete[] sigma;
            delete[] p_val;
            delete[] p_prev;
            delete[] p_last;
            delete[] order;
            #endif
        }

        // GPU
        #ifdef GPU
        #pragma omp section
        {
            cudaDeviceReset();

            /*
            cudaFree(d_result_block);
            cudaFree(d_result);

            cudaFree(d_queue);
            cudaFree(d_delta);
            cudaFree(d_visited);
            cudaFree(d_depth);
            cudaFree(d_sigma);
            cudaFree(d_p_val);
            cudaFree(d_p_prev);
            cudaFree(d_p_last);
            cudaFree(d_order);

            cudaFree(d_indices);
            cudaFree(d_ends);
            cudaFree(d_edges);
            cudaFree(d_vertices);
            */
        }
        #endif
    }

    timer_end("Cleaning");


    //
    // Debug / Timings
    //

    #ifdef DEBUG_SOLUTION

    cout << endl;

    for (auto &it : timings) {
        cout << " - " << it.first << ": ";
        cout << setprecision(5) << fixed << it.second << " sec.";
        cout << endl;
    }

    timers.clear();
    timings.clear();

    #endif
}
