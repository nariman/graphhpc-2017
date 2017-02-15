/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <cstring>       // std::memset
#include <iostream>      // std::cout
#include <unordered_map> // std::unordered_map

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

    edge_size_t* factor;


    // Debug
    // Edges before filtering

    #ifdef DEBUG_SOLUTION
    cout << endl << endl << "BEFORE FILTERING (edges: " << edges << ")" << endl;

    for (vertex_size_t v = 0; v < vertices; v++) {
        cout << "Vertex #" << v << ": ";

        for (edge_size_t e = indices[v]; e < indices[v + 1]; e++) {
            cout << ends[e] << " | ";
        }

        cout << endl;
    }

    cout << endl;
    #endif // DEBUG_SOLUTION


    // Filtering edges

    {
        // Allocations

        unordered_map<vertex_size_t, edge_size_t>** filtered_ends_set =
            new unordered_map<vertex_size_t, edge_size_t>*[vertices];

        #pragma omp parallel for
        for (vertex_size_t v = 0; v < vertices; v++) {
            filtered_ends_set[v] =
                new unordered_map<vertex_size_t, edge_size_t>;
        }

        // Counting equal edges
        // Edge is equal, when:
        // for v in V,
        // for i, j, where indices[v] <= i, j < indices[v + 1],
        // ends[i] == ends[j].
        // In other words, counting multiple edges for each vertex

        #pragma omp parallel for
        for (vertex_size_t v = 0; v < vertices; v++) {
            for (edge_size_t e = indices[v]; e < indices[v + 1]; e++) {
                if (filtered_ends_set[v]->count(ends[e]) > 0) {
                    filtered_ends_set[v]->at(ends[e])++;
                } else {
                    filtered_ends_set[v]->insert({ends[e], 1});
                }
            }
        }

        // Allocations

        edge_size_t* filtered_indices    = new edge_size_t[vertices + 1];
                     filtered_indices[0] = 0;

        // Counting new edges, and number of edges for each vertex

        for (vertex_size_t v = 0; v < vertices; v++) {
            filtered_indices[v + 1] = 
                filtered_indices[v] + filtered_ends_set[v]->size();
        }


        // Allocations

        edge_size_t    filtered_edges = filtered_indices[vertices];

        vertex_size_t* filtered_ends  = new vertex_size_t[filtered_edges];
                       factor         = new   edge_size_t[filtered_edges];

        // New edges

        #pragma omp parallel for
        for (vertex_size_t v = 0; v < vertices; v++) {
            edge_size_t e = filtered_indices[v];

            for (auto i = filtered_ends_set[v]->begin();
                 i != filtered_ends_set[v]->end();
                 i++) {
                filtered_ends[e] = i->first;
                       factor[e] = i->second;
                e++;
            }
        }

        // Swap it!

        edges   = filtered_edges;
        ends    = filtered_ends;
        indices = filtered_indices;

        // Deal with it!

        #pragma omp parallel for
        for (vertex_size_t v = 0; v < vertices; v++) {
            delete filtered_ends_set[v];
        }

        delete[] filtered_ends_set;
    }


    // Debug
    // Edges after filtering

    #ifdef DEBUG_SOLUTION
    cout << endl << "AFTER FILTERING (edges: " << edges << ")" << endl;

    for (vertex_size_t v = 0; v < vertices; v++) {
        cout << "Vertex #" << v << ": ";

        for (edge_size_t e = indices[v]; e < indices[v + 1]; e++) {
            cout << ends[e] << "." << factor[e] << " | ";
        }

        cout << endl;
    }

    cout << endl << endl;
    #endif // DEBUG_SOLUTION


    // Allocations

    vertex_size_t* order   = new vertex_size_t[threads * vertices];
    edge_size_t*   p_last  = new   edge_size_t[threads * vertices];
    edge_size_t*   p_prev  = new   edge_size_t[threads * edges + 1];
    vertex_size_t* p_val   = new vertex_size_t[threads * edges + 1];
    edge_size_t*   p_fact  = new   edge_size_t[threads * edges + 1];
    edge_size_t*   sigma   = new   edge_size_t[threads * vertices];
    vertex_size_t* depth   = new vertex_size_t[threads * vertices];
    bool*          visited = new          bool[threads * vertices];
    double*        delta   = new        double[threads * vertices];
    vertex_size_t* queue   = new vertex_size_t[threads * vertices];

    double*        res     = new        double[threads * vertices];

    p_prev++;
    p_val++;
    p_fact++;

    fill(res, res + threads * vertices, 0.0);


    // Main logic

    #pragma omp parallel firstprivate(order, p_last, p_prev, p_val, p_fact, sigma, depth, visited, delta, queue, res)
    {
        unsigned tid = omp_get_thread_num();

        order   += vertices * tid;
        p_last  += vertices * tid;
        p_prev  += edges    * tid;
        p_val   += edges    * tid;
        p_fact  += edges    * tid;
        sigma   += vertices * tid;
        depth   += vertices * tid;
        visited += vertices * tid;
        delta   += vertices * tid;
        queue   += vertices * tid;

        res     += vertices * tid;

        #pragma omp for
        for (vertex_size_t s = 0; s < vertices; s++) {
            memset(p_last, 0, sizeof(edge_size_t) * vertices);
            memset(sigma,  0, sizeof(edge_size_t) * vertices);

            fill(visited, visited + vertices, false);
            fill(delta,   delta   + vertices, 0.0);

            vertex_size_t* order_pos    = order;
            edge_size_t    p_pos        = 0; // achtung!
            vertex_size_t* queue_front  = queue;
            vertex_size_t* queue_back   = queue;

            sigma[s]      = 1;
            depth[s]      = 0;
            visited[s]    = true;
            *queue_back++ = s;

            while (queue_front != queue_back) {
                vertex_size_t v = *queue_front++;
                *order_pos++ = v;

                for (edge_size_t e = indices[v]; e < indices[v + 1]; e++) {
                    vertex_size_t t = ends[e];

                    if (!visited[t]) {
                        depth[t] = depth[v] + 1;
                        *queue_back++ = t;
                        visited[t] = true;
                    }

                    if (depth[t] > depth[v]) {
                        sigma[t] += sigma[v] * factor[e];

                        p_prev[p_pos] = p_last[t];
                        p_val[p_pos] = v;
                        p_fact[p_pos] = factor[e];
                        p_last[t] = p_pos++;
                    }
                }
            }

            while (--order_pos != order)  {
                vertex_size_t v = *order_pos;
                edge_size_t i = p_last[v];
                double d = (1 + delta[v]) / (double) sigma[v];

                while (i != 0) {
                    delta[p_val[i]] += sigma[p_val[i]] * d * p_fact[i];
                    i = p_prev[i];
                }

                res[v] += delta[v] / 2;
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
    delete[] visited;
    delete[] depth;
    delete[] sigma;
    delete[] --p_fact;
    delete[] --p_val;
    delete[] --p_prev;
    delete[] p_last;
    delete[] order;

    delete[] factor;
    delete[] indices;
    delete[] ends;
}
