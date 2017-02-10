/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <cassert>
#include <cstring>  // std::memset
#include <iostream> // std::cout
#include <set>      // std::set

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

    // restructuring vertices

    set<vertex_size_t>* vei = new set<vertex_size_t>;

    for (vertex_size_t v = 0; v < vertices; v++) {
        vei->insert(v);
    }


    // restructuring edges

    multiset<vertex_size_t>** nei = new multiset<vertex_size_t>*[vertices];

    #pragma omp parallel for
    for (vertex_size_t v = 0; v < vertices; v++) {
        nei[v] = new multiset<vertex_size_t>;
    }

    #pragma omp parallel for
    for (vertex_size_t v = 0; v < vertices; v++) {
        for (vertex_size_t i = indices[v]; i < indices[v + 1]; i++) {
            nei[v]->insert(ends[i]);
        }
    }


    // allocations

    edge_size_t*  partial           = new edge_size_t[vertices];
    vertex_size_t vertices_adaptive = vertices;
    edge_size_t   edges_adaptive    = edges;

    memset(partial, 0, sizeof(edge_size_t) * vertices);


    // leafs precomputation

    {
        // allocations

        vertex_size_t* queue   = new vertex_size_t[vertices];
        vertex_size_t* queue_front = queue;
        vertex_size_t* queue_back  = queue;

        // first step

        #pragma omp parallel for
        for (vertex_size_t v = 0; v < vertices; v++) {
            if (nei[v]->size() < 2) {
                #pragma omp critical
                {
                    *queue_back++ = v;
                }
            }
        }

        // next steps

        while (queue_front != queue_back) {
            vertex_size_t v = *queue_front++;
            vertex_size_t u = *(nei[v]->begin());

            result[u] += 2 *
                         (vertices - partial[v] - partial[u] - 2) *
                         (partial[v] + 1);

            partial[u] += partial[v] + 1;

            vei->erase(v);
            nei[u]->erase(v);
            edges_adaptive -= 2;

            if (nei[u]->size() < 2) {
                *queue_back++ = u;
            }
        }

        delete[] queue;
    }

    vertices_adaptive = vei->size();


    // main computation

    if (vertices > 1) {
        unsigned threads = omp_get_max_threads();

        // constructing list of vertices

        vertex_size_t* vertices_list_adaptive = 
            new vertex_size_t[vertices_adaptive];
        
        {
            vertex_size_t* pos = vertices_list_adaptive;

            for (auto i = vei->begin(); i != vei->end(); i++) {
                *pos++ = *i;
            }
        }

        // allocations

        vertex_size_t* order   = new vertex_size_t[threads * vertices_adaptive];
        edge_size_t*   p_last  = new   edge_size_t[threads * vertices_adaptive];
        edge_size_t*   p_prev  = new   edge_size_t[threads * edges_adaptive + 1];
        vertex_size_t* p_val   = new vertex_size_t[threads * edges_adaptive + 1];
        edge_size_t*   sigma   = new   edge_size_t[threads * vertices_adaptive];
        vertex_size_t* depth   = new vertex_size_t[threads * vertices_adaptive];
        bool*          visited = new          bool[threads * vertices_adaptive];
        double*        delta   = new        double[threads * vertices_adaptive];
        vertex_size_t* queue   = new vertex_size_t[threads * vertices_adaptive];

        p_prev++;
        p_val++;

        // main logic

        #pragma omp parallel firstprivate(order, p_last, p_prev, p_val, sigma, depth, visited, delta, queue)
        {
            // offset calculation

            unsigned tid = omp_get_thread_num();

            order   += vertices_adaptive * tid;
            p_last  += vertices_adaptive * tid;
            p_prev  += edges_adaptive    * tid;
            p_val   += edges_adaptive    * tid;
            sigma   += vertices_adaptive * tid;
            depth   += vertices_adaptive * tid;
            visited += vertices_adaptive * tid;
            delta   += vertices_adaptive * tid;
            queue   += vertices_adaptive * tid;

            // main cycle

            #pragma omp for nowait
            for (vertex_size_t i = 0; i < vertices_adaptive; i++) {
                vertex_size_t s = vertices_list_adaptive[i];

                memset(p_last, 0, sizeof(edge_size_t) * vertices_adaptive);
                memset(sigma,  0, sizeof(edge_size_t) * vertices_adaptive);

                fill(visited, visited + vertices_adaptive, false);
                fill(delta,   delta   + vertices_adaptive, 0.0);

                vertex_size_t* order_pos   = order;
                vertex_size_t  p_pos       = 0;
                vertex_size_t* queue_front = queue;
                vertex_size_t* queue_back  = queue;

                sigma[s]      = 1;
                depth[s]      = 0;
                visited[s]    = true;
                *queue_back++ = s;

                while (queue_front != queue_back) {
                    vertex_size_t v = *queue_front++;
                    *order_pos++ = v;

                    for (auto j = nei[v]->begin(); j != nei[v]->end(); j++) {
                        vertex_size_t t = *j;

                        if (!visited[t]) {
                            depth[t] = depth[v] + 1;
                            *queue_back++ = t;
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

                while (--order_pos != order)  {
                    vertex_size_t v = *order_pos;
                    edge_size_t i = p_last[v];
                    double d = (1 + delta[v] + partial[v]) / (double) sigma[v];

                    while (i != 0) {
                        vertex_size_t t = p_val[i];
                        delta[t] += sigma[t] * d;
                        i = p_prev[i];
                    }

                    #pragma omp atomic
                    result[v] += delta[v] * (partial[s] + 1) / 2;
                }
            }
        }

        // de

        delete[] queue;
        delete[] delta;
        delete[] visited;
        delete[] depth;
        delete[] sigma;
        delete[] --p_val;
        delete[] --p_prev;
        delete[] p_last;
        delete[] order;

        delete[] vertices_list_adaptive;
    }


    // de

    delete[] partial;

    #pragma omp parallel for
    for (vertex_size_t v = 0; v < vertices; v++) {
        delete nei[v];
    }

    delete[] nei;
    delete vei;
}
