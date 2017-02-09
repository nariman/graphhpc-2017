/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#ifndef GRAPHHPC_2017_MAIN_H
#define GRAPHHPC_2017_MAIN_H


typedef unsigned vertex_size_t;
typedef unsigned long long edge_size_t;

/**
 * Calculates betweenness centrality by the algorithm.
 * 
 * @param vertices Vertices number
 * @param edges    Edges number
 * @param ends     Edge ends
 * @param indices  Edge indices
 * @param result   Result array
 */
void run(vertex_size_t  vertices, 
         edge_size_t    edges, 
         vertex_size_t* ends, 
         edge_size_t*   indices, 
         double*        result);

#endif // GRAPHHPC_2017_MAIN_H
