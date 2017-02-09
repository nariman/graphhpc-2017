/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <algorithm> // std::log
#include <cassert>   // std::assert
#include <cstring>   // std::memset, std::strcmp
#include <iomanip>   // std::setprecision
#include <iostream>  // std::cout, std::fixed

#include "main.h"

#if defined(CLOCK_MONOTONIC)
#define CLOCK CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define CLOCK CLOCK_REALTIME
#else
#error "Failed to find a timing clock."
#endif

#define DEFAULT_NUMBER_OF_ITERATIONS 4

using namespace std;

/**
 * Prints usage to the standart i/o.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void usage(int argc, char** argv)
{
    cout << "Usage:" << endl;
    cout << "    " << argv[0] << " -in <input filename> [options...]" << endl;
    cout << endl;

    cout << "Options:" << endl;
    cout << "    -in     <input filename>  (required) -- input filename with defined graph" << endl;
    cout << "    -out    <output filename>            -- output filename for result of calculation of betweenness centrality" << endl;
    cout << "    -nIters <iterations>                 -- number of solution iterations, positive number" << endl;
    cout << endl;
}
#pragma GCC diagnostic pop

int main(int argc, char** argv) {
    //
    // Initialization
    //

    if (argc == 1) {
        usage(argc, argv); // first argument is an executable filename
        return 1;
    }

    /**
     * Input filename with given graph.
     */
    string               input_filename;

    /**
     * Output filename for result of calculation.
     */
    string               output_filename;

    /**
     * Number of solution iterations.
     */
    int                  number_of_iterations = DEFAULT_NUMBER_OF_ITERATIONS;

    /**
     * Scale of the graph (log2 of number of vertices).
     */
    vertex_size_t        scale;

    /**
     * Number of vertices.
     */
    vertex_size_t        vertices;

    /**
     * Number of edges.
     */
    edge_size_t          edges;

    /**
     * Edge ends.
     */
    vertex_size_t*       ends;

    /**
     * Edge indices.
     * Rule is: for each v (vertex in graph)
     *  - ends[index[v]] is a first edge (edge end) from v,
     *  - ends[index[v + 1] - 1] is a last edge (edge end) from v.
     */
    edge_size_t*         indices;

    /**
     * Array with times of each iteration of the solution.
     */
    double*              iteration_times = new double[number_of_iterations];

    /**
     * Result of calculation.
     */
    double*              result;

    //
    // Configuration
    //
    
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-in")) {
            input_filename = argv[++i];
        }

        if (!strcmp(argv[i], "-out")) {
            output_filename = argv[++i];
        }

        if (!strcmp(argv[i], "-nIters")) {
            number_of_iterations = (int) atoi(argv[++i]);
        }
    }

    // input filename is required, number of iterations must be positive
    if (!input_filename.length() || number_of_iterations < 1) {
        usage(argc, argv);
        return 1;
    }

    // output filename can be not specified
    if (!output_filename.length()) {
        output_filename = input_filename + ".res";
    }

    //
    // Reading
    //

    {
        // #pragma GCC diagnostic push
        // #pragma GCC diagnostic ignored "-Wunused-result"

        FILE *f = fopen(input_filename.c_str(), "rb");
        assert(f != NULL);

        assert(fread(&vertices, sizeof(vertex_size_t), 1, f) == 1);
        assert(fread(&edges, sizeof(edge_size_t), 1, f) == 1);

        scale = log(vertices) / log(2);
        
        unsigned char align;
        assert(fread(&align, sizeof(unsigned char), 1, f) == 1);
        
        indices = new edge_size_t[vertices + 1];
        ends = new vertex_size_t[edges];

        assert(fread(indices, sizeof(edge_size_t), vertices + 1, f) == 
               vertices + 1);
        assert(fread(ends, sizeof(vertex_size_t), edges, f) == edges);

        fclose(f);

        // #pragma GCC diagnostic pop
    }

    //
    // Calculation information
    //

    cout << "Input filename:       " << input_filename << endl;
    cout << "Output filename:      " << output_filename << endl;
    cout << "Number of iterations: " << number_of_iterations << endl;
    cout << endl;

    cout << "Scale:    " << scale << endl;
    cout << "Vertices: " << vertices << endl;
    cout << "Edges:    " << edges << endl;

    cout << endl;

    //
    // Calculation
    //

    result = new double[vertices];

    cout << "Calculation:" << endl;

    for (int it = 1; it <= number_of_iterations; it++) {
        cout << " - " << it << " of " << number_of_iterations;

        struct timespec start_time, finish_time;

        fill(result, result + vertices, (double) 0.0);

        clock_gettime(CLOCK, &start_time);
        run(vertices, edges, ends, indices, result);
        clock_gettime(CLOCK, &finish_time);

        double time = finish_time.tv_sec - (double) start_time.tv_sec +
            (finish_time.tv_nsec - (double) start_time.tv_nsec) * 1.0e-9;
        iteration_times[it - 1] = time;

        cout << ": " << setprecision(5) << fixed << time << " sec." << endl;
    }

    cout << endl;

    //
    // Result writing
    //

    {
        FILE *f = fopen(output_filename.c_str(), "wb");
        assert(f != NULL);

        assert(fwrite(result, sizeof(double), vertices, f) == vertices);

        fclose(f);
    }

    //
    // Performance information
    //

    {
        double min_time, max_time, all_time, avg_time;
        min_time = max_time = all_time = avg_time = iteration_times[0];

        for (int i = 1; i < number_of_iterations; i++) {
            all_time += iteration_times[i];

            if (iteration_times[i] < min_time) {
                min_time = iteration_times[i];
            }

            if (iteration_times[i] > max_time) {
                max_time = iteration_times[i];
            }
        }
        
        avg_time = all_time / number_of_iterations;
        
        cout << "Iteration times:" << endl;

        cout << " - Min: " << min_time << " sec." << endl;
        cout << " - Avg: " << avg_time << " sec." << endl;
        cout << " - Max: " << max_time << " sec." << endl;

        cout << endl;

        // For the contest judge system
        cout << "avg = " << avg_time << " sec." << endl;
    }

    //
    // Collecting
    //

    delete[] result;
    delete[] iteration_times;
    delete[] ends;
    delete[] indices;

    // 
    // Done
    //

    return 0;
}
