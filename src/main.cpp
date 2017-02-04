/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <algorithm> // std::log
#include <cassert>   // std::assert
#include <cstring>   // std::strcmp
#include <iomanip>   // std::setprecision
#include <iostream>  // std::cout, std::fixed

#if defined(CLOCK_MONOTONIC)
#define CLOCK CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define CLOCK CLOCK_REALTIME
#else
#error "Failed to find a timing clock."
#endif

#define DEFAULT_NUMBER_OF_ITERATIONS 4

using namespace std;


typedef unsigned vertex_size_t;
typedef unsigned long long edge_size_t;

class Graph {
public:
    /**
     * Scale of the graph (log2 of number of vertices).
     */
    vertex_size_t scale;

    /**
     * Number of vertices.
     */
    vertex_size_t vertices;

    /**
     * Number of edges.
     */
    edge_size_t   edges;

    /**
     * Edge ends.
     */
    vertex_size_t* ends;

    /**
     * Edge indices.
     * Rule is: for each v (vertex in graph)
     *  - ends[index[v]] is a first edge (edge end) from v,
     *  - ends[index[v + 1] - 1] is a last edge (edge end) from v.
     */
    edge_size_t* indices;

    /**
     * Create a new graph with provided numbers of vertices and edges.
     *
     * @param vertices Nubmer of vertices
     * @param edges    Number of edges
     * @param ends     Edge ends
     * @param indices  Edge indices
     */
    Graph(vertex_size_t vertices, edge_size_t edges, vertex_size_t* ends,
          edge_size_t* indices) {
        this->scale = log(vertices) / log(2);

        this->vertices = vertices;
        this->edges = edges;

        this->ends = ends;
        this->indices = indices;
    }

    /**
     * Delete a graph and free memory, allocated to the edges.
     */
    ~Graph() {
        delete[] this->ends;
        delete[] this->indices;
    }
};

/**
 * Prints usage to the standart i/o.
 */
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

int main(int argc, char** argv) {
    //
    // Initialization
    //

    if (argc == 1) {
        usage(argc, argv); // first arguments is a executable filename
        return 1;
    }

    string input_filename;
    string output_filename;
    int    number_of_iterations = DEFAULT_NUMBER_OF_ITERATIONS;

    Graph* graph;

    double* iteration_times = new double[number_of_iterations];
    double* result;

    assert(iteration_times != NULL);

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
        FILE *f = fopen(input_filename.c_str(), "rb");
        assert(f != NULL);

        vertex_size_t vertices;
        edge_size_t edges;
        
        fread(&vertices, sizeof(vertex_size_t), 1, f);
        fread(&edges, sizeof(edge_size_t), 1, f);
        
        unsigned char align; fread(&align, sizeof(unsigned char), 1, f);
        
        edge_size_t* indices = new edge_size_t[vertices + 1];
        vertex_size_t* ends = new vertex_size_t[edges];

        fread(indices, sizeof(edge_size_t), vertices + 1, f);
        fread(ends, sizeof(vertex_size_t), edges, f);

        fclose(f);

        graph = new Graph(vertices, edges, ends, indices);
    }

    //
    // Calculation information
    //

    cout << "Input filename:       " << input_filename << endl;
    cout << "Output filename:      " << output_filename << endl;
    cout << "Number of iterations: " << number_of_iterations << endl;
    cout << endl;

    cout << "Scale:    " << graph->scale << endl;
    cout << "Vertices: " << graph->vertices << endl;
    cout << "Edges:    " << graph->edges << endl;

    cout << endl;

    //
    // Calculation
    //

    cout << "Calculation:" << endl;

    for (int i = 1; i <= number_of_iterations; i++) {
        cout << " - " << i << " of " << number_of_iterations;

        struct timespec start_time, finish_time;

        result = new double[graph->vertices];
        for (vertex_size_t i = 0; i < graph->vertices; result[i++] = 0);

        clock_gettime(CLOCK, &start_time);

        {
            // algorithm
        }

        clock_gettime(CLOCK, &finish_time);

        double time = finish_time.tv_sec - (double) start_time.tv_sec;
        time = time + time * 1.0e-9;
        iteration_times[i - 1] = time;

        cout << ": " << setprecision(5) << fixed << time << " sec." << endl;
    }

    cout << endl;

    //
    // Answer writing
    //

    {
        FILE *f = fopen(output_filename.c_str(), "wb");
        assert(f != NULL);

        fwrite(result, sizeof(double), graph->vertices, f);

        fclose(f);
    }

    //
    // Performance information
    //

    {
        double min_time, max_time, all_time, avg_time;
        min_time = max_time = avg_time = iteration_times[0];

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

        cout << "Time = " << avg_time << " sec." << endl;
    }

    //
    // Collecting
    //

    delete[] result;
    delete[] iteration_times;
    delete graph;

    // 
    // Done
    //

    return 0;
}
