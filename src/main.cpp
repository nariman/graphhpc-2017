/*
 * GraphHPC-2017 Contest
 * Betweenness Centrality Problem
 */

#include <iostream> // std::cout, std::fixed
#include <iomanip>  // std::setprecision
#include <cstring>  // std::strcmp
#include <cassert>  // std::assert

#if defined(CLOCK_MONOTONIC)
#define CLOCK CLOCK_MONOTONIC
#elif defined(CLOCK_REALTIME)
#define CLOCK CLOCK_REALTIME
#else
#error "Failed to find a timing clock."
#endif

#define DEFAULT_NUMBER_OF_ITERATIONS 4

using namespace std;


void usage(int argc, char **argv)
{
    cout << "Usage:" << endl;
    cout << "    " << argv[0] << " -in <input filename> [options...]" << endl;
    cout << endl;

    cout << "Options:" << endl;
    cout << "    -in     <input filename>  (required) -- input filename with defined graph" << endl;
    cout << "    -out    <output filename>            -- output filename for results of calculation of betweenness centrality" << endl;
    cout << "    -nIters <iterations>                 -- number of solution iterations" << endl;
    cout << endl;
}

int main(int argc, char **argv) {
    //
    // Initialization
    //

    if (argc == 1) {
        usage(argc, argv);
        return 1;
    }

    string input_filename = "";
    string output_filename = "";
    int    number_of_iterations = DEFAULT_NUMBER_OF_ITERATIONS;

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
    
    if (!input_filename.length() || number_of_iterations < 1) {
        usage(argc, argv);
        return 1;
    }
    
    if (!output_filename.length()) {
        output_filename = input_filename + ".res";
    }

    //
    // Reading
    //

    //
    // Calculation information
    //

    cout << "Input filename:       " << input_filename << endl;
    cout << "Output filename:      " << output_filename << endl;
    cout << "Number of iterations: " << number_of_iterations << endl;
    cout << endl;

    cout << "Vertices: " << 1 << endl;
    cout << "Edges:    " << 1 << endl;
    cout << endl;

    //
    // Calculation
    //

    struct timespec start_time, finish_time;
    double* iteration_times = new double[number_of_iterations];
    assert(iteration_times != NULL);

    cout << "Starting calculation:" << endl;

    for (int i = 1; i <= number_of_iterations; i++) {
        cout << " - " << i << " of " << number_of_iterations;

        clock_gettime(CLOCK, &start_time);
        // algorithm
        clock_gettime(CLOCK, &finish_time);

        double time = finish_time.tv_sec - (double) start_time.tv_sec;
        time = time + time * 1.0e-9;
        iteration_times[i - 1] = time;

        cout << ": " << setprecision(5) << fixed << time << " sec." << endl;
    }

    cout << endl;

    //
    // Performance information
    //

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

    //
    // Collecting
    //

    delete[] iteration_times;

    // 
    // Done
    //

    return 0;
}
