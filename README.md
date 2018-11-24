# GraphHPC-2017 Contest
Implementation for the GraphHPC-2017 Contest  
Betweenness Centrality Problem

### Links
* [GraphHPC Conference Page](http://www.dislab.org/GraphHPC-2017/en/)
* [GraphHPC Contests Page](http://contest.dislab.org/)
  * [GraphHPC-2017 Contest Page](http://contest.dislab.org/archive/2017)

## Disclaimer

This code does not claim to be the fastest. I had no much time as well as
experience in OpenMP and CUDA to write performance efficient solution.

## Reference files

You can find the reference files (testgen, sample solution) in the `reference`
folder.
Take a note, that the reference files are licensed under GNU GPLv2.

## Data format

Sorry for the data format, but it was defined by contest organizers. If you need
other format than that is described below, you have to code your own parser or
converter.

#### Input data

An input format is a binary-encoded file with vertices and edges numbers, align,
and lists of indices and edge ends.  
A graph structure is pretty usual for contests: there is an array of out
edges' ends (other vertices), and array of indices. Vertex's out edges can be
determined by the first array's slice with the following rule:

```
Rule is: for each v (vertex in graph)
 - ends[index[v]] is a first edge (edge end) from v,
 - ends[index[v + 1] - 1] is a last edge (edge end) from v.
```

The graph is directed.

If the binary-encoded input format is not you are looking for, you can write
your own parser. See the `Reading` section in the `main.cpp` file.

#### Output data

An output format is a binary-encoded file too with an array of calculated
betweenness centrality value for each vertex.

See the `Result writing` section in the `main.cpp` file, if you need to change
this behaviour.

## Running

### Compilation

Depending on which solution you need to run, you can choose the suitable
command to compile:

```bash
# This solution uses only OpenMP interface
$ make solution-openmp

# Only CUDA
$ make solution-cuda

# This is aт attempt to use OpenMP and CUDA in parallel
$ make solution-mixed
# Alias
$ make solution

# You can compile all of them at once
$ make all
```

To compile CUDA and mixed versions, make sure you've installed CUDA Toolkit.

#### Test generation and validation tools

_You can skip this section, if you have your own tests._

To compile these tools, make sure you have OpenMPI libraries installed on your
system.

Compile the test generation tools:

```bash
# For RMAT graph generation algorithm
$ make gen_rmat

# For random generation
$ make gen_random
```

There's a validation tool, you can use it too, but it is not supposed to be used
on large graphs due to inefficient algorithm.

```bash
$ make validation
```

### Tests generation

_You can skip this section, if you have your own tests._

To run the solution, you need to generate some tests first:

```bash
# We'll place tests in the /tests folder
$ mkdir tests
$ cd tests

# Option -s is an exponent of the base 2 - the graph size
$ ../bin/gen_rmat -s 5
$ ../bin/gen_random -s 10
```

Generated tests will be saved in the current working directory.

### Running

The solution binaries is located in the `bin` directory.  
To run the solution, provide the input file path with an `-in` argument:

```bash
# You can run solution-cuda or solution-mixed too.
# We placed our test in the /tests folder
$ ./bin/solution-openmp -nIters 5 -in ./tests/random-10
```

Answer will be saved in the same directory where test is located under the same
name with a `.res` suffix.

**You can run any solution or reference tool without any arguments to print the
help.**

That's all.

## License

MIT.  
See LICENSE file for more information.
