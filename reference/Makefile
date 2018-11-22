# defines
MPICXX=mpicxx
SPEEDFLAGS=-O3
CXXFLAGS=-Wall $(SPEEDFLAGS)
LDFLAGS=-lrt $(SPEEDFLAGS)

TARGET = validation gen_RMAT gen_random reference reference_bfs reference_mpi solution solution_mpi

all: $(TARGET)

# your own implementation
solution: main.o solution.o graph_tools.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

# your own distributed implementation
solution_mpi: main_mpi.o solution_mpi.o graph_tools.o gen_RMAT_mpi.o gen_random_mpi.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

# reference implementation with matrices
reference: main.o reference.o graph_tools.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

# reference implementation with bfs
reference_bfs: main.o reference_bfs.o graph_tools.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

# distributed reference implementation
reference_mpi: main_mpi.o reference_mpi.o graph_tools.o gen_RMAT_mpi.o gen_random_mpi.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

# RMAT generator
gen_RMAT: gen_RMAT.o graph_tools.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

# Erdos-Renyi (random) graph generator
gen_random: gen_random.o graph_tools.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

# validation
validation: validation.o reference_bfs.o graph_tools.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

.cpp.o:
	$(MPICXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)
