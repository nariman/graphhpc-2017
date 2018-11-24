CC          := g++
CCMPI       := mpicxx
CCCUDA      := nvcc
CFLAGS      := -fopenmp -O3 -std=c++11 -Wall -Wextra # -DDEBUG_SOLUTION
CFLAGSCUDA  := -x cu --gpu-architecture=compute_20 --Wno-deprecated-gpu-targets # -DDEBUG_SOLUTION

SRC         := src
REFERENCE   := reference
BIN         := bin
BUILD       := build

INCS        := -I $(SRC)
REFINCS     := -I $(REFERENCE)
LIBS        := -L/usr/local/cuda/lib64 -lcudart -lrt

EXT         := cpp
EXTCUDA     := cu


all: solution-openmp solution-cuda solution-mixed
solution: solution-mixed


solution-mixed: $(BUILD)/main.o $(BUILD)/mixed.o
	@mkdir -p $(BIN)
	$(CC) $^ -o $(BIN)/solution-mixed $(CFLAGS) $(LIBS) $(INCS)

solution-openmp: $(BUILD)/main.o $(BUILD)/openmp.o
	@mkdir -p $(BIN)
	$(CC) $^ -o $(BIN)/solution-openmp $(CFLAGS) $(INCS)

solution-cuda: $(BUILD)/main.o $(BUILD)/cuda.o
	@mkdir -p $(BIN)
	$(CC) $^ -o $(BIN)/solution-cuda $(CFLAGS) $(LIBS) $(INCS)


gen_rmat: $(BUILD)/gen_rmat.o $(BUILD)/graph_tools.o
	@mkdir -p $(BIN)
	$(CCMPI) $^ -o $(BIN)/gen_rmat $(CFLAGS) $(INCS)

gen_random: $(BUILD)/gen_random.o $(BUILD)/graph_tools.o
	@mkdir -p $(BIN)
	$(CCMPI) $^ -o $(BIN)/gen_random $(CFLAGS) $(INCS)

validation: $(BUILD)/validation.o $(BUILD)/reference_bfs.o $(BUILD)/graph_tools.o
	@mkdir -p $(BIN)
	$(CCMPI) $^ -o $(BIN)/validation $(CFLAGS) $(INCS)


$(BUILD)/main.o: $(SRC)/main.$(EXT)
	@mkdir -p $(BUILD)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

$(BUILD)/mixed.o: $(SRC)/mixed.$(EXTCUDA)
	@mkdir -p $(BUILD)
	$(CCCUDA) -c -o $@ $< -std=c++11 --compiler-options "$(CFLAGS)" $(CFLAGSCUDA) $(INCS)

$(BUILD)/openmp.o: $(SRC)/openmp.$(EXT)
	@mkdir -p $(BUILD)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

$(BUILD)/cuda.o: $(SRC)/cuda.$(EXTCUDA)
	@mkdir -p $(BUILD)
	$(CCCUDA) -c -o $@ $< $(CFLAGSCUDA) $(INCS)


$(BUILD)/graph_tools.o: $(REFERENCE)/graph_tools.$(EXT)
	@mkdir -p $(BUILD)
	$(CCMPI) -c -o $@ $< $(CFLAGS) $(REFINCS)

$(BUILD)/gen_rmat.o: $(REFERENCE)/gen_RMAT.$(EXT)
	@mkdir -p $(BUILD)
	$(CCMPI) -c -o $@ $< $(CFLAGS) $(REFINCS)

$(BUILD)/gen_random.o: $(REFERENCE)/gen_random.$(EXT)
	@mkdir -p $(BUILD)
	$(CCMPI) -c -o $@ $< $(CFLAGS) $(REFINCS)

$(BUILD)/validation.o: $(REFERENCE)/validation.$(EXT)
	@mkdir -p $(BUILD)
	$(CCMPI) -c -o $@ $< $(CFLAGS) $(REFINCS)

$(BUILD)/reference_bfs.o: $(REFERENCE)/reference_bfs.$(EXT)
	@mkdir -p $(BUILD)
	$(CCMPI) -c -o $@ $< $(CFLAGS) $(REFINCS)


clean:
	@rm -r $(BUILD)/*
	@rm $(BIN)/*

.PHONY: clean
