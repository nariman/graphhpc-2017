CC			:= g++-4.8
CCCUDA		:= nvcc
CFLAGS		:= -fopenmp -O3 -std=c++11 -Wall -Wextra # -DDEBUG_SOLUTION
CFLAGSCUDA	:= -x cu --gpu-architecture=compute_20 --Wno-deprecated-gpu-targets # -DDEBUG_SOLUTION

SRC			:= src
BIN			:= bin
BUILD		:= build

INCS		:= -I $(SRC)
LIBS		:= -L/usr/local/cuda/lib64 -lcudart -lrt

EXT			:= cpp
EXTCUDA		:= cu

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

clean:
	@rm -r $(BUILD)/*
	@rm $(BIN)/solution-mixed
	@rm $(BIN)/solution-openmp
	@rm $(BIN)/solution-cuda

.PHONY: clean
