CC			:= g++-4.8
CCCUDA		:= nvcc
CFLAGS		:= -fopenmp -O3 -std=c++11 -Wall -Wextra
CFLAGSCUDA	:= -x cu

SRC			:= src
BIN			:= bin
BUILD		:= build

INCS		:= -I $(SRC)

EXT			:= cpp
EXTCUDA		:= cu

all: solution-openmp solution-cuda
solution: solution-openmp

solution-openmp: $(BUILD)/main.o $(BUILD)/openmp.o
	@mkdir -p $(BIN)
	$(CC) $^ -o $(BIN)/solution-openmp $(CFLAGS) $(INCS)

solution-cuda: $(BUILD)/main.o $(BUILD)/cuda.o
	@mkdir -p $(BIN)
	$(CC) $^ -o $(BIN)/solution-cuda -L/usr/local/cuda/lib64 -lcudart -lrt $(CFLAGS) $(INCS)

$(BUILD)/main.o: $(SRC)/main.cpp
	@mkdir -p $(BUILD)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

$(BUILD)/openmp.o: $(SRC)/openmp.cpp
	@mkdir -p $(BUILD)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS)

$(BUILD)/cuda.o: $(SRC)/cuda.cu
	@mkdir -p $(BUILD)
	$(CCCUDA) -c -o $@ $< $(CFLAGSCUDA) $(INCS)

clean:
	@rm -r $(BUILD)/*
	@rm $(BIN)/solution-openmp
	@rm $(BIN)/solution-cuda

.PHONY: clean
