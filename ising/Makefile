#SHELL := /bin/bash
CC=gcc-7
NVCC=nvcc
CFLAGS=  -O3  -std=c99 -lm
CUDAFLAGS= --compiler-options -Wall
EXE=$(SDIR)/isingV0 $(SDIR)/isingV1 $(SDIR)/isingV2 $(SDIR)/isingV3
LIBS=$(LDIR)/isingV0.a $(LDIR)/isingV1.a $(LDIR)/isingV2.a $(LDIR)/isingV3.a
SDIR=./src
LDIR=./lib
IDIR=./inc
MAIN=tester.cu
INCLUDES= -I$(IDIR)
LIBRARIES= -L$(LDIR)


all:	$(EXE)

lib:	$(LIBS)

$(SDIR)/isingV0: $(SDIR)/isingV0.c
	$(CC) $(CFLAGS) -o $@ $^

$(SDIR)/isingV1: $(LDIR)/isingV1.a $(SDIR)/$(MAIN)
	$(NVCC) $(INCLUDES) -o $@ $^ $(CUDAFLAGS)

$(SDIR)/isingV2: $(LDIR)/isingV2.a $(SDIR)/$(MAIN)
	$(NVCC) $(INCLUDES) -o $@ $^ $(CUDAFLAGS)


$(SDIR)/isingV3: $(LDIR)/isingV3.a $(SDIR)/$(MAIN)
	$(NVCC) $(INCLUDES) -o $@ $^ $(CUDAFLAGS)




$(LDIR)/%.a: $(SDIR)/%.o
	ar rcs $@ $^


$(SDIR)/isingV1.o: $(SDIR)/isingV1.cu
	$(NVCC) -o $@ -c $^ $(CUDAFLAGS)

$(SDIR)/isingV2.o: $(SDIR)/isingV2.cu
	$(NVCC) -o $@ -c $^ $(CUDAFLAGS)

$(SDIR)/isingV3.o: $(SDIR)/isingV3.cu
	$(NVCC) -o $@ -c $^ $(CUDAFLAGS)




clean:
	rm -f $(SDIR)/*.o $(EXE) $(LIBS)
