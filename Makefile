CC=g++
CFLAGS=-c -I dep/dist/include/openbabel-2.0 -O3 -fopenmp
LDFLAGS=-lopenbabel -L dep/dist/lib/ -lqhull -lgomp

all: smiles_to_sssa

smiles_to_sssa: smiles_to_sssa.o
	$(CC) $(LDFLAGS) smiles_to_sssa.o -o smiles_to_sssa

smiles_to_sssa.o: smiles_to_sssa.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) smiles_to_sssa.cpp

clean:
	rm -rf smiles_to_sssa.o smiles_to_sssa
