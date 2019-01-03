CC=g++
CFLAGS=-Wall -g -c -I dep/dist/include/libqhull -I dep/dist/include/openbabel-2.0/ -O3 -fopenmp
LDFLAGS=-L dep/dist/lib/ -lqhull -lgomp -lopenbabel

all: smiles_to_sssa sdf_to_sssa

smiles_to_sssa: smiles_to_sssa.o
	$(CC) $(LDFLAGS) smiles_to_sssa.o -o smiles_to_sssa

smiles_to_sssa.o: smiles_to_sssa.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) smiles_to_sssa.cpp

sdf_to_sssa: sdf_to_sssa.o
	$(CC) $(LDFLAGS) sdf_to_sssa.o -o sdf_to_sssa

sdf_to_sssa.o: sdf_to_sssa.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) sdf_to_sssa.cpp

clean:
	rm -rf smiles_to_sssa.o smiles_to_sssa sdf_to_sssa sdf_to_sssa.o
