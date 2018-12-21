CC=g++
CFLAGS=-c -I /usr/include/openbabel-2.0 -O3
LDFLAGS=-lopenbabel -L /usr/lib/openbabel-2.3.2/lib -lqhull

all: smiles_to_sssa

smiles_to_sssa: smiles_to_sssa.o
	$(CC) $(LDFLAGS) smiles_to_sssa.o -o smiles_to_sssa

smiles_to_sssa.o: smiles_to_sssa.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) smiles_to_sssa.cpp

clean:
	rm -rf smiles_to_sssa.o smiles_to_sssa
