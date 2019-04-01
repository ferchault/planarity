#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <signal.h>
#include <omp.h>

// if present, use rdkit instead of OpenBabel
#define USE_RDKIT

// if present, calculate moments of inertia instead
#define INERTIA

// Set to 1 to disable pipelining
#define REDIS_PIPELINING 50

#ifdef USE_OBABEL
	#include <openbabel/mol.h>
	#include <openbabel/op.h>
	#include <openbabel/obconversion.h>
#endif

#ifdef USE_RDKIT
	#include <GraphMol/GraphMol.h>
	#include <GraphMol/MolOps.h>
	#include <GraphMol/SmilesParse/SmilesParse.h>
	#include <GraphMol/FileParsers/MolSupplier.h>
	#include <GraphMol/FileParsers/MolWriters.h>
	#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
	#include <GraphMol/FileParsers/FileParsers.h>
	#include <GraphMol/DistGeomHelpers/Embedder.h>
#endif

#ifndef INERTIA
extern "C" {
	#include <qhull_a.h>
}
#endif
extern "C" { 	
	#include "hiredis/hiredis.h" 
};

using namespace std;

redisContext *context;
int remaining_replies = 0;

redisContext * redis_connect() {
	redisContext * context = redisConnect(getenv("CONFORMERREDIS"), 80);
	if (context != NULL && context->err) {
		printf("Error: %s\n", context->errstr);
	}
	redisReply *reply;
	reply = (redisReply*)redisCommand(context, "AUTH chemspacelab");
	freeReplyObject(reply);

	return context;
}

// reads a single SMILES
int redis_fetch(redisContext * context, string * line) {
	redisReply *reply;
	if (REDIS_PIPELINING > 1) {
		if (remaining_replies == 0) {
			for (int i=0; i < REDIS_PIPELINING; ++i) {
				redisAppendCommand(context, "RPOP SSSASMILES");
			}
			remaining_replies = REDIS_PIPELINING;
		}
		redisGetReply(context, (void**)&reply);
		remaining_replies--;
	} else {
		reply = (redisReply*)redisCommand(context, "RPOP SSSASMILES");
	}
	
	if (reply == NULL || reply->type != REDIS_REPLY_STRING) {
		// invalid reply || no further work or server-side issue
		freeReplyObject(reply);
		return 2;
	}

	// valid reply
	if (reply->len == 0) {
		return 1;
	}
	line->assign(reply->str, reply->len);
	freeReplyObject(reply);
	return 0;
}

void cleanup(int s) {
	redisFree(context);
	exit(1);
}

int main(int argc,char **argv)
{
	struct sigaction sigIntHandler;

	sigIntHandler.sa_handler = cleanup;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);

	omp_set_num_threads(1);
	context = redis_connect();

	#ifdef USE_OBABEL
		OpenBabel::OBConversion conv;
		conv.SetInFormat("SMI");
		OpenBabel::OBMol mol;
		OpenBabel::OBOp* gen3d = OpenBabel::OBOp::FindType("gen3D");
	#endif

	int numatoms;
	double * coords;
	double in[3];
	double com[3];
	double mass, totalmass, x, y, z;
	string line;
	while (redis_fetch(context, &line) == 0) {
		#ifdef USE_OBABEL
			conv.ReadString(&mol, line);
			gen3d->Do(dynamic_cast<OpenBabel::OBBase*>(&mol), "4");
			qh_new_qhull(3, mol.NumAtoms(), mol.GetCoordinates(), 0, "qhull s FA QJ Pp", NULL, NULL);
			numatoms = mol.NumAtoms():
		#endif
		#ifdef USE_RDKIT
			try {
				RDKit::ROMol *mol_ro = RDKit::SmilesToMol(line);
				RDKit::RWMOL_SPTR mol( new RDKit::RWMol( *mol_ro ) );
				RDKit::MolOps::addHs(*mol);
				numatoms = mol->getNumAtoms();
				RDKit::DGeomHelpers::EmbedMolecule( *mol );
				RDKit::MMFF::MMFFOptimizeMolecule(*mol, 100, "MMFF94s");
				const RDGeom::POINT3D_VECT &vec = mol->getConformer(0).getPositions();
				#ifndef INERTIA
					coords = new double[numatoms*3];
					for (int i = 0; i < numatoms; i++) {
						coords[i*3 + 0]  = vec[i].x;
						coords[i*3 + 1]  = vec[i].y;
						coords[i*3 + 2]  = vec[i].z;
					}
					qh_new_qhull(3, numatoms, coords, 0, "qhull s FA QJ Pp", NULL, NULL);
				#endif
				#ifdef INERTIA
					totalmass = 0.;
					com[0] = 0.;
					com[1] = 0.;
					com[2] = 0.;
					for (int i = 0; i < numatoms; i++) {
						mass = mol->getAtomWithIdx()->getMass();
						totalmass += mass;
						com[0] += vec[i].x*mass;
						com[1] += vec[i].y*mass;
						com[2] += vec[i].z*mass;
					}	
					com[0] /= totalmass;
					com[1] /= totalmass;
					com[2] /= totalmass;
					in[0] = 0;
					in[1] = 0;
					in[2] = 0;
					for (int i = 0; i < numatoms; i++) {
						mass = mol->getAtomWithIdx()->getMass();
						x = vec[i].x - com[0];
						y = vec[i].y - com[1];
						z = vec[i].z - com[2];
						in[0] += mass * (y*y  + z*z);
						in[1] += mass * (z*z  + x*x);
						in[2] += mass * (y*y  + x*x);
					}
				#endif
				delete[] coords;
				delete mol_ro;
			} catch ( ... ) {
				continue;
			}
		#endif
		
		#ifndef INERTIA
			qh_getarea(qh facet_list);
			cout << qh totvol << " " << qh totarea << " " << " " << numatoms << endl;

			qh_freeqhull(!qh_ALL);
		#endif
		#ifdef INERTIA
			cout << in[0] << " " << in[1] << " " << in[2] << endl;
		#endif
	}

	redisFree(context);
	return 0;
}
