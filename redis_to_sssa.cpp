#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <signal.h>
#include <omp.h>

#define USE_RDKIT

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

extern "C" {
	#include <qhull_a.h>
}
extern "C" { 	
	#include "hiredis/hiredis.h" 
};

using namespace std;

redisContext *context;

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
	reply = (redisReply*)redisCommand(context, "RPOP SSSASMILES");
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
				coords = new double[numatoms*3];
				const RDGeom::POINT3D_VECT &vec = mol->getConformer(0).getPositions();
				for (int i = 0; i < numatoms; i++) {
					coords[i*3 + 0]  = vec[i].x;
					coords[i*3 + 1]  = vec[i].y;
					coords[i*3 + 2]  = vec[i].z;
				}
				qh_new_qhull(3, numatoms, coords, 0, "qhull s FA QJ Pp", NULL, NULL);
				delete[] coords;
				delete mol_ro;
			} catch ( ... ) {
				continue;
			}
		#endif

		qh_getarea(qh facet_list);
		cout << qh totvol << " " << qh totarea << " " << " " << numatoms << endl;

		qh_freeqhull(!qh_ALL);
	}

	redisFree(context);
	return 0;
}
