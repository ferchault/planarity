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
  #include <GraphMol/FileParsers/MolSupplier.h>
  #include <GraphMol/FileParsers/MolWriters.h>
  #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
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
        RDKit::MMFF::MMFFOptimizeMolecule(*mol, 100, "MMFF94s");
        numatoms = mol->getNumAtoms();
        qh_new_qhull(3, mol.NumAtoms(), mol->getConformer()->getPositions(), 0, "qhull s FA QJ Pp", NULL, NULL);
      } catch( RDKit::MolSanitizeException &e ) {
        continue;
      }

      delete mol;
    #endif

    qh_getarea(qh facet_list);
    cout << qh totvol << " " << qh totarea << " " << " " << numatoms << endl;

    #ifdef USE_RDKIT
      delete mol;
      delete mol_ro;
    #endif

    qh_freeqhull(!qh_ALL);
  }
  
  redisFree(context);
  return 0;
}
