#include <iostream>
#include <fstream>
#include <openbabel/mol.h>
#include <openbabel/op.h>
#include <openbabel/obconversion.h>
extern "C" {
    #include <qhull/qhull_a.h>
}
using namespace std;

int main(int argc,char **argv)
{
  if(argc < 2) {
    cout << "Usage:" << argv[0] << "SMIfile";
    return 1;
  }

  ifstream infile(argv[1]);

  OpenBabel::OBConversion conv;
  conv.SetInFormat("SMI");
  OpenBabel::OBMol mol;
  OpenBabel::OBOp* gen3d = OpenBabel::OBOp::FindType("gen3D");

  string line;
  while (getline(infile, line)) {
    conv.ReadString(&mol, line);
    gen3d->Do(dynamic_cast<OpenBabel::OBBase*>(&mol), "4");
    qh_new_qhull(3, mol.NumAtoms(), mol.GetCoordinates(), 0, "qhull s FA", NULL, NULL);
    qh_getarea(qh facet_list);
    cout << line << "  " << qh totvol << " " << qh totarea << endl;
    qh_freeqhull(!qh_ALL);
  }
  
  return 0;
}
