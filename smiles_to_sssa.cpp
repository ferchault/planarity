#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/op.h>
#include <openbabel/obconversion.h>
extern "C" {
    #include <qhull/qhull_a.h>
}
using namespace std;

int main(int argc,char **argv)
{
  if(argc<3) {
    cout << "Usage: ProgrameName InputFileName OutputFileName";
    return 1;
  }

  OpenBabel::OBConversion conv;
  conv.SetInFormat("SMI");
  OpenBabel::OBMol mol;

  conv.ReadFile(&mol, argv[1]);
  OpenBabel::OBOp* pOp = OpenBabel::OBOp::FindType("gen3D");
  pOp->Do(dynamic_cast<OpenBabel::OBBase*>(&mol), "4");
  qh_new_qhull(3, mol.NumAtoms(), mol.GetCoordinates(), 0, "qhull s FA", NULL, NULL);
  qh_getarea(qh facet_list);
  cout << qh totvol << endl;
  cout << qh totarea << endl;
  qh_freeqhull(!qh_ALL);
  return 0;
}
