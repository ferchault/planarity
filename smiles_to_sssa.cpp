#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/op.h>
#include <openbabel/obconversion.h>
using namespace std;

int main(int argc,char **argv)
{
  if(argc<3)
  {
    cout << "Usage: ProgrameName InputFileName OutputFileName";
    return 1;
  }
  
  OpenBabel::OBConversion conv;
  conv.SetInFormat("SMI");
  OpenBabel::OBMol mol;

  conv.ReadFile(&mol, argv[1]);
  OpenBabel::OBOp* pOp = OpenBabel::OBOp::FindType("gen3D");
  pOp->Do(dynamic_cast<OpenBabel::OBBase*>(&mol), "4");
  cout << mol.GetCoordinates() << endl;
  return 0;
}
