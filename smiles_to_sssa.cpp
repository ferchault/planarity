#include <iostream>
#include <fstream>
#include <cmath>
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
  vertexT *vertexA;
  vertexT *vertexB;
  double distance, maxdistance;

  string line;
  while (getline(infile, line)) {
    conv.ReadString(&mol, line);
    gen3d->Do(dynamic_cast<OpenBabel::OBBase*>(&mol), "4");
    qh_new_qhull(3, mol.NumAtoms(), mol.GetCoordinates(), 0, "qhull s FA", NULL, NULL);
    qh_getarea(qh facet_list);

    // Calculate largest distance
    maxdistance = 0;
    for (vertexA = qh vertex_list; vertexA && vertexA->next; vertexA = vertexA->next) {
      for (vertexB = vertexA->next; vertexB && vertexB->next; vertexB = vertexB->next) {
        distance = 0;
        for (int dim = 0; dim < 3; dim++) {
          distance += (vertexA->point[dim] - vertexB->point[dim])*(vertexA->point[dim] - vertexB->point[dim]);
        }
        if (distance > maxdistance) {
          maxdistance = distance;
        }
      }
    }
    maxdistance = sqrt(maxdistance);

    cout << line << "  " << qh totvol << " " << qh totarea << " " << maxdistance << endl;
    qh_freeqhull(!qh_ALL);
  }
  
  return 0;
}
