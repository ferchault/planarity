#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
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

  omp_set_num_threads(1);

  OpenBabel::OBConversion conv;
  conv.SetInFormat("SDF");
  OpenBabel::OBMol mol;
  OpenBabel::OBOp* gen3d = OpenBabel::OBOp::FindType("gen3D");
  vertexT *vertexA;
  vertexT *vertexB;
  double distance, maxdistance;

  bool notatend = conv.ReadFile(&mol, argv[1]);
  while (notatend) {
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

    cout << "X " << qh totvol << " " << qh totarea << " " << maxdistance << " " << mol.NumAtoms() << endl;
    qh_freeqhull(!qh_ALL);
    notatend = conv.Read(&mol);
  }
  
  return 0;
}
