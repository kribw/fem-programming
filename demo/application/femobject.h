#ifndef FEMOBJECT_H
#define FEMOBJECT_H

#include "node.h"
#include <trianglesystem/gmtrianglesystem.h>

using namespace GMlib;

class FEMObject : TriangleFacets<float> {
private:
  ArrayLX<Node> _nodes; // Array of nodes
  DMatrix<float> _A;    // Stiffness matrix
  DMatrix<float> _b;    // Load vector
  float _f;             // External force

public:
  FEMObject();
  Vector<Vector<float, 2>, 3> findVectors(TSEdge<float> *e);
  Vector<Vector<float, 2>, 3> findVectors(TSTriangle<float> *tr, Node *n);
  void randomTriangulation(int n, float r);
  void regularTriangulation(int m, int n, float r);
  void computation();
  void updateHeight(float f);
};

#endif // FEMOBJECT_H
