#ifndef NODE_H
#define NODE_H

#include <trianglesystem/gmtrianglesystem.h>

using namespace GMlib;

class Node {
private:
  TSVertex<float> *_v; // Node vertex

public:
  Node();
  Node(TSVertex<float> *v);

  void setZ(float z);
  Array<TSTriangle<float> *> getTriangles();
  TSEdge<float> *getNeighbor(Node &n);
  bool isThis(TSVertex<float> *v);
};

#endif // NODE_H
