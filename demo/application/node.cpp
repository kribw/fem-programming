#include "node.h"

Node::Node() {}
Node::Node(TSVertex<float> *v) { _v = v; }

void Node::setZ(float z) { _v->setZ(z); }
Array<TSTriangle<float> *> Node::getTriangles() { return _v->getTriangles(); }

TSEdge<float> *Node::getNeighbor(Node &n) {
  Array<TSEdge<float> *> edge = _v->getEdges();
  for (int i = 0; i < edge.size(); i++) {
    if (n.isThis(edge[i]->getOtherVertex(*_v))) {
      return edge[i];
    }
  }
  return NULL;
}

bool Node::isThis(TSVertex<float> *v) { return v == _v; }
