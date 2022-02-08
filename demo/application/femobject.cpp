#include "femobject.h"

FEMObject::FEMObject() {}

// 1: Intersection where two nodes have a common edge
Vector<Vector<float, 2>, 3> FEMObject::findVectors(TSEdge<float> *e) {
  Array<TSTriangle<float> *> tr = e->getTriangle();   // 2 triangles
  Array<TSVertex<float> *> v1 = tr[0]->getVertices(); // 3 vertices
  Array<TSVertex<float> *> v2 = tr[1]->getVertices(); // 3 vertices
  Point<float, 2> p0, p1, p2, p3;
  p0 = e->getFirstVertex()->getParameter(); // Returns vertex coordinates
  p1 = e->getLastVertex()->getParameter();

  for (int i = 0; i < 3; i++) {
    if (v1[i] != e->getFirstVertex() && v1[i] != e->getLastVertex()) {
      p2 = v1[i]->getParameter();
    }

    if (v2[i] != e->getFirstVertex() && v2[i] != e->getLastVertex()) {
      p3 = v2[i]->getParameter();
    }
  }

  Vector<Vector<float, 2>, 3> d;
  d[0] = p1 - p0;
  d[1] = p2 - p0;
  d[2] = p3 - p0;
  return d;
}

// 2: Intersection of the node with itself
Vector<Vector<float, 2>, 3> FEMObject::findVectors(TSTriangle<float> *tr,
                                                   Node *n) {
  Point<float, 2> p0, p1, p2;
  Array<TSVertex<float> *> v = tr->getVertices();

  if (n->isThis(v[1])) {
    std::swap(v[0], v[1]);
    std::swap(v[1], v[2]);
  }

  if (n->isThis(v[2])) {
    std::swap(v[0], v[1]);
    std::swap(v[1], v[2]);
    // Need to check this
  }

  p0 = v[0]->getParameter();
  p1 = v[1]->getParameter();
  p2 = v[2]->getParameter();

  Vector<Vector<float, 2>, 3> d;
  d[0] = p1 - p0;
  d[1] = p2 - p0;
  d[2] = p2 - p1;
  return d;
}

void FEMObject::randomTriangulation(int n, float r) {
  // n = number of points within boundary
  // r = radius

  // Based on regular triangulation
  // Boundary points are preserved
  // Avoid empty areas

  // std::swap(v[i], v[j]);
  // random "indicies"

  this->insertAlways(TSVertex<float>(Point<float, 2>(0.0, 0.0)));

  // Foreach ring
  //  for (int j = 0; j < m; j++) {

  //    // Foreach point in each ring
  //    for (int i = 0; i < n * (j + 1); i++) {
  //      Angle alpha = (i * M_2PI) / (n * (j + 1));

  //      // Rotation matrix
  //      SqMatrix<float, 2> R = SqMatrix<float, 2>(alpha);

  //      Point<float, 2> p = R * Vector<float, 2>(r * (j + 1) / m + 1);
  //      this->insertAlways(TSVertex<float>(p));
  //    }
  //  }

  triangulateDelaunay();
  // Not sure if correct to call this here
}

void FEMObject::regularTriangulation(int m, int n, float r) {
  // m = number of rings
  // n = number of points within boundary
  // r = radius
  // Insert outer ring
  this->insertAlways(TSVertex<float>(Point<float, 2>(0.0, 0.0)));

  // Foreach ring
  for (int j = 0; j < m; j++) {

    // Foreach point in each ring
    for (int i = 0; i < n * (j + 1); i++) {
      Angle alpha = (i * M_2PI) / (n * (j + 1));

      // Rotation matrix
      SqMatrix<float, 2> R = SqMatrix<float, 2>(alpha);

      Point<float, 2> p = R * Vector<float, 2>(r * (j + 1) / m + 1);
      this->insertAlways(TSVertex<float>(p));
    }
  }

  triangulateDelaunay();
  // Not sure if correct to call this here
}

void FEMObject::computation() {
  for (int i = 0; i < _nodes.size(); i++) {

    // If the point does not belong to the boundary
    //      if(Node::isThis(_no))
    //    _nodes += Node((*this[i]));

    _A = DMatrix<float>(0, 0);
  }

  for (int i = 0; i < _nodes.size(); i++) {
    for (int j = 0; j < i; j++) {
      TSEdge<float> *edge = _nodes[i].getNeighbor(_nodes[j]);

      if (edge != NULL) {

        // compute non-diagonal element of the stiffness matrix
      }
    }
  }

  //  Array<TSTriangle<float> *> triangles = _nodes[i].getTriangles();
  //  for (int i = 0; i < triangles.size(); i++) {
  //    // compute diagonal element of the stiffness matrix
  //    //  Then compute the load vector in the similar way.
  //  }
}

void FEMObject::updateHeight(float f) {
  DVector<float> x = DMatrix<float>(_A.invert() * (f * _b)).toDVector();

  for (int i = 0; i < _nodes.size(); i++) {
    _nodes[i].setZ(x[i]);
  }
}
