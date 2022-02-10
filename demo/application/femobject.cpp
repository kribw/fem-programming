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
    std::swap(v[0], v[2]);
    std::swap(v[1], v[2]);
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

  this->triangulateDelaunay();
  // Not sure if correct to call this here
}

void FEMObject::regularTriangulation(int n, int m, float r) {
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

      Point<float, 2> p = R * Vector<float, 2>(r * (j + 1) / m, 0.0);
      this->insertAlways(TSVertex<float>(p));
    }
  }

  this->triangulateDelaunay();
}

void FEMObject::computation() {
  // Fill in the array of _nodes
  for (int i = 0; i < this->getSize(); i++) {

    // If the point does not belong to the boundary
    if (!(*this)[i].boundary()) {
      TSVertex<float> *vertex = this->getVertex(i);
      _nodes += Node(vertex);
    }
  }

  // Create [n x n] zero matrix (_A)
  _A = DMatrix<float>(this->getSize());

  // Compute elements of the stiffness matrix
  for (int i = 0; i < _nodes.size(); i++) {
    for (int j = 0; j < i; j++) {
      TSEdge<float> *edge = _nodes[i].getNeighbor(_nodes[j]);

      if (edge == NULL) {
        _A[i][j] = _A[j][i] = 0;
        continue;
      }

      Vector<Vector<float, 2>, 3> vectors = findVectors(edge);
      auto dd = 1 / (vectors[0] * vectors[0]);
      auto dh1 = dd * vectors[1] * vectors[0];
      auto dh2 = dd * vectors[2] * vectors[0];
      auto area1 = std::abs(vectors[0] ^ vectors[1]);
      auto area2 = std::abs(vectors[0] ^ vectors[2]);
      auto h1 = dd * area1 * area1;
      auto h2 = dd * area2 * area2;
      _A[i][j] = _A[j][i] = (dh1 * (1 - dh1) / h1 - dd) * area1 / 2 +
                            (dh2 * (1 - dh2) / h2 - dd) * area2 / 2;
    }
  }

  // Diagonal element
  Array<TSTriangle<float> *> triangles = _nodes[i].getTriangles();
  float s_tr = 0;
  for (int i = 0; i < triangles.size(); i++) {
    Vector<Vector<float, 2>, 3> vectors = findVectors(triangles[k], &_nodes[i]);
    _A[i][i] +=
        (vectors[2] * vectors[2]) / (2 * std::abs(vectors[0] * vectors[1]));
    s_tr += triangles[i]->getArea2D();
  }

  // Compute load vector
  for (int i = 0; i < _b.getDim1(); ++i) {
    _b[i] = (1 / 3) *;
  }

  // Update heights using the coefficients (xi)
  updateHeight(_f);
}

void FEMObject::updateHeight(float f) {
  DVector<float> x = DMatrix<float>(_A.invert() * (f * _b)).toDVector();

  for (int i = 0; i < _nodes.size(); i++) {
    _nodes[i].setZ(x[i]);
  }
}
