#include "femobject.h"
#include <cstdlib>
#include <random>

FEMObject::FEMObject() {}

// 1: Intersection where two nodes have a common edge
Vector<Vector<float, 2>, 3> FEMObject::findVectors(TSEdge<float> *e) {
  Array<TSTriangle<float> *> tr = e->getTriangle();   // 2 triangles
  Array<TSVertex<float> *> v1 = tr[0]->getVertices(); // 3 vertices
  Array<TSVertex<float> *> v2 = tr[1]->getVertices(); // 3 vertices
  Point<float, 2> p0, p1, p2, p3;
  p0 = e->getFirstVertex()->getParameter(); // Returns vertex coordinates
  p1 = e->getLastVertex()->getParameter();

  for (int i = 0; i < 3; ++i) {
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
  // !! Minimum 2x2 for this to even run
  if (n < 2)
    n = 2;

  // m = number of rings
  int m = n + 2;
  if (n <= 3)
    m = 4;

  // ------------- !! Triangulation start !! -------------

  // Insert center point
  insertAlways(TSVertex<float>(Point<float, 2>(0.0, 0.0)));

  std::vector<Point<float, 2>> points{};
  std::vector<Point<float, 2>> boundary_points{};

  // Foreach ring
  for (int j = 0; j < m; j++) {

    // Foreach point in each ring
    for (int i = 0; i < n * (j + 1); i++) {
      Angle alpha = (i * M_2PI) / (n * (j + 1));

      // Rotation matrix
      SqMatrix<float, 2> R = SqMatrix<float, 2>(alpha);
      Point<float, 2> p = R * Vector<float, 2>(r * (j + 1) / m, 0.0);

      points.emplace_back(p);
    }
  }
  // ------------- !! Triangulation end !! -------------

  for (int i = points.size() - 1; boundary_points.size() < n * m; --i) {
    boundary_points.emplace_back(points[i]);
    points.pop_back();
  }

  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  std::uniform_int_distribution<> distr(0, points.size() - 1);

  int t = (n * m * 100);
  for (int i = 0; i < t; ++i) {
    int j = distr(gen);
    int k = distr(gen);
    std::swap(points[j], points[k]);
  }

  int internal_pts = 0;

  for (int i = (n * m) - n; i > 0; i -= n) {
    internal_pts += n;
  }

  for (int i = 0; i < internal_pts; ++i) {
    insertAlways(TSVertex<float>(points[i]));
  }

  for (auto const &point : boundary_points) {
    insertAlways(TSVertex<float>(point));
  }

  triangulateDelaunay();
}

void FEMObject::regularTriangulation(int n, int m, float r) {
  // n = number of boundary points
  // m = number of rings
  // r = radius

  // Insert center point
  insertAlways(TSVertex<float>(Point<float, 2>(0.0, 0.0)));

  // Foreach ring
  for (int j = 0; j < m; j++) {

    // Foreach point in each ring
    for (int i = 0; i < n * (j + 1.0); i++) {
      Angle alpha = (i * M_2PI) / (n * (j + 1.0));

      // Rotation matrix
      SqMatrix<float, 2> R = SqMatrix<float, 2>(alpha);

      Point<float, 2> p = R * Vector<float, 2>(r * (j + 1.0) / m, 0.0);
      insertAlways(TSVertex<float>(p));
    }
  }

  triangulateDelaunay();
}

void FEMObject::computation() {
  // Fill in the array of _nodes
  for (int i = 0; i < size(); i++) {

    // Insert node if point is not a boundary point
    TSVertex<float> *vertex = getVertex(i);
    if (!vertex->boundary()) {
      _nodes += Node(vertex);
    }
  }

  // Create [n x n] zero matrix (_A) and load vector (_b)
  int size = _nodes.size();
  _A = DMatrix<float>(size, size);
  _b = DMatrix<float>(size, 1);

  // Compute non-diagonal elements
  for (int i = 0; i < _nodes.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      TSEdge<float> *edge = _nodes[i].getNeighbor(_nodes[j]);

      if (edge == NULL) {
        _A[i][j] = _A[j][i] = 0.0;
        continue;
      }

      auto const vectors = findVectors(edge);
      auto const &d = vectors[0];
      auto const dd = 1.0 / (d * d);
      auto const &a1 = vectors[1];
      auto const &a2 = vectors[2];
      auto const dh1 = dd * (a1 * d);
      auto const dh2 = dd * (a2 * d);
      auto const area1 = std::abs(d ^ a1);
      auto const area2 = std::abs(d ^ a2);
      auto const h1 = dd * (area1 * area1);
      auto const h2 = dd * (area2 * area2);

      _A[i][j] = _A[j][i] = (dh1 * (1.0 - dh1) / h1 - dd) * area1 / 2.0 +
                            (dh2 * (1.0 - dh2) / h2 - dd) * area2 / 2.0;
    }
  }

  // Compute diagonal elements
  for (int i = 0; i < _nodes.size(); ++i) {
    Array<TSTriangle<float> *> triangles = _nodes[i].getTriangles();
    float s_tk = 0.0;

    for (int k = 0; k < triangles.size(); ++k) {
      Vector<Vector<float, 2>, 3> vectors =
          findVectors(triangles[k], &_nodes[i]);
      s_tk +=
          (vectors[2] * vectors[2]) / (2.0 * std::abs(vectors[0] ^ vectors[1]));
    }

    _A[i][i] = s_tk;
    continue;
  }

  // Compute load vector
  for (int i = 0; i < _nodes.size(); ++i) {
    Array<TSTriangle<float> *> triangles = _nodes[i].getTriangles();
    float s_tr = 0.0;

    for (int k = 0; k < triangles.size(); ++k) {
      s_tr += triangles[k]->getArea2D();
    }

    *_b[i].getPtr() = _f * (1.0 / 3.0) * s_tr;
  }
}

void FEMObject::updateHeight(float f) {
  DVector<float> x = (_A.invert() * (f * _b)).toDVector();

  for (int i = 0; i < _nodes.size(); ++i) {
    _nodes[i].setZ(x[i]);
  }
}

void FEMObject::setForce(float f) { _f = f; }

void FEMObject::demoRegular() {
  setForce(2.5);
  regularTriangulation(5, 4, 1.0);
  computation();
}

void FEMObject::demoRandom() {
  setForce(2.5);
  randomTriangulation(5, 1);
  computation();
}

void FEMObject::localSimulate(double dt) {
  static double t = 0;
  t += dt;

  GMlib::Vector<float, 3> vec(sin(t), cos(t), 5 * dt);
  vec *= 0.05;
  GMlib::Angle angle = GMlib::Angle(1);
  //    this->move(vec);
  //    this->rotate(GMlib::Angle(1), GMlib::Vector<float, 3>(0.0f, 1.0f,
  //    0.0f)); this->turn(angle);
  //    this->tilt(angle);
  //    this->roll(angle);

  static float f = 0.1;
  if (f <= 1.0) {
    this->updateHeight(f * _f);
    f += 0.01;
  } else {
    this->roll(angle);
  }
}
