
#ifndef PointSet_h
#define PointSet_h 1

#include "Piont3D.h"
#include "Global.h"
#include <set>
#include <list>
using namespace std;

typedef set<Point3D*, greater<Point3D*> > PointSet;
typedef list<Point3D*> PointList;
typedef list<struct Vertax> VertexList;
typedef list<struct Edge> EdgeList;
typedef list<struct Facet> FaceList;
typedef list<struct WaveCoeffi> WaveCoeffiList;

#endif