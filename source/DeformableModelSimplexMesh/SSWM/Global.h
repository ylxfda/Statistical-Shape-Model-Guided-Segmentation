#ifndef Gloabal_h  
#define Gloabal_h       1

#include "Piont3D.h"

/********* NEIGHBORS DIRECTIONS ************/

  static const char  NORTH     = 4;
  static const char  SOUTH     = 2;
  static const char  WEST      = 1;
  static const char  EAST      = 3;

  static const char  TL     = 5;
  static const char  TR     = 6;
  static const char  BL      = 7;
  static const char  BR      = 8;

  /********* NEIGHBORS DIRECTIONS ************/

  static const char  _LEFT      = 0;
  static const char  _RIGHT     = 1;
  static const char  _NORTH     = 2;
  static const char  _SOUTH     = 3;
  static const char  _WEST      = 4;
  static const char  _EAST      = 5;
  static const char  _NW        = 6;
  static const char  _NE        = 7;
  static const char  _SE        = 8;
  static const char  _SW        = 9;
  static const char  _NW_OL     = 10;
  static const char  _NE_OL     = 11;
  static const char  _SE_OL     = 12;
  static const char  _SW_OL     = 13;

  /********* Types of Coefficients ************/

  static const char  NorthEdgeCoeffi = 0;
  static const char  SouthEdgeCoeffi = 1;
  static const char  WestEdgeCoeffi = 2;
  static const char  EastEdgeCoeffi = 3;
  static const char  FacetCoeffi = 4;
  static const char  ScaleCoffi =5;

  struct  Vertax
  {
	  int TreeIndex;
	  int NodeIndex;
	  char Position;
  };		

  struct Facet
  {
	  int TreeIndex;
	  int NodeIndex;
  };

  struct Edge
  {
	  int TreeIndex;
	  int NodeIndex;
	  char EdgeID;
  };

  struct WaveCoeffi   //for the wavelet coefficients list
  {
	  int TreeIndex;  //the tree this coefficient belongs to
	  int NodeIndex;  //the parent quadriangle contains this child wavelet coefficient
	  Point3D* Coeffi;  //pointer to where the coefficient stores.
	  char TypeID;      //the type of the coefficient it contains.   
  };

#endif