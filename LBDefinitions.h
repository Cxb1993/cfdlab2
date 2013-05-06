#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q           19
#define FLUID       0
#define NO_SLIP     1
#define MOVING_WALL 2

  static const int LATTICEVELOCITIES[19][3] = {
		  {0,-1,-1}, {-1,0,-1}, {0,0,-1}, {1,0,-1}, {0,1,-1},
  	  	  {-1,-1,0}, {0,-1,0}, {1,-1,0}, {-1,0,0}, {0,0,0},
  	  	  {1,0,0}, {-1,1,0}, {0,1,0}, {1,1,0}, {0,-1,1},
  	  	  {-1,0,1}, {0,0,1}, {1,0,1}, {0,1,1}};
  static const double LATTICEWEIGHTS[19] = {
		  1/36, 1/36, 2/36, 1/36, 1/36,
  	  	  1/36, 2/36, 1/36, 2/36, 12/36,
  	  	  2/36, 1/36, 2/36, 1/36, 1/36,
  	  	  1/36, 2/36, 1/36, 1/36};
  static const double C_S = 1.0/1.732050808;

#endif

