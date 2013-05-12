#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

  /* read velocity values */
  double vw1,vw2,vw3;
  read_double( argv[1], "vw1", &vw1 );
  read_double( argv[1], "vw2", &vw2 );
  read_double( argv[1], "vw3", &vw3 );
  velocityWall[0] = vw1;
  velocityWall[1] = vw2;
  velocityWall[2] = vw3;

  READ_INT   ( argv[1], *xlength);
  READ_INT   ( argv[1], *timesteps );
  READ_INT   ( argv[1], *timestepsPerPlotting );
  READ_DOUBLE( argv[1], *tau );
  return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){

  int x,y,z,i, idx = 0;

  for(z = 0; z < xlength + 2; ++z)
  for(y = 0; y < xlength + 2; ++y)
  for(x = 0; x < xlength + 2; ++x)
  {
    /* index in the field */
    idx = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x );

    /* NO_SLIP at  boundary, MOVING_WALL on top, else FLUID */
    if (z == xlength + 1)
    	flagField[idx] = MOVING_WALL;
    else if (x == 0 || x == xlength + 1 || y == 0 || y == xlength + 1 || z == 0)
    	flagField[idx] = NO_SLIP;
    else
    	flagField[idx] = FLUID;

    for(i = 0; i < Q; ++i)
    {
    	collideField[Q * idx + i] = LATTICEWEIGHTS[i];
    	streamField [Q * idx + i] = LATTICEWEIGHTS[i];
    }
  }
}

