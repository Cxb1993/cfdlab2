#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"

int main(int argc, char *argv[]){

  double  *collideField = NULL;
  double  *streamField = NULL;
  int      *flagField = NULL;
  double  tau;
  double  velocityWall[3];
  int     timesteps, timestepsPerPlotting, xlength;

  char* filename;
  int t;

  /* check if the argument is passed */
   if(argc < 2)
   {
       printf("pass the configuration file as first argument\n");
       return 0;
   }

  /* read the parameters*/
  readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);

  filename = "./Out/Output";
  t = 0;

  initialiseFields(collideField, streamField, flagField, xlength);

  for(t = 0; t < timesteps; t++)
  {
	  double *swap = NULL;
	  doStreaming(collideField, streamField, flagField, xlength);
	  swap = collideField;
	  collideField = streamField;
	  streamField = swap;

	  doCollision(collideField, flagField, &tau, xlength);
	  treatBoundary(collideField, flagField, velocityWall, xlength);

	  if(t%timestepsPerPlotting == 0)
	  {
		  writeVtkOutput(collideField, flagField, filename, t, xlength);
	  }
  }

  return 0;
}

#endif

