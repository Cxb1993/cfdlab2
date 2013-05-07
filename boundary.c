#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){

 int x,y,z,i,idx,fIdx = 0;
 double dens = 0;
 double ci = 0;
 double uwall = 0;

 /* loop over all cells */
 for(z = 0; z < xlength + 2; ++z)
 for(y = 0; y < xlength + 2; ++y)
 for(x = 0; x < xlength + 2; ++x)
 {
	 /* calculate index */
	 idx = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x );

	 /* only border */
	 if (flagField[idx] != FLUID)
	 {
		 for (i = 0; i < Q; ++i)
    	 {
	         /* is the neighbor is in the domain */
	         if ( (z+LATTICEVELOCITIES[i][2] > 0 && z+LATTICEVELOCITIES[i][2] < xlength + 1)
	           && (y+LATTICEVELOCITIES[i][1] > 0 && y+LATTICEVELOCITIES[i][1] < xlength + 1)
	           && (x+LATTICEVELOCITIES[i][0] > 0 && x+LATTICEVELOCITIES[i][0] < xlength + 1) )
	         {
	        	 /* get index of neighbor cell */
	        	 fIdx = ((z+LATTICEVELOCITIES[i][2])*(xlength+2)*(xlength+2)
	                   + (y+LATTICEVELOCITIES[i][1])*(xlength+2)+ (x+LATTICEVELOCITIES[i][0]));

	             uwall = 0;
	             /* calculate second term of boundary treatment */
	             if( flagField[idx] == MOVING_WALL)
	             {
	            	 computeDensity(&collideField[Q*fIdx],&dens);
	                 ci = LATTICEVELOCITIES[i][0]*wallVelocity[0]
	                    + LATTICEVELOCITIES[i][1]*wallVelocity[1]
	                    + LATTICEVELOCITIES[i][2]*wallVelocity[2];

	                 /* formula 18 goes here */
	                 uwall = 2*LATTICEWEIGHTS[i]*dens*ci/(C_S*C_S);
	              }

	              /* formula 16 + 18 goes here */
	              collideField[Q * idx + i] = collideField[Q * fIdx + (Q-1-i)] + uwall;
	          }
    	 }
     }
 }
}

