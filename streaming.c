#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){

    int x,y,z,i,idx,fIdx = 0;

    for (z = 0; z < xlength + 2; ++z)
    for (y = 0; y < xlength + 2; ++y)
    for (x = 0; x < xlength + 2; ++x)
    {
        /* index in the field */
        idx = ( z*(xlength+2)*(xlength+2) + y*(xlength+2) + x );
        if (flagField[idx] == FLUID)
        {
        	for (i = 0; i < Q; ++i)
            {
        		/* index in the field */
        		fIdx = ( (    z+LATTICEVELOCITIES[i][2])*(xlength+2)*(xlength+2)
                           + (y+LATTICEVELOCITIES[i][1])*(xlength+2)
                           + (x+LATTICEVELOCITIES[i][0]));
        		/* add the inverted vector*/
        	    streamField[Q * idx + (Q-1-i)] = collideField[Q * fIdx + (Q-1-i)];

           }
        }
    }
}

