#include "collision.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq)
{
	/* See formula 13 */
	int i;
	for(i=0; i<Q; i++)
	{
		currentCell[i] = currentCell[i] - 1.0/ *tau * (currentCell[i] - feq[i]);
	}
}

void doCollision(double *collideField, int *flagField, const double * const tau,int xlength){
  double density;
  double velocity[3];
  double feq[Q];
  double* currentCell;
  int currentIndex;

  /* Set pointer to first cell that is not a boundary */
  currentIndex = (xlength+3);

  /* Last row is a border */
  while(currentIndex < (xlength+2)*(xlength+2)*(xlength+2) - (xlength+2))
  {
	  /* No left and right border, only inner cells */
	  if(flagField[currentIndex] == FLUID)
	  {
			currentCell = &collideField[currentIndex*Q];
			computeDensity(currentCell, &density);
			computeVelocity(currentCell, &density, velocity);
			computeFeq(&density, velocity, feq);
			computePostCollisionDistributions(currentCell, tau, feq);
	  }
	  currentIndex++;
  }
}

