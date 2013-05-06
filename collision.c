#include "collision.h"
#include "computeCellValues.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq)
{
	/* See formula 13 */
	int i;
	for(i=0; i<19; i++)
	{
		currentCell[i] = currentCell[i] - 1.0/ *tau * (currentCell[i] - feq[i]);
	}
}

void doCollision(double *collideField, int *flagField, const double * const tau,int xlength){
  double density;
  double velocity[3];
  double feq[19];
  double* currentCell;
  int currentIndex;

  /* Set pointer to first cell that is not a boundary */
  currentIndex = xlength+3;

  /* Last row is a border */
  while(currentIndex < 19*(xlength+2)*(xlength+2)*(xlength+2) - 19*(xlength+2))
  {
	  /* No left and right border, only inner cells */
	  if(currentIndex % (xlength+2) != 0 && currentIndex % (xlength+2) != xlength+1)
	  {
			currentCell = &collideField[currentIndex];
			computeDensity(currentCell, &density);
			computeVelocity(currentCell, &density, velocity);
			computeFeq(&density, velocity, feq);
			computePostCollisionDistributions(currentCell, tau, feq);
			/* Next cell begins after 19 f-values */
			currentIndex += 19;
	  }
  }
}

