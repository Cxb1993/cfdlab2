#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double * const currentCell, double *density)
{
	/* See formula 9 */
	double sum;
	int i;
	sum = 0;
	for(i=0; i<Q; i++)
	{
		sum += currentCell[i];
	}
	*density = sum;
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity)
{
	/* Compute each dimension separated, see figure 1 and table 1 */
	int i;
	velocity[0] = 0;
	velocity[1] = 0;
	velocity[2] = 0;
	for(i = 0; i<Q; i++)
	{
		velocity[0] += currentCell[i]*LATTICEVELOCITIES[i][0];
		velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1];
		velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2];
	}


	velocity[0] = velocity[0] / *density;
	velocity[1] = velocity[1] / *density;
	velocity[2] = velocity[2] / *density;
}

void computeFeq(const double * const density, const double * const velocity, double *feq)
{
	/* See formula 10 */
	int i;
	for(i = 0; i < Q; i++)
	{
		feq[i] = LATTICEWEIGHTS[i] * *density * (1
				+ (LATTICEVELOCITIES[i][0]*velocity[0] +
					LATTICEVELOCITIES[i][1]*velocity[1] +
					LATTICEVELOCITIES[i][2]*velocity[2])/(C_S*C_S)
				+ (LATTICEVELOCITIES[i][0]*velocity[0] +
						LATTICEVELOCITIES[i][1]*velocity[1] +
						LATTICEVELOCITIES[i][2]*velocity[2])*
						(LATTICEVELOCITIES[i][0]*velocity[0] +
						LATTICEVELOCITIES[i][1]*velocity[1] +
						LATTICEVELOCITIES[i][2]*velocity[2])/(2*C_S*C_S*C_S*C_S)
				- (velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*C_S*C_S)
				);
	}
}
