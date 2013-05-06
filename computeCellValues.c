#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <math.h>

void computeDensity(const double const * currentCell, double *density)
{
	/* See formula 9 */
	double sum;
	int i;
	for(i=0; i<19; i++)
	{
		sum += currentCell[i];
	}
	density = sum;
}

void computeVelocity(const double const * currentCell, const double * const density, double *velocity)
{
	/* Compute each dimension separated, see figure 1 and table 1 */
	velocity[0] = -currentCell[1] + currentCell[3] - currentCell[5] +
			currentCell[7] - currentCell[8] + currentCell[10] - currentCell[11] +
			currentCell[13] - currentCell[15] + currentCell[17];
	velocity[1] = -currentCell[0] + currentCell[4] - currentCell[5] -
			currentCell[6] - currentCell[7] + currentCell[11] + currentCell[12] +
			currentCell[13] - currentCell[14] + currentCell[18];
	velocity[2] = -currentCell[0] - currentCell[1] - currentCell[2] -
			currentCell[3] - currentCell[4] + currentCell[14] + currentCell[15] +
			currentCell[16] + currentCell[17] + currentCell[18];

	velocity[0] /= density;
	velocity[1] /= density;
	velocity[2] /= density;
}

void computeFeq(const double const * density, const double const * velocity, double *feq)
{
	/* See formula 10 */
	int i;
	for(i = 0; i < 19; i++)
	{
		feq[i] = LATTICEWEIGHTS[i] * density * (1
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

/* OLD CODE
	double cs;
	cs = 1 / 1.732050808;
	feq[0] = 1/36 * density * (1 + (-velocity[1]-velocity[2])/(cs*cs) +
			(-velocity[1]-velocity[2])*(-velocity[1]-velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[1] = 1/36 * density * (1 + (-velocity[0]-velocity[2])/(cs*cs) +
			(-velocity[0]-velocity[2])*(-velocity[0]-velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[2] = 2/36 * density * (1 + (-velocity[2])/(cs*cs) +
			(-velocity[2])*(-velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[3] = 1/36 * density * (1 + (velocity[0]-velocity[2])/(cs*cs) +
			(velocity[0]-velocity[2])*(velocity[0]-velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[4] = 1/36 * density * (1 + (velocity[1]-velocity[2])/(cs*cs) +
			(velocity[1]-velocity[2])*(velocity[1]-velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[5] = 1/36 * density * (1 + (-velocity[0]-velocity[1])/(cs*cs) +
			(-velocity[0]-velocity[1])*(-velocity[0]-velocity[1])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[6] = 2/36 * density * (1 + (-velocity[1])/(cs*cs) +
			(-velocity[1])*(-velocity[1])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[7] = 1/36 * density * (1 + (velocity[0]-velocity[1])/(cs*cs) +
			(velocity[0]-velocity[1])*(velocity[0]-velocity[1])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[8] = 2/36 * density * (1 + (-velocity[0])/(cs*cs) +
			(-velocity[0])*(-velocity[0])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[9] = 12/36 * density * (1 -(velocity[0]*velocity[0] + velocity[1]*velocity[1] +
			velocity[2]*velocity[2])/(2*cs*cs));
	feq[10] = 2/36 * density * (1 + (velocity[0])/(cs*cs) +
			(velocity[0])*(velocity[0])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[11] = 1/36 * density * (1 + (-velocity[0]+velocity[1])/(cs*cs) +
			(-velocity[0]+velocity[1])*(-velocity[0]+velocity[1])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[12] = 2/36 * density * (1 + (velocity[1])/(cs*cs) +
			(velocity[1])*(velocity[1])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[13] = 1/36 * density * (1 + (velocity[0]+velocity[1])/(cs*cs) +
			(velocity[0]+velocity[1])*(velocity[0]+velocity[1])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[14] = 1/36 * density * (1 + (-velocity[1]+velocity[2])/(cs*cs) +
			(-velocity[1]+velocity[2])*(-velocity[1]+velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[15] = 1/36 * density * (1 + (-velocity[0]+velocity[2])/(cs*cs) +
			(-velocity[0]+velocity[2])*(-velocity[0]+velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[16] = 2/36 * density * (1 + (velocity[2])/(cs*cs) +
			(velocity[2])*(velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[17] = 1/36 * density * (1 + (velocity[0]+velocity[2])/(cs*cs) +
			(velocity[0]+velocity[2])*(velocity[0]+velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));
	feq[18] = 1/36 * density * (1 + (velocity[1]+velocity[2])/(cs*cs) +
			(velocity[1]+velocity[2])*(velocity[1]+velocity[2])/(2*cs*cs*cs*cs) -
			(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])/(2*cs*cs));*/
}
