#include "visualLB.h"
#include "computeCellValues.h"
#include "helper.h"


void write_vtkHeader( FILE *fp, int xlength)
{
	if( fp == NULL )
	{
		char szBuff[80];
		sprintf( szBuff, "Null pointer in write_vtkHeader" );
		ERROR( szBuff );
		return;
	}

	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"generated by CFD-lab course output (originally written by Tobias Neckel, modified by Philipp Huebner and Leonard Hoess) \n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"\n");
	fprintf(fp,"DATASET STRUCTURED_GRID\n");
	fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength+2, xlength+2, xlength+2);
	fprintf(fp,"POINTS %i float\n", (xlength+2)*(xlength+2)*(xlength+2) );
	fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int xlength)
{
	double originX = 0.0;
	double originY = 0.0;
	double originZ = 0.0;

	int i = 0;
	int j = 0;
	int k = 0;

	for(i = 0; i < xlength+2; i++)
	{
		for(j = 0; j < xlength+2; j++)
		{
			for(k = 0; k < xlength+2; k++)
			{
				/* dx = dy = dz = 1 */
				fprintf(fp, "%f %f %f\n", originX+i, originY+j, originZ + k);
			}
		}
	}
}


void writeVtkOutput(const double * const collideField, const int * const flagField, const char* filename, unsigned int t,int xlength)
{
	double* currentCell;
	double density;
	double velocity[3];

	int x, y, z;

	FILE *fp=NULL;
	char szFileName[80];
	sprintf( szFileName, "%s.%i.vtk", filename, t);
	fp = fopen( szFileName, "w");
	if(fp == NULL)
	{
		char szBuff[80];
		sprintf( szBuff, "Failed to open %s", szFileName);
		ERROR( szBuff );
		return;
	}

	write_vtkHeader(fp, xlength);
	write_vtkPointCoordinates(fp, xlength);

	fprintf(fp,"CELL_DATA %i \n", (xlength+2)*(xlength+2)*(xlength+2) );

	fprintf(fp,"\n");
	fprintf(fp, "VECTORS velocity float\n");
	/* Output only inner cells from 1 to xlength - 1 */
	for(x = 1; x < xlength+1; x++)
	{
		for(y = 1; y < xlength+1; y++)
		{
			for(z = 1; z < xlength+1; z++)
			{
				/* Compute adress of current cell */
				currentCell = collideField + 19*(z*(xlength*xlength) + y*xlength + x);
				computeDensity(currentCell, &density);
				computeVelocity(currentCell, &density, velocity);
				fprintf(fp, "%f %f %f\n", velocity[0], velocity[1], velocity[2]);
			}
		}
	}

	fprintf(fp,"\n");
	fprintf(fp,"CELL_DATA %i \n", (xlength+2)*(xlength+2)*(xlength+2) );
	fprintf(fp, "SCALARS pressure float 1 \n");
	fprintf(fp, "LOOKUP_TABLE default \n");
	for(x = 1;  x < xlength+1; x++)
	{
		for(y = 1; y < xlength+1; y++)
		{
			for(z = 1; z < xlength+1; z++)
			{
				currentCell = &collideField[19*(z*(xlength*xlength) + y*xlength + x)];
				computeDensity(currentCell, &density);
				fprintf(fp, "%f\n", density);
			}
		}
	}

	if( fclose(fp) )
	{
		char szBuff[80];
		sprintf( szBuff, "Failed to close %s", szFileName );
		ERROR( szBuff );
	}
}
