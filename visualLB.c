#include "visualLB.h"
#include "computeCellValues.h"
#include "helper.h"
#include "LBDefinitions.h"


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
	fprintf(fp,"generated by CFD-lab course output (originally written by Tobias Neckel, modified by Philipp Huebner and Leonard Hoess)\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"\n");
	fprintf(fp,"DATASET STRUCTURED_GRID\n");
	/*
	fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength+1, xlength+1, xlength+1);
	fprintf(fp,"POINTS %i float\n", (xlength+1)*(xlength+1)*(xlength+1) );
	*/
	fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength, xlength, xlength);
	fprintf(fp,"POINTS %i float\n", (xlength)*(xlength)*(xlength) );
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

	for(i = 0; i < xlength+1; i++)
	{
		for(j = 0; j < xlength+1; j++)
		{
			for(k = 0; k < xlength+1; k++)
			{

				fprintf(fp, "%f %f %f\n", originX+i, originY+j, originZ + k);
			}
		}
	}
}


void writeVtkOutput(const double * const collideField, const int * const flagField, const char* filename, unsigned int t,int xlength)
{
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

	fprintf(fp,"\n");
	fprintf(fp,"POINT_DATA %i \n", (xlength)*(xlength)*(xlength) );
	fprintf(fp, "VECTORS velocity float\n");
	/* Output only inner cells from 1 to xlength + 1 */
	for(x = 1; x < xlength+1; x++)
	{
		for(y = 1; y < xlength+1; y++)
		{
			for(z = 1; z < xlength+1; z++)
			{
				/* Compute adress of current cell */
				computeDensity(&collideField[Q*(z*(xlength+2)*(xlength+2) + y*(xlength+2) + x)], &density);
				computeVelocity(&collideField[Q*(z*(xlength+2)*(xlength+2) + y*(xlength+2) + x)], &density, velocity);
				fprintf(fp, "%f %f %f\n", velocity[0], velocity[1], velocity[2]);
			}
		}
	}

	fprintf(fp,"\n");
	fprintf(fp, "SCALARS pressure float 1 \n");
	fprintf(fp, "LOOKUP_TABLE default \n");
	for(x = 1;  x < xlength+1; x++)
	{
		for(y = 1; y < xlength+1; y++)
		{
			for(z = 1; z < xlength+1; z++)
			{
				computeDensity(&collideField[Q*(z*(xlength+2)*(xlength+2) + y*(xlength+2) + x)], &density);
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

