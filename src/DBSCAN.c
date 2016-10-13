#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.h"
#include "countPoints.h"
#include "clusters.h"

int main(int argc, char ** argv) {
	
	if(argc != 7) {
		printf("ERROR! Incorrect number of input arguments\n");
		printf("DBSCAN inputEvents output searchRadius minPts minCorPointsInEachCluster nonCorePoints\n");
		return 1;
	}

	double xMin = 999999999, yMin = 999999999, xMax = -999999999, yMax = -999999999;
	
	FILE * input;
	FILE * output;

	double radius = atof(argv[3]);
	int minPts = atoi(argv[4]);
	double minCore = atof(argv[5]);
	bool nonCorePoints = true;
	if(atoi(argv[6]) == 0)
		nonCorePoints = false;
		

	if(NULL == (input = fopen(argv[1], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}

	int count = getCount(input, xMin, xMax, yMin, yMax);

	double * x;
	double * y;
	
	if(NULL == (x = (double *)malloc(sizeof(double) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (y = (double *)malloc(sizeof(double) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	readPoints(input, x, y);

	fclose(input);
	
	int nBlockX = ceil((xMax - xMin) / radius);
	int nBlockY = ceil((yMax - yMin) / radius);

	int * index;

	index = indexPoints(x, y, count, xMin, yMin, nBlockX, nBlockY, radius);
	
	int * countPoints = countInDistance_Single(x, y, index, nBlockX, nBlockY, radius);

	int * clusters = doClusterDBSCAN(x, y, index, nBlockX, nBlockY, radius, minPts, xMin, yMin, countPoints, minCore, nonCorePoints);
	
	//Output 
	if(NULL == (output = fopen(argv[2], "w"))) {
		printf("ERROR: Can't open the output file.\n");
		exit(1);
	}
	
	for(int i = 0; i < count; i++)
	{
		fprintf(output, "%lf,%lf,%d\n", x[i], y[i], clusters[i]);
	}

	fclose(output);

	free(clusters);



	free(x);
	free(y);

	free(index);
	free(countPoints);

	return 0;
}
