#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.h"
#include "countPoints.h"
#include "clusters.h"

int main(int argc, char ** argv) {

	if(argc != 9) {
		printf("ERROR! Incorrect number of input arguments\n");
		printf("ESCIB_Poisson inputBackground inputEvents output searchRadius significance(alpha) baselineRatio minCorPointsInEachCluster nonCorePoints\n");
		return 1;
	}

	double xMin = 999999999, yMin = 999999999, xMax = -999999999, yMax = -999999999;

	FILE * inputB;
	FILE * inputE;
	FILE * output;

	double radius = atof(argv[4]);
	double significance = atof(argv[5]);

	double baseLineRatio = atof(argv[6]);
	double minCore = atof(argv[7]);
	bool nonCorePoints = true;
	if(atoi(argv[8]) == 0)
		nonCorePoints = false;

	if(NULL == (inputB = fopen(argv[1], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}
	if(NULL == (inputE = fopen(argv[2], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}

		
	int countB = getCount(inputB, xMin, xMax, yMin, yMax);
	int countE = getCount(inputE, xMin, xMax, yMin, yMax);

	double * xB;
	double * yB;
	double * xE;
	double * yE;

	if(NULL == (xB = (double *)malloc(sizeof(double) * countB)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (yB = (double *)malloc(sizeof(double) * countB)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (xE = (double *)malloc(sizeof(double) * countE)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (yE = (double *)malloc(sizeof(double) * countE)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	printf("Number of background points: %d\n", countB);
	printf("Number of event points: %d\n", countE);
	printf("X Range: %lf - %lf\n", xMin, xMax);
	printf("Y Range: %lf - %lf\n", yMin, yMax);
	printf("Search radius %lf\n", radius);

	readPoints(inputB, xB, yB);
	readPoints(inputE, xE, yE);

	int nBlockX = ceil((xMax - xMin) / radius);
	int nBlockY = ceil((yMax - yMin) / radius);

//	printf("Index blocks: %d * %d\n", nBlockX, nBlockY);

	int * indexB;
	int * indexE;


	indexB = indexPoints(xB, yB, countB, xMin, yMin, nBlockX, nBlockY, radius);
	indexE = indexPoints(xE, yE, countE, xMin, yMin, nBlockX, nBlockY, radius);

	fclose(inputB);
	fclose(inputE);

	int * countPointsE = countInDistance_Single(xE, yE, indexE, nBlockX, nBlockY, radius);
	int * countPointsB = countInDistance_Double(xE, yE, xB, yB, indexE, indexB, nBlockX, nBlockY, radius);

	free(xB);
	free(yB);
	free(indexB);

	double * lambda;
	if(NULL == (lambda = (double *)malloc(sizeof(double) * countE)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < countE; i++)
	{
		//lambda[i] = (double)(countPointsB[i]) * countE / countB;
		lambda[i] = (double)(countPointsB[i]) * countE * baseLineRatio / countB;
	}
	
	free(countPointsB);


	int * clusters =  doClusterPoi(xE, yE, indexE, nBlockX, nBlockY, radius, xMin, yMin, countPointsE, lambda, significance, minCore, nonCorePoints);
	//Output 
	if(NULL == (output = fopen(argv[3], "w"))) {
		printf("ERROR: Can't open the output file.\n");
		exit(1);
	}


	for(int i = 0; i < countE; i++) {
		fprintf(output, "%lf,%lf,%d\n", xE[i], yE[i], clusters[i]);
	}

	fclose(output);

	free(countPointsE);
	free(xE);
	free(yE);
	free(indexE);
	free(lambda);

	free(clusters);

	return 0;
}
