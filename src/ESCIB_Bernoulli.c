#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io.h"
#include "countPoints.h"
#include "clusters.h"

int main(int argc, char ** argv) {

	if(argc != 9) {
		printf("ERROR! Incorrect number of input arguments\n");
		printf("ESCIB_Bernoulli inputCase inputControl output searchRadius significance(alpha) baselineRatio minCorPointsInEachCluster nonCorePoints\n");
		return 1;
	}

	double xMin = 999999999, yMin = 999999999, xMax = -999999999, yMax = -999999999;

	FILE * inputCas;
	FILE * inputCon;
	FILE * output;

	double radius = atof(argv[4]);
	double significance = atof(argv[5]);

	double baseLineRatio = atof(argv[6]);
	double minCore = atof(argv[7]);
	bool nonCorePoints = true;
	if(atoi(argv[8]) == 0)
		nonCorePoints = false;

	if(NULL == (inputCas = fopen(argv[1], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}
	if(NULL == (inputCon = fopen(argv[2], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}

		
	int countCas = getCount(inputCas, xMin, xMax, yMin, yMax);
	int countCon = getCount(inputCon, xMin, xMax, yMin, yMax);

	double * xCas;
	double * yCas;
	double * xCon;
	double * yCon;

	if(NULL == (xCas = (double *)malloc(sizeof(double) * countCas)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (yCas = (double *)malloc(sizeof(double) * countCas)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (xCon = (double *)malloc(sizeof(double) * countCon)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (yCon = (double *)malloc(sizeof(double) * countCon)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	printf("Number of cases: %d\n", countCas);
	printf("Number of controls: %d\n", countCon);
	printf("X Range: %lf - %lf\n", xMin, xMax);
	printf("Y Range: %lf - %lf\n", yMin, yMax);

	readPoints(inputCas, xCas, yCas);
	readPoints(inputCon, xCon, yCon);

	int nBlockX = ceil((xMax - xMin) / radius);
	int nBlockY = ceil((yMax - yMin) / radius);

//	printf("Index blocks: %d * %d\n", nBlockX, nBlockY);

	int * indexCas;
	int * indexCon;


	indexCas = indexPoints(xCas, yCas, countCas, xMin, yMin, nBlockX, nBlockY, radius);
	indexCon = indexPoints(xCon, yCon, countCon, xMin, yMin, nBlockX, nBlockY, radius);

	fclose(inputCas);
	fclose(inputCon);

	int * countPointsCas = countInDistance_Single(xCas, yCas, indexCas, nBlockX, nBlockY, radius);
	int * countPointsCon = countInDistance_Double(xCas, yCas, xCon, yCon, indexCas, indexCon, nBlockX, nBlockY, radius);

	double p = baseLineRatio * countCas / (countCas + countCon); 

	int * clusters = doClusterBer(xCas, yCas, indexCas, xCon, yCon, indexCon, nBlockX, nBlockY, radius, xMin, yMin, countPointsCas, countPointsCon, p, significance, minCore, nonCorePoints);
	//Output 
	if(NULL == (output = fopen(argv[3], "w"))) {
		printf("ERROR: Can't open the output file.\n");
		exit(1);
	}

	for(int i = 0; i < countCas; i++) {
		fprintf(output, "%lf,%lf,1,%d\n", xCas[i], yCas[i], clusters[i]);
	}

	if(nonCorePoints) {
		for(int i = 0; i < countCon; i++)
		{
			fprintf(output, "%lf,%lf,0,%d\n", xCon[i], yCon[i], clusters[countCas + i]);
		}
	}

	fclose(output);


	free(xCas);
	free(yCas);
	free(indexCas);
	free(countPointsCas);

	free(xCon);
	free(yCon);
	free(indexCon);
	free(countPointsCon);


	free(clusters);

	return 0;
}
