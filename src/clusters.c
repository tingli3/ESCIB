#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * NAME:	PossionTest
 * DESCRIPTION:	calculate the probability to get a value equal or larger than nP under a Poisson (lambda) distribution
 * PARAMETERS:
 * 	int nP:	the value from Poisson distribution
 * 	double lambda: the mean of Poisson distribution
 * RETURN:
 * 	TYPE:	double
 * 	VALUE:	the probability to get a value equal or larger than nP
 */
double PossionTest(int nP, double lambda)
{
	double sum = 1.0;
	double element = 1;
	for(int i = 1; i < nP; i++)
	{
		element = element * lambda / i;
		sum += element;
	}

	sum *= exp(-lambda);
	return 1 - sum;
}

/**
 * NAME:	BinomialTest
 * DESCRIPTION:	calculate the probability to get equal or more cases than nCas under a Binomial (nCas, (nCas+nCon), p) distribution
 * PARAMETERS:
 * 	int nCas: the number of cases
 * 	int nCon: the number of controls
 * 	double p: the p of Binomial distribution (e.g., the probability of any point to be a case)
 * RETURN:
 * 	TYPE:	double
 * 	VALUE:	the probability to get a value equal or larger than nCas
 */
double BinomialTest(int nCas, int nCon, double p)
{
	double q = 1 - p;
	int n = nCas + nCon;
	double logElement = n * log(q);
	 
	double sum = exp(logElement);

	for(int i = 1; i < nCas; i++)
	{
		logElement = logElement + log(n+1-i) + log(p) - log(i) - log(q);
		sum += exp(logElement);
	}
	return 1 - sum;
}

/**
 * NAME:	doClusterPoi
 * DESCRIPTION:	cluster all event points based on a Possion Test
 * PARAMETERS:
 * 	double * x: 		the array of event points' X values
 * 	double * y: 		the array of event points' Y values
 * 	int * index:		the index of all event points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 *	double radius:		the search radius, which is also the block size
 *	double xMin:		the minimum X of all points
 *	double yMin:		the minimum Y of all points
 *	int * eC:		the number of events points (within radius) near each event points
 *	double * lambda:	the local lambda of Possion distribution of each event points
 *	double significance: 	the significane level to tell a cluste core point
 *	int minCore:		the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:	whether a cluster include non-core points
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	the cluster ID of each event point
 */
int * doClusterPoi(double * x, double * y, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int * eC, double * lambda, double significance, int minCore, bool nonCorePoints)
{
	int count = index[nBlockX * nBlockY];

	int * clusterID;
	if(NULL == (clusterID = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++)
	{
		if(PossionTest(eC[i], lambda[i]) < significance)
			clusterID[i] = 0;
		else
			clusterID[i] = -1;
	}

	int * pointsToDo;
	if(NULL == (pointsToDo = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int nPToDo = 0;
	int cID = 0;

	double dist2 = radius * radius;

	double cX, cY;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	int iNb;
	int iNbEnd;

	int coreCount;

	for(int i = 0; i < count; i++)
	{
		if(clusterID[i] != 0)
			continue;
		pointsToDo[0] = i;
		nPToDo = 1;
		cID ++;
		clusterID[i] = cID;
		
		coreCount = 1;	

		while(nPToDo > 0) {
			nPToDo --;
			cX = x[pointsToDo[nPToDo]];		
			cY = y[pointsToDo[nPToDo]];

			colID = (int)((cX - xMin) / radius);
			rowID = (int)((cY - yMin) / radius);

			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);

			for(int row = rowMin; row <= rowMax; row ++)
			{
				for(iNb = index[row * nBlockX + colMin]; iNb < index[row * nBlockX + colMax + 1]; iNb ++)
				{
					if(clusterID[iNb] < 1)
					{
						if(dist2 >= ((x[iNb] - cX) * (x[iNb] - cX) + (y[iNb] - cY) * (y[iNb] - cY)))
						{
							if(clusterID[iNb] != -1)
							{
								pointsToDo[nPToDo] = iNb;
								nPToDo ++;
								coreCount ++;
								clusterID[iNb] = cID;
							}
							else if(nonCorePoints)
								clusterID[iNb] = cID;
						}
					}
				}
			}
		
		}

		if(coreCount <= minCore)
		{
			for(int j = 0; j < count; j++)
			{
				if(clusterID[j] == cID)
					clusterID[j] = -1;
			}
			cID --;
		}
		
	}

	free(pointsToDo);
	return clusterID; 
}


/**
 * NAME:	doClusterBer
 * DESCRIPTION:	cluster all event points based on a Binomial Test
 * PARAMETERS:
 * 	double * xCas: 		the array of case points' X values
 * 	double * yCas: 		the array of case points' Y values
 * 	int * indexCas:		the index of all case points
 * 	double * xCon: 		the array of control points' X values
 * 	double * yCon: 		the array of control points' Y values
 * 	int * indexCon:		the index of all control points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 *	double radius:		the search radius, which is also the block size
 *	double xMin:		the minimum X of all points
 *	double yMin:		the minimum Y of all points
 *	int * casC:		the number of case points (within radius) near each case points
 *	int * conC:		the number of control points (within radius) near each case points
 *	double p:		the p of Possion distribution
 *	double significance: 	the significane level to tell a cluste core point
 *	int minCore:		the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:	whether a cluster include non-core points
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	an array of length (countCas + countCon): the cluster ID of each case and control point
 */
int * doClusterBer(double * xCas, double * yCas, int * indexCas, double * xCon, double * yCon, int * indexCon, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int * casC, int * conC, double p, double significance, int minCore, bool nonCorePoints)
{
	int countCas = indexCas[nBlockX * nBlockY];
	int countCon = indexCon[nBlockX * nBlockY];

	int * clusterID;
	if(NULL == (clusterID = (int *)malloc(sizeof(int) * (countCas + countCon))))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < countCas; i++)
	{
		if(BinomialTest(casC[i], conC[i], p) < significance)
			clusterID[i] = 0;
		else
			clusterID[i] = -1;
	}

	int * pointsToDo;
	if(NULL == (pointsToDo = (int *)malloc(sizeof(int) * countCas)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int nPToDo = 0;
	int cID = 0;

	double dist2 = radius * radius;

	double cX, cY;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	int iNb;

	int coreCount;

	for(int i = 0; i < countCas; i++)
	{
		if(clusterID[i] != 0)
			continue;
		pointsToDo[0] = i;
		nPToDo = 1;
		cID ++;
		clusterID[i] = cID;
		
		coreCount = 1;	

		while(nPToDo > 0) {
			nPToDo --;
			cX = xCas[pointsToDo[nPToDo]];		
			cY = yCas[pointsToDo[nPToDo]];

			colID = (int)((cX - xMin) / radius);
			rowID = (int)((cY - yMin) / radius);

			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);

			for(int row = rowMin; row <= rowMax; row ++)
			{
				for(iNb = indexCas[row * nBlockX + colMin]; iNb < indexCas[row * nBlockX + colMax + 1]; iNb ++)
				{
					if(clusterID[iNb] < 1)
					{
						if(dist2 >= ((xCas[iNb] - cX) * (xCas[iNb] - cX) + (yCas[iNb] - cY) * (yCas[iNb] - cY)))
						{
							if(clusterID[iNb] != -1)
							{
								pointsToDo[nPToDo] = iNb;
								nPToDo ++;
								coreCount ++;
								clusterID[iNb] = cID;
							}
							else if(nonCorePoints)
								clusterID[iNb] = cID;
						}
					}
				}

				if(nonCorePoints) {
					for(iNb = indexCon[row * nBlockX + colMin]; iNb < indexCon[row * nBlockX + colMax + 1]; iNb ++)
					{
						if(clusterID[countCas + iNb] < 1)
						{
							if(dist2 >= ((xCon[iNb] - cX) * (xCon[iNb] - cX) + (yCon[iNb] - cY) * (yCon[iNb] - cY)))
							{
								clusterID[countCas + iNb] = cID;
							}
						}
					}
				}
			}
		
		}

		if(coreCount <= minCore)
		{
			for(int j = 0; j < (countCas + countCon); j++)
			{
				if(clusterID[j] == cID)
					clusterID[j] = -1;
			}
			cID --;
		}
		
	}

	for(int j = countCas; j < (countCas + countCon); j++)
	{
		if(clusterID[j] == 0)
			clusterID[j] = -1;
	}

	free(pointsToDo);
	return clusterID; 
}


/**
 * NAME:	doClusterDBSCAN
 * DESCRIPTION:	cluster all event points using DBSCAN algorithm
 * PARAMETERS:
 * 	double * x: 		the array of events' X values
 * 	double * y: 		the array of events' Y values
 * 	int * index:		the index of all event points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 *	double radius:		the search radius, which is also the block size
 *	int minPts:		the minimum points to form a core points
 *	double xMin:		the minimum X of all points
 *	double yMin:		the minimum Y of all points
 *	int * eC:		the number of event points (within radius) near each event points
 *	int minCore:		the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:	whether a cluster include non-core points
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	an array of length count: the cluster ID of each case and control point
 */
int * doClusterDBSCAN(double * x, double * y, int * index, int nBlockX, int nBlockY, double radius, int minPts, double xMin, double yMin, int * eC, int minCore, bool nonCorePoints) {

	int count = index[nBlockX * nBlockY];

	int * clusterID;
	if(NULL == (clusterID = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++)
	{
		if(eC[i] >= minPts)
			clusterID[i] = 0;
		else
			clusterID[i] = -1;
	}

	int * pointsToDo;
	if(NULL == (pointsToDo = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int nPToDo = 0;
	int cID = 0;

	double dist2 = radius * radius;

	double cX, cY;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	int iNb;

	int coreCount;

	for(int i = 0; i < count; i++)
	{
		if(clusterID[i] != 0)
			continue;
		pointsToDo[0] = i;
		nPToDo = 1;
		cID ++;
		clusterID[i] = cID;
		
		coreCount = 1;	

		while(nPToDo > 0) {
			nPToDo --;
			cX = x[pointsToDo[nPToDo]];		
			cY = y[pointsToDo[nPToDo]];

			colID = (int)((cX - xMin) / radius);
			rowID = (int)((cY - yMin) / radius);

			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);

			for(int row = rowMin; row <= rowMax; row ++)
			{
				for(iNb = index[row * nBlockX + colMin]; iNb < index[row * nBlockX + colMax + 1]; iNb ++)
				{
					if(clusterID[iNb] < 1)
					{
						if(dist2 >= ((x[iNb] - cX) * (x[iNb] - cX) + (y[iNb] - cY) * (y[iNb] - cY)))
						{
							if(clusterID[iNb] != -1)
							{
								pointsToDo[nPToDo] = iNb;
								nPToDo ++;
								coreCount ++;
								clusterID[iNb] = cID;
							}
							else if(nonCorePoints)
								clusterID[iNb] = cID;
						}
					}
				}
			}
		
		}

		if(coreCount <= minCore)
		{
			for(int j = 0; j < count; j++)
			{
				if(clusterID[j] == cID)
					clusterID[j] = -1;
			}
			cID --;
		}
		
	}


	free(pointsToDo);
	return clusterID;
}
