#include <stdio.h>
#include <stdlib.h>

/**
 * NAME:	getCount
 * DESCRIPTION:	get the number of points in a input file; update the bounding box of all points accordingly
 * PARAMETERS:
 * 	FILE * file: the input file 
 * 	double &xMin: the Mininum X of all points, can be updated in this function if necessary
 * 	double &xMax: the Maximum X of all points, can be updated in this function if necessary
 * 	double &yMin: the Minimum Y of all points, can be updated in this function if necessary
 * 	double &xMin: the Maxinum Y of all points, can be updated in this function if necessary
 * RETURN:
 * 	TYPE:	int
 * 	VALUE:	the number of points in the file
 */
int getCount(FILE * file, double &xMin, double &xMax, double &yMin, double &yMax)
{
	int count = 0;
	double x, y;
	int t;
	rewind(file);
	
	while(fscanf(file, "%lf,%lf,%d\n", &x, &y, &t) != EOF) {
		count ++;
		if(x < xMin)
			xMin = x;
		if(x > xMax)
			xMax = x;
		if(y < yMin)
			yMin = y;
		if(y > yMax)
			yMax = y;
	}

	return count;
	
}

/**
 * NAME:	readPoints
 * DESCRIPTION:	read all points (X, Y) in a file
 * PARAMETERS:
 * 	FILE * file: the input file
 * 	double * x: the array to store points' X values
 * 	double * y: the array to store points' Y values
 * RETURN: none
 */
void readPoints(FILE * file, double * x, double * y)
{
	int t;
        rewind(file);
	int count = 0;

	while(fscanf(file, "%lf,%lf,%d\n", x + count, y + count, &t) != EOF) {
		count ++;
	}
}

/**
 * NAME:	indexPoints
 * DESCRIPTION:	index all points based on the block they falls in. the points will be re-ordered based on their blocksIDs accendingly. a seperate index table is created to store the ending array index (in the re-ordered array x and y) of points in each block.
 * PARAMETERS:
 * 	double * &x: 		array points' X values, will be changed to a new array of ordered points
 * 	double * &y: 		array points' Y values, will be changed to a new array of ordered points
 * 	int:			the total number of points
 * 	double xMin:		the minimum X of all points, used to calculate the blockID of each point
 * 	double yMin:		the minimum Y of all points, used to calculate the blockID of each point
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 * 	double blockSize:	the size (side length) of each index block
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	an array with a length equal to (the total number of index blocks + 1), storing the starting and ending array index of points in each block
 */
int * indexPoints(double * &x, double * &y, int count, double xMin, double yMin, int nBlockX, int nBlockY, double blockSize)
{
	int * index;
	int * pointsInB;

	double * newX;
	double * newY;
	
	if(NULL == (index = (int *)malloc(sizeof(int) * (nBlockY * nBlockX + 1))))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (pointsInB = (int *)malloc(sizeof(int) * nBlockY * nBlockX)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (newX = (double *)malloc(sizeof(double) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (newY = (double *)malloc(sizeof(double) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	//Read all points the 1st time to get the number of points in each block
	
	for(int i = 0; i < nBlockY * nBlockX; i++)
	{
		pointsInB[i] = 0;
	}

	int rowID, colID;
	int blockID;
	for(int i = 0; i < count; i++)
	{
		colID = (int)((x[i] - xMin) / blockSize);
		rowID = (int)((y[i] - yMin) / blockSize);
		blockID = colID + rowID * nBlockX;

		pointsInB[blockID] ++;
	}

	index[0] = 0;
	for(int i = 1; i < nBlockY * nBlockX + 1; i++)
	{
		index[i] = index[i - 1] + pointsInB[i - 1];
	}

	//Read all points the 2nd time to fill these points in the new array

	//From this time, pointsInB is used to store the index of next-to-fill points in each block
	pointsInB[0] = 0;
	for(int i = 1; i < nBlockY * nBlockX; i++)
	{
		pointsInB[i] = index[i];
	}

	for(int i = 0; i < count; i++)
	{
		colID = (int)((x[i] - xMin) / blockSize);
		rowID = (int)((y[i] - yMin) / blockSize);
		blockID = colID + rowID * nBlockX;
		newX[pointsInB[blockID]] = x[i];
		newY[pointsInB[blockID]] = y[i];
		pointsInB[blockID] ++;
	}
	

	free(pointsInB);
	free(x);
	free(y);


	x = newX;
	y = newY;


	return index;
}
