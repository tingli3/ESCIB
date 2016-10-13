#include <stdio.h>
#include <stdlib.h>

/**
 * NAME:	countInDistance_Single
 * DESCRIPTION:	get the number of type A points within a distance of each type A point
 * PARAMETERS:
 * 	double * xE:		type A points' X values 
 * 	double * yE:		type A points' Y values 
 * 	int * indexE:		the index of type A points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 * 	double distance:	the distance, which is also the size (side length) of each index block
 * RETURN:
 * 	TYPE:	int * 
 * 	VALUE:	an array of the numbers of points within the distance, ordered the same as xE and yE
 */

int * countInDistance_Single(double * xE, double * yE, int * indexE, int nBlockX, int nBlockY, double distance)
{
	int countE = indexE[nBlockX * nBlockY];
	int * count;
	
	if(NULL == (count = (int *)malloc(sizeof(int) * countE)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	double x, y;
	double dis2 = distance * distance;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;
	int iC, iP;
	int pCol, pRow;
	for(rowID = 0; rowID < nBlockY; rowID ++)
	{
		for(colID = 0; colID < nBlockX; colID ++)
		{
			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);
			for(iC = indexE[rowID * nBlockX + colID]; iC < indexE[rowID * nBlockX + colID + 1]; iC++)
			{
				x = xE[iC];
				y = yE[iC];
				count[iC] = 0;
				for(int row = rowMin; row <= rowMax; row ++)
				{
					for(iP = indexE[row * nBlockX + colMin]; iP < indexE[row * nBlockX + colMax + 1]; iP ++)
					{
						if(dis2 >= ((xE[iP] - x) * (xE[iP] - x) + (yE[iP] - y) * (yE[iP] - y)))
							count[iC] ++;
					}

				}

			}
		}
	}
	return count;
}

/**
 * NAME:	countInDistance_Double
 * DESCRIPTION:	get the number of type B points within a distance of each type A point
 * PARAMETERS:
 * 	double * xE:		type A points' X values 
 * 	double * yE:		type A points' Y values 
 * 	double * xB:		type B points' X values 
 * 	double * yB:		type B points' Y values 
 * 	int * indexE:		the index of type A points
 * 	int * indexB:		the index of type B points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 * 	double distance:	the distance, which is also the size (side length) of each index block
 * RETURN:
 * 	TYPE:	int * 
 * 	VALUE:	an array of the numbers of points within the distance
 */

int * countInDistance_Double(double * xE, double * yE, double * xB, double * yB, int * indexE, int * indexB, int nBlockX, int nBlockY, double distance)
{
	int countE = indexE[nBlockX * nBlockY];
	int countB = indexB[nBlockX * nBlockY];

	int * count;
	
	if(NULL == (count = (int *)malloc(sizeof(int) * countE)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	double x, y;
	double dis2 = distance * distance;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;
	int iC, iP;
	int pCol, pRow;
	for(rowID = 0; rowID < nBlockY; rowID ++)
	{
		for(colID = 0; colID < nBlockX; colID ++)
		{
			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);
			for(iC = indexE[rowID * nBlockX + colID]; iC < indexE[rowID * nBlockX + colID + 1]; iC++)
			{
				x = xE[iC];
				y = yE[iC];
				count[iC] = 0;
				for(int row = rowMin; row <= rowMax; row ++)
				{
					for(iP = indexB[row * nBlockX + colMin]; iP < indexB[row * nBlockX + colMax + 1]; iP ++)
					{
						if(dis2 >= ((xB[iP] - x) * (xB[iP] - x) + (yB[iP] - y) * (yB[iP] - y)))
							count[iC] ++;
					}

				}

			}
		}
	}
	return count;
}

