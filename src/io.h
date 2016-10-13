#ifndef IOH
#define IOH

int getCount(FILE * file, double &xMin, double &xMax, double &yMin, double &yMax);
void readPoints(FILE * file, double * x, double * y);
int * indexPoints(double * &x, double * &y, int count, double xMin, double yMin, int nBlockX, int nBlockY, double blockSize);

#endif
