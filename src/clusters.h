#ifndef CH
#define CH
//Poisson
int * doClusterPoi(double * x, double * y, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int * eC, double * lambda, double significance, int minCores, bool nonCorePoints);
//Bernoulli
int * doClusterBer(double * xCas, double * yCas, int * indexCas, double * xCon, double * yCon, int * indexCon, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int * casC, int * conC, double p, double significance, int minCore, bool nonCorePoints);
//DBSCAN
int * doClusterDBSCAN(double * x, double * y, int * index, int nBlockX, int nBlockY, double radius, int minPts, double xMin, double yMin, int * eC, int minCore, bool nonCorePoints);

#endif
