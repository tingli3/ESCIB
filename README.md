# ESCIB
**E**xpansion-based **S**patial **C**lustering over **I**nhomogeneous **B**ackground

## ESCIB_Bernoulli
ESCIB with a Bernoulli model, used for case-control study
### To execute:
  ESCIB_Bernoulli inputCase inputControl output searchRadius significance(alpha) baselineRatio minCorPointsInEachCluster nonCorePoints
### Arguments:
1. inputCase: input file of case points, a csv without header with two columns: x and y
2. inputControl: input file of control points, a csv without header with two columns: x and y
3. output: output file name
4. searchRadius: search radius to check significance and to expand clusters
5. significance(alpha): significance level to decide core points
6. baselineRatio: the ratio null hypothesis to complete randomness baseline 1 means the same as baseline, 2 means twice the baseline
7. minCorPointsInEachCluster: minimum number of core points in each cluster
8. nonCorePoints: whether clusters should keep non-core points
  * 0: not keeping
  * 1: keeping
  
## ESCIB_Poisson
ESCIB with a (inhomogeneous Poisson) model, used for detecting spatial clusters over a changing background intensity
### To execute:
  ESCIB_Poisson inputBackground inputEvents output searchRadius significance(alpha) baselineRatio minCorPointsInEachCluster nonCorePoints
1. inputBackground: input file of background points, a csv without header with two columns: x and y
2. inputEvents: input file of event points, a csv without header with two columns: x and y
3. output: output file name
4. searchRadius: search radius to check significance and to expand clusters
5. significance(alpha): significance level to decide core points
6. baselineRatio: the ratio null hypothesis to complete randomness baseline 1 means the same as baseline, 2 means twice the baseline
7. minCorPointsInEachCluster: minimum number of core points in each cluster
8. nonCorePoints: whether clusters should keep non-core points
  * 0: not keeping
  * 1: keeping

## DBSCAN
An implementation of DBSCAN algroithm for comparison purpose
### To execute:
  DBSCAN inputEvents output searchRadius minPts minCorPointsInEachCluster nonCorePoints
### Arguments:
1. inputEvents: input file of control points, a csv without header with two columns: x and y
2. output: output file name
3. searchRadius: search radius to check density and to expand clusters
4. minPts: minimum points to decide core points
5. minCorPointsInEachCluster: minimum number of core points in each cluster
6. nonCorePoints: whether clusters should keep non-core points
  * 0: not keeping
  * 1: keeping

