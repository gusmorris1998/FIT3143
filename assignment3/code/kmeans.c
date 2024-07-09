#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define min(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define ITERATIONS 300
#define TOLERANCE 0

int globalPointCount, globalIterations, globalK;
int *globalDataPoints;
int *globalClusterID;
float *centroidsGlobal;

void runKmeans()
{

    int clusterChangeCount; // Used for calculating the tolerance
    int i , j;
    double minDist, currentDist;
    float toleranceValue = 100;

	// Cluster id associated with each point
    int *pointClusterID = (int *)malloc(globalPointCount * sizeof(int));

	// Sum of the cluster values
    float *clusterSum = (float *)malloc(globalK * 3 * sizeof(float));

	// No. of points in a cluster.
    int *clusterCount = (int *)malloc(globalK * sizeof(int));

	// Start of loop
    int iterationCount = 0;

    printf("Start K-means algorithm\n");

    while ((toleranceValue > TOLERANCE) && (iterationCount < ITERATIONS)) {
		// Reset the sum array at each iteration
        for (i = 0; i < globalK * 3; i++) {
            clusterSum[i] = 0.0;
        }

		// Reset clusterCount array
        for (i = 0; i < globalK; i++) {
            clusterCount[i] = 0;
        }

        clusterChangeCount = 0;
        // Allocate each data point to a cluster dependant on proximity to each clusters centroid.

        //Iterate through all points
        for (i = 0; i < globalPointCount; i++) {
            minDist = DBL_MAX;
            int initialCluster = pointClusterID[i];
            // For all data points iterate through all centroids and determine distance between them
            for (j = 0; j < globalK; j++) {
                int indexStart = (iterationCount * globalK + j) * 3;
                currentDist = pow((double)(centroidsGlobal[indexStart] - (float)globalDataPoints[i * 3]), 2.0) +
                               pow((double)(centroidsGlobal[indexStart + 1] - (float)globalDataPoints[i * 3 + 1]), 2.0) +
                               pow((double)(centroidsGlobal[indexStart + 2] - (float)globalDataPoints[i * 3 + 2]), 2.0);
                // ClusterID and distance will be designated accordindly by shortest distance.
                if (currentDist < minDist) {
                    minDist = currentDist;
                    pointClusterID[i] = j;

                }
            }
            // Update the count depending if the point has changed clusters.
            if (pointClusterID[i] != initialCluster) {
                clusterChangeCount += 1;
            }

            // Dictates convergence.
            toleranceValue = (float)clusterChangeCount / (float)globalPointCount;

            // Number of points within the cluster.
            clusterCount[pointClusterID[i]] += 1;

			// Next three lines compute the total sum of the points within the cluster.
            clusterSum[pointClusterID[i] * 3] += (float)globalDataPoints[i * 3];
            clusterSum[pointClusterID[i] * 3 + 1] += (float)globalDataPoints[i * 3 + 1];
            clusterSum[pointClusterID[i] * 3 + 2] += (float)globalDataPoints[i * 3 + 2];
        }

        printf("%d / %d\n", clusterChangeCount, globalPointCount);
        printf("%f\n", toleranceValue);

        // Compute the mean point of the cluster by taking the sum of the cluster by the amount of points in the cluster. This is the
            // new centroid
        for (i = 0; i < globalK; i++) {
            int indexStart = ((iterationCount + 1) * globalK + i) * 3;
            centroidsGlobal[indexStart] = clusterSum[i * 3] / clusterCount[i];
            centroidsGlobal[indexStart + 1] = clusterSum[i * 3 + 1] / clusterCount[i];
            centroidsGlobal[indexStart + 2] = clusterSum[i * 3 + 2] / clusterCount[i];
        }

        iterationCount++;
    }

    globalIterations = iterationCount;

    // Assign points to final choice for cluster centroids
    for (i = 0; i < globalPointCount; i++) {
        // Assign points to clusters
        globalClusterID[i * 4] = globalDataPoints[i * 3];
        globalClusterID[i * 4 + 1] = globalDataPoints[i * 3 + 1];
        globalClusterID[i * 4 + 2] = globalDataPoints[i * 3 + 2];
        globalClusterID[i * 4 + 3] = pointClusterID[i];
    }
}

void SequentialKmeans(int N, int K, int* dataPoints, int** ClusterID,float** stepCentroids,int* iterations) {

    // Initialize global variables, these can be used across functions
    int i;
    globalPointCount = N;
    globalIterations = *iterations;
    globalK = K;
    globalDataPoints = dataPoints;

	// Each 4 bit block is made up of (x, y, z, clusterID)
    *ClusterID = (int *)malloc(N * 4 * sizeof(int));
    globalClusterID = *ClusterID;

   // Space allocation for centroids at each iteration. 
    centroidsGlobal = (float *)calloc((ITERATIONS + 1) * K * 3, sizeof(float));

    // Randomly assign centroids. Random data set hence assign first k-data points
        // This also allows comparison between parallel and sequential implementation.
    for (i = 0; i < K; i++) {
        centroidsGlobal[i * 3] = dataPoints[i * 3];
        centroidsGlobal[i * 3 + 1] = dataPoints[i * 3 + 1];
        centroidsGlobal[i * 3 + 2] = dataPoints[i * 3 + 2];
    }

    // Run function, this is the point outside initilisations in which parallization can occur.
    runKmeans();

    *iterations = globalIterations;
    int centroids = (*iterations + 1) * K * 3;
    *stepCentroids = (float *)calloc(centroids, sizeof(float));

    for (i = 0; i < centroids; i++) {
        (*stepCentroids)[i] = centroidsGlobal[i];
    }

    printf("ITERATIONS:%d\n", globalIterations);
}

void readData(const char *fileName, int *N, int **points) {
    int i;
	FILE *fin = fopen(fileName, "r");
	fscanf(fin, "%d", N);
	*points = (int *)malloc(sizeof(int) * ((*N) * 3));
	for (i = 0; i < (*N) * 3; i++) {
		fscanf(fin, "%d", (*points + i));
	}
	fclose(fin);
}

/*
    Cluster points file printing function.
    Line has form pointi.x,     pointi.y,   pointi.z,   clusterID
*/
void printClusters(const char *fileName, int N, int *clusterPoints) {
	FILE *fout = fopen(fileName, "w");
    int i;
	for (i = 0; i < N; i++) {
		fprintf(fout, "%d %d %d %d\n", *(clusterPoints + (i * 4)), *(clusterPoints + (i * 4) + 1), *(clusterPoints + (i * 4) + 2), *(clusterPoints + (i * 4) + 3));
	}
	fclose(fout);
}

/*
    Centroid points file printing function.
    Line has form iterationi.centroid1.x iterationi.centroid1.y iterationi.centroid1.z  ...  iterationi.centroidj.x iterationi.centroidj.y iterationi.centroidj.z  ... iterationi.centroidk.x iterationi.centroidk.y iterationi.centroidk.z
*/
void printCentroids(const char *centroidsFile, int K, int iterations, float *centroids) {
    int i, j;
	FILE *fout = fopen(centroidsFile, "w");
	for (i = 0; i < iterations + 1; i++) {
		for (j = 0; j < K; j++) {
			fprintf(fout, "%f %f %f, ", *(centroids + (i * K + j) * 3), *(centroids + (i * K + j) * 3 + 1), *(centroids + (i * K + j) * 3 + 2));
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
}

int main() {

	int N, iterations; // N: number of data points
	int *dataPoints, *clusterPoints; // Array of data points, Array of cluster points
	float* centroids;

	double intital, final, time;

    // Number of clusters in which we will organise the data.
    int K = 30;

    // Dataset and function to read the data
	char *fileName = "dataset.txt";
	readData(fileName, &N, &dataPoints);

	intital = omp_get_wtime();
    // K-means function.
	SequentialKmeans(N, K, dataPoints, &clusterPoints, &centroids, &iterations);
	final = omp_get_wtime();	

    char clusterFile[105] = "clusterPoints.txt";
    char centroidFile[105] = "centroidPoints.txt";

	printClusters (clusterFile, N, clusterPoints);
	printCentroids (centroidFile, K, iterations, centroids);

   	time = final - intital;
	printf("Time Taken: %lf \n", time);
	
	return EXIT_SUCCESS;
}