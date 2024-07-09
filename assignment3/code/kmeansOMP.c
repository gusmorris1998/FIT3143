#include <omp.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#define TOLERANCE 0

int globalPointCount, globalThreads, globalIterations, globalK, globalClusterChangeCount;
int *globalDataPoints;
float *centroidsGlobal;
float globalToleranceValue;
int *globalClusterID;
int **globalClusterCount;

void runKmeansThread(int *threadID)
{
    int *id = (int *)threadID;

    // Assigning data points range to each thread
    float toleranceValue = 100;
    int threadClusterChangeCount;
    int threadLength = globalPointCount / globalThreads;
    int start = (*id) * threadLength;
    int end = start + threadLength;
    if (end + threadLength > globalPointCount) {
        //To assign last undistributed points to this thread for computation, change end index to number_of_points_global
        end = globalPointCount;
        threadLength = globalPointCount - start;
    }

    int i, j;
    double minDist, currentDist;
    int *pointClusterID = (int *)malloc(threadLength * sizeof(int));

    float *clusterSum = (float *)malloc(globalK * 3 * sizeof(float));
    int *clusterCount = (int *)malloc(globalK * sizeof(int));

    // Start of loop
    int iterationCount = 0;
    while ((toleranceValue > TOLERANCE) && (iterationCount < 200)) {
        globalClusterChangeCount = 0;
        // Initialize cluster_points_sum or centroid to 0.0
        for (i = 0; i < globalK * 3; i++) {
            clusterSum[i] = 0.0;
        }

        // Initialize number of points for each cluster to 0
        for (i = 0; i < globalK; i++) {
            clusterCount[i] = 0;
        }

        threadClusterChangeCount = 0;
        for (i = start; i < end; i++) {
            //Assign these points to their nearest cluster
            minDist = DBL_MAX;
            int initialCluster = pointClusterID[i - start];
            for (j = 0; j < globalK; j++) {
                int indexStart = (iterationCount * globalK + j) * 3;
                currentDist = pow((double)(centroidsGlobal[indexStart] - (float)globalDataPoints[i * 3]), 2.0) + pow((double)(centroidsGlobal[indexStart + 1] - (float)globalDataPoints[i * 3 + 1]), 2.0) + pow((double)(centroidsGlobal[indexStart + 2] - (float)globalDataPoints[i * 3 + 2]), 2.0);
                if (currentDist < minDist) {
                    minDist = currentDist;
                    pointClusterID[i - start] = j;
                }
            }
            if (pointClusterID[i - start] != initialCluster) {
                threadClusterChangeCount += 1;
            }
            
            clusterCount[pointClusterID[i - start]] += 1;
            clusterSum[pointClusterID[i - start] * 3] += (float)globalDataPoints[i * 3];
            clusterSum[pointClusterID[i - start] * 3 + 1] += (float)globalDataPoints[i * 3 + 1];
            clusterSum[pointClusterID[i - start] * 3 + 2] += (float)globalDataPoints[i * 3 + 2];
        }

        #pragma omp critical
        {
            for (i = 0; i < globalK; i++)
            {
                centroidsGlobal[((iterationCount + 1) * globalK + i) * 3] = (centroidsGlobal[((iterationCount + 1) * globalK + i) * 3] * globalClusterCount[iterationCount][i] + clusterSum[i * 3]) / (float)(globalClusterCount[iterationCount][i] + clusterCount[i]);
                centroidsGlobal[((iterationCount + 1) * globalK + i) * 3 + 1] = (centroidsGlobal[((iterationCount + 1) * globalK + i) * 3 + 1] * globalClusterCount[iterationCount][i] + clusterSum[i * 3 + 1]) / (float)(globalClusterCount[iterationCount][i] + clusterCount[i]);
                centroidsGlobal[((iterationCount + 1) * globalK + i) * 3 + 2] = (centroidsGlobal[((iterationCount + 1) * globalK + i) * 3 + 2] * globalClusterCount[iterationCount][i] + clusterSum[i * 3 + 2]) / (float)(globalClusterCount[iterationCount][i] + clusterCount[i]);
                
                globalClusterCount[iterationCount][i] += clusterCount[i];
            }
            globalClusterChangeCount += threadClusterChangeCount;
        }

        #pragma omp barrier
        if (*id == 0)
        {
            globalIterations++;
            globalToleranceValue = (float)globalClusterChangeCount / (float)globalPointCount;
            // printf("%d / %d = %f\n", globalClusterChangeCount, globalPointCount, toleranceValue);
        }

    #pragma omp barrier
        toleranceValue = globalToleranceValue;
        iterationCount++;
    }



    for (i = start; i < end; i++) {
        globalClusterID[i * 4] = globalDataPoints[i * 3];
        globalClusterID[i * 4 + 1] = globalDataPoints[i * 3 + 1];
        globalClusterID[i * 4 + 2] = globalDataPoints[i * 3 + 2];
        globalClusterID[i * 4 + 3] = pointClusterID[i - start];
    }
}

void ompKmeans(int threadCount,int N,int K,int *dataPoints,int **clusterID,float **centroids,int *iterations){

    globalPointCount = N;
    globalThreads = threadCount;
    globalIterations = 0;
    globalK = K;
    globalDataPoints = dataPoints;

    *clusterID = (int *)malloc(N * 4 * sizeof(int));
    globalClusterID = *clusterID;

    centroidsGlobal = (float *)calloc((201) * K * 3, sizeof(float));

    int i;
    for (i = 0; i < K; i++) {
        centroidsGlobal[i * 3] = dataPoints[i * 3];
        centroidsGlobal[i * 3 + 1] = dataPoints[i * 3 + 1];
        centroidsGlobal[i * 3 + 2] = dataPoints[i * 3 + 2];
    }

    globalClusterCount = (int **)malloc(200 * sizeof(int *));
    for (i = 0; i < 200; i++) {
        globalClusterCount[i] = (int *)calloc(K, sizeof(int));
    }

    omp_set_num_threads(threadCount);
    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        runKmeansThread(&ID);
    }

    *iterations = globalIterations;

    int centroidSize = (*iterations + 1) * K * 3;
    printf("ITERATIONS:%d\n", *iterations);
    *centroids = (float *)calloc(centroidSize, sizeof(float));
    for (i = 0; i < centroidSize; i++) {
        (*centroids)[i] = centroidsGlobal[i];
    }

}

void readData(const char *fileName, int *N, int **points) {\
    int i;
	FILE *fin = fopen(fileName, "r");
	fscanf(fin, "%d", N);
	*points = (int *)malloc(sizeof(int) * ((*N) * 3));
	for (i = 0; i < (*N) * 3; i++) {
		fscanf(fin, "%d", (*points + i));
	}
	fclose(fin);
}

void printClusters(const char *cluster_filename, int N, int *cluster_points)
{
	FILE *fout = fopen(cluster_filename, "w");
    int i = 0;
	for (i = 0; i < N; i++) {
		fprintf(fout, "%d %d %d %d\n",*(cluster_points + (i * 4)), *(cluster_points + (i * 4) + 1), *(cluster_points + (i * 4) + 2), *(cluster_points + (i * 4) + 3));
	}
	fclose(fout);
}

void printCentroids(const char *fileName, int K, int iterationCount, float *centroids)
{   
    int i, j;
	FILE *fout = fopen(fileName, "w");
	for (i = 0; i < iterationCount + 1; i++) {
		for (j = 0; j < K; j++) {
			fprintf(fout, "%f %f %f, ", *(centroids + (i * K + j) * 3), *(centroids + (i * K + j) * 3 + 1),  *(centroids + (i * K + j) * 3 + 2));
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
}

int main()
{

	int N, K, iterations;
	int threadCount;
	int *dataPoints, *clusterPoints;
	float* centroids;

    K = 30;
    threadCount = 8;
    char *fileName = "dataset.txt";

	double initial, final, time;

	readData(fileName, &N, &dataPoints);

	initial = omp_get_wtime();
	ompKmeans(threadCount, N, K, dataPoints, &clusterPoints, &centroids, &iterations);
	final = omp_get_wtime();

    char clusterFile[105] = "clusterPointsOMP.txt";
    char centroidFile[105] = "centroidPointsOMP.txt";

	printClusters (clusterFile, N, clusterPoints);
	printCentroids (centroidFile, K, iterations, centroids);

   	time = final - initial;
	printf("Time Taken: %lf \n", time);

    return EXIT_SUCCESS;
	
}