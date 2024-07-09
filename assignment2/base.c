#include <pthread.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "shared.h"
#include <unistd.h>

int finalise = 0;
FILE *file;

int iteration = 0;
int nodeMessages = 0;

int *reportedNodes;

double cTime = 0;

int calculateNearest(MPI_Status status, baseReport report, int nodeRank) {
    int i, j;
    int coordinates[DIMENSIONS], neighbourCount, neighbours[NEIGHBOURS], freePorts[cartSize];

    // Iterate through neighbours
    for (i = 0; i < report.neighbourCount; i++) {
        coordinates[0] = report.neighbours[i] / dims[1];
        coordinates[1] = report.neighbours[i] % dims[1];

        neighbourCount = getSize(report.neighbours[i], dims, neighbours);

        fprintf(file, "Neighbour Node: %d, with Coordinates: (%d,%d) has: %d free ports.\n",
                report.neighbours[i], coordinates[0], coordinates[1], report.neighbourFreePorts[i]);

        // Itererate through neighbours, neighbours
        for (j = 0; j < neighbourCount; j++) {
            fprintf(file, "Node: %d, within proximity\n", neighbours[j]);
            // Sets the relevant nodes to whether they have availability
            if (reportedNodes[neighbours[j]] != 0) {
                freePorts[neighbours[j]] = 1;
            }
        }
    }

    // Prints to logs, where the cars should be sent to.
    fprintf(file, "Send Cars to: ");
    for (i = 0; i < cartSize; i++) {
        if (freePorts[i] == 1) {
            fprintf(file, " %d, ", i);
            // Also send a message to the node, so that it can print to console.
            MPI_Isend(&freePorts[i], 1, MPI_INT, nodeRank, BASETONODE, MPI_COMM_WORLD, &status);
            nodeMessages++;
        }
    }

    return EXIT_SUCCESS;
}

void *nodeCommunication(void *arg) {
    // The node for which this thread is running for.
    int node = (int)arg;

    baseReport report;
    time_t initial;

    MPI_Status status;
    char *buffer;
    int pos, gotMail, i, j;
    int nodeCoords[DIMENSIONS];
    int bSize, bFreePorts, bNeighbourCount, bNeighbours, bNeighbourFreePorts, bInitial;
    double programInitial, programEnd, programTime;

    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &bFreePorts);
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &bNeighbourCount);
    MPI_Pack_size(4, MPI_INT, MPI_COMM_WORLD, &bNeighbours);
    MPI_Pack_size(4, MPI_INT, MPI_COMM_WORLD, &bNeighbourFreePorts);
    MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &bInitial);

    bSize = bFreePorts + 
            bNeighbourCount + 
            bNeighbours + 
            bNeighbourFreePorts + 
            bInitial;
    buffer = (char *)malloc(bSize);

    gotMail = 0;
    // Consistent running loop similar to nodes.c loop
    while (SIMULATION) {
        // Enters if the iterations have run there course.
        if (finalise) {
            break;
        }

        // Non-blocking probe from the relevant node corresponding to this thread, if there exists a message waiting for receptions
            // ticks the flag for which it can then enter.
        MPI_Iprobe(node, NODETOBASE, MPI_COMM_WORLD, &gotMail, &status);
        // Flag ticked => enter
        if (gotMail) {
            initial = time(NULL);
            pos = 0;
            // Update the array such that this nodes has sent a report => low availability.
            reportedNodes[node] = 1;
            // Recieve the buffer from the node.
            MPI_Recv(buffer, bSize, MPI_PACKED, node, NODETOBASE, MPI_COMM_WORLD, &status);
            struct tm *myTime = localtime(&initial);
            char timeB[100];
            strftime(timeB, sizeof(timeB), "%a %Y-%m-%d %H:%M:%S", myTime);
            
            fprintf(file, "----------- Base Report from Node: %d at Iteration: %d ----------- \n",
                    status.MPI_SOURCE, iteration);
            fprintf(file, "Reporting Node: %d at Coordinates: (%d,%d) at time: %s\n", 
                    status.MPI_SOURCE, status.MPI_SOURCE / dims[1], status.MPI_SOURCE % dims[1], timeB);

            // Unpacks the buffer
            MPI_Unpack(buffer, bSize, &pos, &report.freePorts, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, bSize, &pos, &report.neighbourCount, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, bSize, &pos, &report.neighbours, 4, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, bSize, &pos, &report.neighbourFreePorts, 4, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, bSize, &pos, &report.inital, 1, MPI_DOUBLE, MPI_COMM_WORLD);
            

            programInitial = report.inital;
            programEnd = MPI_Wtime();
            programTime = fabs(programEnd - programInitial);
            cTime += programTime;

            printf("Base node processing time: %f\n", cTime);

            fprintf(file, "Reporting Node %d number of adjacent node: %d available port: %d\n", status.MPI_SOURCE, report.neighbourCount, report.freePorts);

            // Where to send the node based on previous information.
            calculateNearest(status, report, node);

            fprintf(file, "\n");
            fprintf(file, "Program Time: %.2f\n", programTime);
        }
        gotMail = 0;
        sleep(1);
        if (iteration % 5 == 0) {
            printf("Awaiting Response: %d\n", node);
        }
    }
    free(buffer);
    return EXIT_SUCCESS;
}

int runBase(char *dirname) {
    int i, final;
    reportedNodes = (int *)malloc(sizeof(int) * cartSize);


    char buffer[200], filename[200];
    snprintf(buffer, sizeof(buffer), "mkdir -p %s", dirname);
    system(buffer);
    snprintf(filename, sizeof(filename), "%s/base.txt", dirname);

    file = fopen(filename, "w");
    pthread_t iThread[cartSize];

    // For each node within the cart, create a seperate thread for communication
    for (i = 0; i < cartSize; i++) {
        pthread_create(&iThread[i], NULL, nodeCommunication, (void *)i);
    }

    // Run to certain iteration
    for (i = 0; i < ITERATIONS; i++) {
        // Prevents console spam
        if (i % 5 == 0) {
            printf("------------- Program Iteration: %d -------------\n", iteration);
        }        
        iteration++;
        sleep(1);
    }
    finalise = 1;

    for (i = 0; i < cartSize; i++) {
        pthread_join(iThread[i], NULL);
    }
    
    // SEND SIGNAL TO NODES TO TERMINATE.
    final = FINALISE;
    for (i = 0; i < baseRank; i++) {
        MPI_Send(&final, 1, MPI_INT, i, FINALISE, MPI_COMM_WORLD);
        nodeMessages++;
    }
    fprintf(file, "Messages between Nodes and Base: %d, Communication time: %.2f\n", nodeMessages, cTime);
    fprintf(file, "--------------------- END OF BASE REPORT ---------------------\n");

    fclose(file);
    
    free(reportedNodes);
    return EXIT_SUCCESS;
}

int getSize(int rank, int *dimensions, int *neighbors) {
    int coordinates[DIMENSIONS];
    int size = 0;

    coordinates[0] = rank / dimensions[1];
    coordinates[1] = rank % dimensions[1];

    if (coordinates[0] < dimensions[0] - 1) {
        coordinates[0] = coordinates[0] + 1;
        neighbors[size++] = coordinates[0] * dimensions[1] + coordinates[1];
    }

    if (coordinates[0] > 0) {
        coordinates[0] = coordinates[0] - 1;
        neighbors[size++] = coordinates[0] * dimensions[1] + coordinates[1];
        coordinates[0] = coordinates[0] + 1;
    }

    if (coordinates[1] < dimensions[1] - 1) {
        coordinates[1] = coordinates[1] + 1;
        neighbors[size++] = coordinates[0] * dimensions[1] + coordinates[1];
        coordinates[1] = coordinates[1] - 1;
    }

    if (coordinates[1] > 0) {
        coordinates[1] = coordinates[1] - 1;
        neighbors[size++] = coordinates[0] * dimensions[1] + coordinates[1];
        coordinates[1] = coordinates[1] + 1;
    }
    return size;
}