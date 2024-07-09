#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "shared.h"
#include <pthread.h>
#define NODE_THRESHOLD 1
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISPLACEMENT 1
#define MAX_NEIGHBOURS 4
#define INTERVAL 2
#define PORTS 5
char *directoryName;

nodeLogs *nodeLog;
pthread_mutex_t lockOne = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lockTwo = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lockThree = PTHREAD_MUTEX_INITIALIZER;

int iteration;
int nodeRank;
int nodeCoords[DIMENSIONS];
int portStatus[PORTS] = {0};

volatile int endAll = 0;
volatile int endThread = 0;

MPI_Comm nodeCart;
MPI_Comm universalComm;

int neighbourMessageCount = 0;
int baseMessageCount = 0;

double commTime = 0;

int initialiseNodes(int rank, MPI_Comm communicator, char *directory)
{

    nodeRank = rank;
    if (nodeRank == 0)
        printf("Comm Size: %d: Grid Dimension =[%d x %d] \n", cartSize, dims[0], dims[1]);

    MPI_Cart_create(communicator, DIMENSIONS, dims, periods, REORDER, &nodeCart);

    char makeDirectory[200];
    snprintf(makeDirectory, sizeof(makeDirectory), "mkdir -p %s", directory);
    system(makeDirectory);
    char filename[200];
    snprintf(filename, sizeof(filename), "%s/node%d.txt", directory, nodeRank);

    FILE *file = fopen(filename, "w");
    fprintf(file, "Initialising Log for Node: %d\n", nodeRank);
    fclose(file);
    return EXIT_SUCCESS;
}
void *portSimulation(void *arg)
{
    while (!endAll)
    {
        srand(time(NULL) + (int)arg);
        int sleep_time = rand() % 5;
        sleep(sleep_time);

        // 0 or 1, determing avalability
        int is_available = rand() % 2;
        portStatus[(int)arg] = is_available;
    }
    printf("terminating threads %d\n", (int)arg);
    return NULL;
}

packBaseReport(baseReport report) {
    report.inital = MPI_Wtime();
    printf("Packing report for Base in Node: %d at time: %lf\n", nodeRank, report.inital);
    
    char *buffer;
    int position = 0;

    int bSize, bFreePorts, bNeighbourCount, bNeighbours, bNeighbourFreePorts, bInitial;
    MPI_Pack_size(1, MPI_INT, universalComm, &bFreePorts);
    MPI_Pack_size(1, MPI_INT, universalComm, &bNeighbourCount);
    MPI_Pack_size(4, MPI_INT, universalComm, &bNeighbours);
    MPI_Pack_size(4, MPI_INT, universalComm, &bNeighbourFreePorts);
    MPI_Pack_size(1, MPI_DOUBLE, universalComm, &bInitial);

    bSize = bFreePorts + 
            bNeighbourCount + 
            bNeighbours + 
            bNeighbourFreePorts + 
            bInitial;

    buffer = (char*)malloc(bSize);

    MPI_Pack(&report.freePorts, 1, MPI_INT, buffer, bSize, &position, universalComm);
    MPI_Pack(&report.neighbourCount, 1, MPI_INT, buffer, bSize, &position, universalComm);
    MPI_Pack(&report.neighbours, 4, MPI_INT, buffer, bSize, &position, universalComm);
    MPI_Pack(&report.neighbourFreePorts, 4, MPI_INT, buffer, bSize, &position, universalComm);
    MPI_Pack(&report.inital, 1, MPI_DOUBLE, buffer, bSize, &position, universalComm);
    MPI_Send(buffer, bSize, MPI_PACKED, baseRank, NODETOBASE, universalComm);

    baseMessageCount++;

    free(buffer);
}

void *neighbourReporting(void *arg)
{
    int north, south, west, east;
    double initial, finish;

    // get neighbour nodes
    MPI_Cart_shift(nodeCart, ROW, DISPLACEMENT, &south, &north);
    MPI_Cart_shift(nodeCart, COLUMN, DISPLACEMENT, &west, &east);


    int neighbourRanks[NEIGHBOURS];
        neighbourRanks[0] = south;
        neighbourRanks[1] = north;
        neighbourRanks[2] = west;
        neighbourRanks[3] = east;

    MPI_Status neighbourStatus[NEIGHBOURS];
    MPI_Request neighbourWaiting[NEIGHBOURS],
                neighbourConfirmed[NEIGHBOURS];
    
    int neighbourFreePorts[NEIGHBOURS];
    int existsPorts = 0;
    int message = 1;

    baseReport report;

    char filename[200];
    snprintf(filename, sizeof(filename), "%s/node%d.txt", directoryName, nodeRank);


    initial = MPI_Wtime();

    FILE *file = fopen(filename, "a");
    // Send send and recieves to all neighbours of the reporting node.
    for (int i = 0; i < NEIGHBOURS; i++) {
        if (neighbourRanks[i] != -2) {
            fprintf(file, "Node: %d requesting to Node: %d for port count.\n", nodeRank, neighbourRanks[i]);
        }
        // Both non-blocking send and recieves. Let's function continue. Block exists at Waitall       
        MPI_Isend(&message, 1, MPI_INT, neighbourRanks[i], RESPONSE, nodeCart, &neighbourWaiting[i]);
        MPI_Irecv(&neighbourFreePorts[i], 1, MPI_INT, neighbourRanks[i], RESPONSE, nodeCart, &neighbourConfirmed[i]);
        neighbourMessageCount += 2;
    }
    fclose(file);

    // Baricade. Program will not continue until all neighbour reports recieved
    MPI_Waitall(NEIGHBOURS, &neighbourConfirmed, &neighbourStatus);
    finish = MPI_Wtime();
    commTime += finish - initial;

    printf("Node: %d, Comunication Time: %.2f\n", nodeRank, commTime);
    report.freePorts = nodeLog->logs[nodeLog->end].freePorts;
    report.neighbourCount = 0;

    
    file = fopen(filename, "a");
    // Update reports
    for (int i = 0; i < NEIGHBOURS; i++) {
        if (neighbourRanks[i] != -2) {
            report.neighbours[report.neighbourCount] = neighbourRanks[i];
            report.neighbourFreePorts[report.neighbourCount] = neighbourFreePorts[i];
            report.neighbourCount++;
            // States the base report necessary, instead print neighbour availability.  
            if (neighbourFreePorts[i] > NODE_THRESHOLD) {
                // Flag so that a base report is not sent.
                existsPorts = 1;
                fprintf(file, "Node: %d redirect to Node: %d with %d free ports\n", nodeRank, neighbourRanks[i], neighbourFreePorts[i]);
            }
        }
    }
    fclose(file);

    // If entered previous loop has dictated that there is no neighbour availability.
    if (!existsPorts)
    {   
        packBaseReport(report);
    }

    // ----- MUTEX LOCK -----
    pthread_mutex_lock(&lockThree);
    endThread = 1;
    pthread_mutex_unlock(&lockThree);
    // ----- MUTEX UNLOCK -----

    return EXIT_SUCCESS;
}

int dataUpdate(nodeLogs *dataLog, int availablity)
{   
    nodeReport node;
    time_t initial = time(NULL);
    struct tm *pTime = gmtime(&initial);
    
    node.freePorts = availablity;

    node.year = pTime->tm_year, node.month = pTime->tm_mon, node.day = pTime->tm_mday;
    node.hour = pTime->tm_hour, node.minute = pTime->tm_min, node.second = pTime->tm_sec;
    
    // Prevents spamming of the console.
    if (iteration % 5 == 0) {
        printf("Iteration: %d Node: %d has %d free ports.\n", iteration, nodeRank, node.freePorts);
    }
    int logSize = dataLog -> size;
    if (logSize == FULL) {
        dataLog->start = (dataLog->start + 1) % FULL;
        dataLog->end = (dataLog->end + 1) % FULL;
        dataLog->logs[dataLog->end] = node;
    } else {
        dataLog->end = (dataLog->end + 1) % FULL;
        dataLog->logs[dataLog->end] = node;
        dataLog->size++;
    }
    return EXIT_SUCCESS;
}

int fileOut(int rank, nodeLogs *nodeLog) {

    int i, size;
    char filename[200];
    snprintf(filename, sizeof(filename), "%s/node%d.txt", directoryName, rank);

    FILE *file = fopen(filename, "a");
    size = nodeLog->size;


    if (file) {
        int start = nodeLog->start;
        fprintf(file, "Node: %d Iteration: %d\n", rank, iteration);
        for (i = 0; i < size; i++) {
            int day = nodeLog->logs[start].day;
            int month = nodeLog->logs[start].month;
            int year = nodeLog->logs[start].year;
            int hour = nodeLog->logs[start].hour;
            int minute = nodeLog->logs[start].minute;
            int second = nodeLog->logs[start].second;
            int ports = nodeLog->logs[start].freePorts;
            fprintf(file, "%d/%d/%d %d:%d:%d, %d\n",
                    day, month, year, hour, minute, second, ports);
            start = (start + 1) % 10;
        }
        fclose(file);
    }
    return EXIT_SUCCESS;
}

int nodeAvailability(pthread_t commThread) {
    int i;
    
    // ----- MUTEX LOCK -----
    pthread_mutex_lock(&lockOne);
    int freePorts = 0;
    for (i = 0; i < PORTS; i++) {
        freePorts += 1 - portStatus[i];
    }
    pthread_mutex_unlock(&lockOne);
    // ----- MUTEX UNLOCK -----

    // ----- MUTEX LOCK -----
    pthread_mutex_lock(&lockTwo);
    // here is where we update the logs for the process.
    dataUpdate(nodeLog, freePorts);
    fileOut(nodeRank, nodeLog);
    pthread_mutex_unlock(&lockTwo);
    // ----- MUTEX UNLOCK -----

    // So if the availability stations is less than a desired theshold, get neighbour reports
    if (freePorts <= NODE_THRESHOLD) {
        pthread_create(&commThread, NULL, neighbourReporting, (void *)nodeRank);
        pthread_join(commThread, NULL);
    }

    return EXIT_SUCCESS;
}

int runNodes(MPI_Comm commWorld, char *inputDir)
{
    int i, flag, commFlag, neighbour;

    directoryName = inputDir;
    universalComm = commWorld;
    // Where all the data for the charging logs is kept

    // Creates the cartesian coordinated Comm
    MPI_Cart_coords(nodeCart, nodeRank, DIMENSIONS, nodeCoords);

    // Logs for which hold history
    nodeLog = (nodeLogs *)malloc(sizeof(nodeLogs));
    nodeLog->size = 0, nodeLog->start = 0, nodeLog->end = -1;;
    pthread_t iThreads[PORTS], commThread;

    // Create a thread for each port to update
    for (int i = 0; i < PORTS; i++)
    {
        pthread_create(&iThreads[i], NULL, portSimulation, (void *)i);
    }


    MPI_Status baseStat, neighbourStat;
    MPI_Request baseReq, neighbourReq;

    iteration = 0, commFlag = 0;
    int message, messageTwo, tag;

    // Non-blocking recieve from base, function continues on
    MPI_Irecv(&message, 1, MPI_INT, baseRank, MPI_ANY_TAG, MPI_COMM_WORLD, &baseReq);

    while (SIMULATION) {
        // Tests for the request from Irecv, setting the flag circumstance dependant
        MPI_Test(&baseReq, &flag, &baseStat);
        // Tag of the message recieved
        tag = baseStat.MPI_TAG;

        // Only enters if flagged following MPI_Test, handling all messages from the base
        if (flag) {
            baseMessageCount++;
            // Tag from recieved status
            switch (tag) {
            
            // Message from base recieved by non-blocking recieve.
            case BASETONODE:
                printf("Next available Station: %d, (%d, %d) \n",
                    message, nodeCoords[0], nodeCoords[1]);
                break;
            
            // Message from base to terminate all.
            case FINALISE:
                printf("Exiting Process: %d.\n", nodeRank);
                endAll = 1;
                for (int i = 0; i < PORTS; i++){
                    pthread_join(iThreads[i], NULL);
                }
                break;
            }
            // Reset recievership for base
            MPI_Irecv(&message, 1, MPI_INT, baseRank, MPI_ANY_TAG, MPI_COMM_WORLD, &baseReq);
        } if (endAll) { break; } // Message from base recieved to terminate => break loop.

        // Probing for relevant messages from neighbours. This is non-blocking. Response tag, signifies to neighbours
        MPI_Iprobe(MPI_ANY_SOURCE, RESPONSE, nodeCart, &commFlag, &neighbourStat);
        // So if this flag is true a neighbour process has probed this neighbour for its availability.
        while (commFlag) {   
            int neighbourNode = neighbourStat.MPI_SOURCE;
            // Blocking call but irrelevant as Iprobe has alread dictated that there exists a message.
            MPI_Recv(&messageTwo, 1, MPI_INT, neighbourStat.MPI_SOURCE, RESPONSE, nodeCart, &neighbourStat);
            printf("Node: %d probed for port availability by Node: %d\n", nodeRank, neighbourNode);
            
            // ----- MUTEX LOCK -----
            pthread_mutex_lock(&lockTwo);
            // So this sends the avalability of this MPI_nodes fuel stations back to the requesting node
            // Be mindful of queue indexing here.
            // Mutex lock required here in order to access 
            MPI_Isend(&nodeLog->logs[nodeLog->end].freePorts, 1, MPI_INT, neighbourStat.MPI_SOURCE, RESPONSE, nodeCart, &neighbourReq);
            neighbourMessageCount += 2;
            pthread_mutex_unlock(&lockTwo);
            // ----- MUTEX UNLOCK -----

            commFlag = 0;
            // Flag set to 0 => will exit if there exists no new messages following probe @ next line
            MPI_Iprobe(MPI_ANY_SOURCE, RESPONSE, nodeCart, &commFlag, &neighbourStat);
        }

        // Checks the availability of the node at each iteration
        nodeAvailability(commThread);

        sleep(1);

        // ----- MUTEX LOCK -----
        pthread_mutex_lock(&lockThree);
        if (endThread) {
            pthread_join(commThread, NULL);
            endThread = 0;
        }
        pthread_mutex_unlock(&lockThree);
        // ----- MUTEX UNLOCK -----

        iteration++;
    }
    free(nodeLog);
    char filename[200];
    snprintf(filename, sizeof(filename), "%s/node%d.txt", directoryName, nodeRank);

    FILE *nodeLogs = fopen(filename, "a");

    fprintf(nodeLogs, "Messages between Node: %d and Base: %d\n", nodeRank, baseMessageCount);
    fprintf(nodeLogs, "Messages between Node: %d and Neighbours: %d\n", nodeRank, neighbourMessageCount);
    fprintf(nodeLogs, "Node: %d to Neighbour Nodes Communication Time: %f\n", nodeRank, commTime);

    MPI_Comm_free(nodeCart);

    return EXIT_SUCCESS;
}