#include <mpi.h>

#define DIMENSIONS 2
#define REORDER 1
#define ROW 0
#define COLUMN 1
#define SIMULATION 1
#define DISPLACEMENT 1
#define NEIGHBOURS 4
#define ITERATIONS 50
#define FULL 10

extern int dims[DIMENSIONS], coord[DIMENSIONS];
extern int periods[DIMENSIONS];
extern int n, m, baseRank;
extern int size, cartSize;

struct
{
    int rank;
    int freePorts;
    MPI_Request sent;
    MPI_Request recieved;
    MPI_Status status;

} typedef neighbourStatus;

struct
{
    int freePorts;
    int neighbourCount;
    int neighbours[4];
    int neighbourFreePorts[4];
    double inital;

} typedef baseReport;

struct
{

    int second;
    int minute;
    int hour;
    int day;
    int month;
    int year;
    int freePorts
    
} typedef nodeReport;

struct
{
    int size;
    int start;
    int end;
    nodeReport logs[10];
} typedef nodeLogs;

enum Tags
{
    REQUEST,
    RESPONSE,
    FINALISE,
    NODETOBASE,
    BASETONODE,
    NODETONODE,
};

int initialiseNodes(int rank, MPI_Comm existing_comm, char *dirname);
int runNodes(MPI_Comm master_comm, char *input_dirname);
int runBase(char *dirname);
void getCoordinates(int rank, int *dims, int *coords);
int getSize(int rank, int *dims, int *neighbors);