#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

typedef struct list ELEMENT;
typedef ELEMENT *LINK;
typedef double *datapoint; 
/* LINK is a pointer to ELEMENT */
/* datapoint is a pointer to a double array */

double** fileToDataPoints(FILE *ifp,int d, int size);
int findDimension(FILE *ifp);
int findInputSize(FILE *ifp);
double **inizializeCentroids(int k, int d, double **datapoints);
double restartClusters(LINK *clusters, int k);
void kMeans(int k, int size, int d, int max_iter, double **cents, LINK *clusters, double **matrix);
void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d);
int updateCentroids(double **cents, LINK *clusters, double** inputMatrix ,int k, int d);
double calculateNorma(double *old, double * new, int d);
double calculateDistance(double *datapoint, double *centroid, int d);

struct list { /*this is an item in a list*/
    //int length;
    int datapoint; 
    struct list *next;   
};

double ** fileToDataPoints(FILE *ifp, int d, int size) {
    double currentCoordinate;
    double *vector = (double *)calloc(size*d, sizeof(double));
    assert(vector != NULL);
    //printf("allocated big contiuous vector\n");
    double **matrix = (double **)calloc(size, sizeof(double *)); 
    assert(matrix != NULL);
    //printf("allocated pointer vector\n");
    for (int i=0; i<size; i++) {
        matrix[i] = vector + i*d;
    }
    //printf("initialized pointers\n");
    for (int i=0; i<size; i++){  
        for(int j=0; j<d; j++){
            fscanf(ifp, "%lf %*[,]", &currentCoordinate);
            matrix[i][j] = currentCoordinate;
        }
    }
    //printf("method finished\n");
    return matrix;
    }

int findDimension(FILE *ifp) {
    char c;
    int commacnt=0;
    int dim=0;
    do {
        c = fgetc(ifp);
        // Taking input single character at a time
        if( feof(ifp) ){
            return -1; 
        }

        if (c == 10) // c==\n 
            {
                dim = commacnt+1;
                break;
            }
        else if (c==44){
            // read first line #of, +1 = d, go to the beggining and create array of arrays (d*sizeof(double))
            // if digit , dot or - => add to the corrent cordinate of (X_i)
            // else if , => add to 
                commacnt++;
            } 
        } while(1); 
    
    rewind(ifp);
    return dim;
}

int findInputSize(FILE *ifp) {
    char c;
    int cnt;
    do {
        c = fgetc(ifp);
        /*printf("%c",c)*/;
        // Taking input single character at a time
        if( feof(ifp) ){
            break; 
        }

        if (c == 10) // c==\n 
            {
                cnt++;
            }
        } while(1); 
    
    rewind(ifp);
    return cnt;
}

double ** inizializeCentroids(int k, int d, double **datapoints){
        /*if there are less vectors in the input than k than print "invalid input" and terminate*/
    double *centroid = (double *)calloc(k*d, sizeof(double));
    assert(centroid != NULL);
    //printf("allocated big contiuous vector\n");
    double **matrix = (double **)calloc(k, sizeof(double *)); 
    assert(matrix != NULL);
    //printf("allocated pointer vector\n");
    for (int i=0; i<k; i++) {
        matrix[i] = centroid + i*d;
    }
    for (int i = 0; i<k; i++){
        for (int j = 0; j<d; j++){
            matrix[i][j] = datapoints[i][j]; 
        }
    }  
    return matrix; 
}

double restartClusters(LINK *clusters, int k) {
    for (int i=0; i<k;i++){
        LINK current = clusters[i];
        while (current != NULL) {
            if (current->next != NULL) {
                LINK next = current->next;
                free(current);
                current = next;
            } else {
                free(current);
            }
        }
        //clusters[i]->length = 0;
        clusters[i]->datapoint = -1;
        clusters[i]->next = NULL;
    }
}

void kMeans(int k, int size, int d, int max_iter, double **cents, LINK *clusters, double **matrix){
    int iter = 0;
    int continue_condition = 1;
    while (iter<max_iter && continue_condition){
        /*
        1) for each item in the input list
            2) for each centroid
                3) calculate distance 
                4) update minimum distance and cluster _i_ from which the distance is minimal
            5) add item to the list in cluster _i_

        6) update centroids + empty cluster list and free the space
        7) update iter
        8) update continue_codition
        */
       assignToCluster(cents, matrix, clusters, k, size, d);
       continue_condition = updateCentroids(cents, clusters, matrix, k, d);
       iter++;
    }
}

void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d) {
    double distance = DBL_MAX;
    double tempDist;
    int clusterIndex;
    /*
        1) for each item in the input list
            2) for each centroid
                3) calculate distance 
                4) update minimum distance and cluster _i_ from which the distance is minimal
            5) add item to the list in cluster _i_
        */
    for(int i=0; i<size; i++){
        for (int j=0; j<k; j++){
            tempDist = calculateDistance(datapoints[i], cents[j], d);
            if (tempDist<distance) {
                distance = tempDist;
                clusterIndex = j;
            }
        }
        if (clusters[clusterIndex]->datapoint == -1){
            clusters[clusterIndex]->datapoint = i;
        }else {
            LINK current = clusters[clusterIndex];
            LINK new = NULL;
            new -> datapoint = i;
            new -> next = current;
            clusters[clusterIndex] = new;
        }

    }
}

int updateCentroids(double **cents, LINK *clusters, double** inputMatrix ,int k, int d) {
    double difference = 0;
    double *sum = (double *)calloc(d, sizeof(int));
    double epsilon = .001;
    LINK current;
    int sizeOfCluster = 0;
    for (int i=0; i<k; i++) {
        current = clusters[i];
        while (current != NULL){
            for (int j=0; j<d; j++) {
                sum[j] += inputMatrix[current->datapoint][j];
            }
            current = current->next;
            sizeOfCluster++;
        }
        for (int m = 0; m<sizeOfCluster; m++){
            sum[m] /= sizeOfCluster;
        }
        int check = (int)(calculateNorma(cents[i], sum, d) < epsilon);
        difference += (calculateNorma(cents[i], sum, d) < epsilon)?1:0;
        cents[i] = sum;
    }
    restartClusters(clusters, k);
    free(sum);
    return difference;
}

double calculateNorma(double *old, double * new, int d){
    double sum = 0;
    for (int i=0; i<d; i++){
        sum += pow(old[i]-new[i],2);
    }
    return sqrtf(sum);
}

double calculateDistance(double *datapoint, double *centroid, int d){
    double distance = 0;
    for (int i = 0; i<d; i++){
        distance+= pow(datapoint[i]-centroid[i],2);
    }
    return distance;
}

int main(int argc, char *argv[]){
    //check validity of inputs- if not "Invalid Input!"
    /*IMPORTANT! - we need to make sure that the input is read corrrectly*/
    int K = atoi(argv[1]);
    printf("K=%d\n",K);
    char *inputFilename=argv[2];
    char *outputFilename=argv[3];
    int max_iter = 200; 
    
    if (argc==5){
        max_iter = atoi(argv[4]);
    }
    printf("max_iter=%d\n",max_iter);
    FILE *ifp=NULL;
    FILE *outfile=NULL;
    int N=200; //change later 

    ifp = fopen(inputFilename, "r");
    assert(ifp!=NULL);
    printf("Opened input file successfuly\n");
    int d = findDimension(ifp);
    printf("d = %d\n",d);
    int size = findInputSize(ifp);
    printf("size = %d\n",size);

    double** datapointMatrix = fileToDataPoints(ifp, d, size);
    printf("%lf\n",datapointMatrix[499][2]);
    fclose(ifp);
    printf("Closed input file successfully\n");
    
    double **centroids = inizializeCentroids(K, d, datapointMatrix); //this array holds K datapoints
    printf("Inititialized centroids\n");
    LINK *clusters = (LINK *)calloc(K, sizeof(LINK));
    assert(clusters != NULL);
    restartClusters(clusters, K);
    printf("Initialized clusters");
    kMeans(K, size, d, max_iter, centroids, clusters, datapointMatrix);
    printf("Successfully ran kMeans\n");
    printf("First index of the first centroid is %lf\n",centroids[0][0]);
    printf("Writing to output file...");
    outfile = fopen(outputFilename,"w");
    assert(outfile != NULL);
    printf("opened outfile successfully\n");
    fclose(outfile);
    printf("closed outfile successfully\n");
    printf("Done!\n");
    return 0;
}