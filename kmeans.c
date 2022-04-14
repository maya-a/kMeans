#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <ctype.h>

typedef struct list ELEMENT;
typedef ELEMENT *LINK;
typedef double *datapoint; 
/* LINK is a pointer to ELEMENT */
/* datapoint is a pointer to a double array */

double** fileToDataPoints(FILE *ifp,int d, int size);
int findDimension(FILE *ifp);
int findInputSize(FILE *ifp);
double **inizializeCentroids(int k, int d, double **datapoints);
void restartClusters(LINK *clusters, int k);
void delete_list(LINK head);
void kMeans(int k, int size, int d, int max_iter, double **cents, LINK *clusters, double **matrix);
void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d);
int updateCentroids(double **cents, LINK *clusters, double** inputMatrix ,int k, int d);
double calculateNorma(double *old, double * new, int d);
double calculateDistance(double *datapoint, double *centroid, int d);
int isNatural(char * s);

struct list { /*this is an item in a list*/
    int datapoint; 
    struct list *next;   
};

double ** fileToDataPoints(FILE *ifp, int d, int size) {
    double currentCoordinate;
    double *vector = NULL;
    double **matrix = NULL;
    int i, j;
    
    vector = (double *)calloc(size*d, sizeof(double));
    assert(vector != NULL && "An Error Has Occurred");
    matrix = (double **)calloc(size, sizeof(double *)); 
    assert(matrix != NULL && "An Error Has Occurred");
    for (i=0; i<size; i++) {
        matrix[i] = vector + i*d;
    }
    for (i=0; i<size; i++){  
        for(j=0; j<d; j++){
            fscanf(ifp, "%lf %*[,]", &currentCoordinate);
            matrix[i][j] = currentCoordinate;
        }
    }
    return matrix;
    }

int findDimension(FILE *ifp) {
    char c;
    int commacnt=0;
    int dim=0;
    do {
        c = fgetc(ifp);
        /* Taking input single character at a time*/
        if( feof(ifp) ){
            return -1; 
        }

        if (c == 10) /* c==\n */
            {
                dim = commacnt+1;
                break;
            }
        else if (c==44){
            /* read first line #of, +1 = d, go to the beggining and create array of arrays (d*sizeof(double))
               if digit , dot or - => add to the corrent cordinate of (X_i)
               else if , => add to */
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
        /* Taking input single character at a time*/
        if( feof(ifp) ){
            break; 
        }

        if (c == 10) /* c==\n */ 
            {
                cnt++;
            }
        } while(1); 
    
    rewind(ifp);
    return cnt;
}

double ** inizializeCentroids(int k, int d, double **datapoints){
    /*if there are less vectors in the input than k than print "invalid input" and terminate*/
    double *centroid = NULL;
    double **matrix = NULL;
    int i, j;
    centroid = (double *)calloc(k*d, sizeof(double));
    assert(centroid != NULL && "An Error Has Occurred");
    matrix = (double **)calloc(k, sizeof(double *)); 
    assert(matrix != NULL && "An Error Has Occurred");
    for (i=0; i<k; i++) {
        matrix[i] = centroid + i*d;
    }
    for (i=0; i<k; i++){
        for (j = 0; j<d; j++){
            matrix[i][j] = datapoints[i][j]; 
        }
    }  
    return matrix; 
}

void restartClusters(LINK *clusters, int k) {
    int i;
    LINK cluster;
    for (i=0; i<k;i++) {
        cluster = clusters[i];
        delete_list(cluster);
        /*while (current != NULL) {
            LINK next = current->next;
            current = NULL;
            free(current);
            current = next;
        }*/
        printf("nullified  cluster %d\n", i);
        clusters[i] = (ELEMENT*)malloc(sizeof(ELEMENT));
        assert(clusters[i] != NULL && "An Error Has Occurred");
        clusters[i]->datapoint = -1;
        clusters[i]->next = NULL;
    }
}
void delete_list(LINK head) {
    if (head != NULL){
        printf("%d",head != NULL);
        delete_list(head->next);
        free(head);
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
       printf("assigning datapoints to clusters: %d/%d\n", (iter+1), max_iter);
       assignToCluster(cents, matrix, clusters, k, size, d);
       printf("updating centroids\n");
       continue_condition = updateCentroids(cents, clusters, matrix, k, d);
       iter++;
    }
}

void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d) {
    double distance = DBL_MAX;
    double tempDist;
    int clusterIndex, i, j;
    LINK current, new;
    /*
        1) for each item in the input list
            2) for each centroid
                3) calculate distance 
                4) update minimum distance and cluster _i_ from which the distance is minimal
            5) add item to the list in cluster _i_
        */
    for(i=0; i<size; i++){
        for (j=0; j<k; j++){
            tempDist = calculateDistance(datapoints[i], cents[j], d);
            if (tempDist<distance) {
                distance = tempDist;
                clusterIndex = j;
            }
        }
        if (clusters[clusterIndex]->datapoint == -1){
            clusters[clusterIndex]->datapoint = i;
        }else {
            current = clusters[clusterIndex];
            new = (ELEMENT*)malloc( sizeof(ELEMENT));
            asser(new !=NULL && "An Error Has Occurred");
            new -> datapoint = i;
            new -> next = current;
            clusters[clusterIndex] = new;
        }
    }
}

int updateCentroids(double **cents, LINK *clusters, double** inputMatrix ,int k, int d) {
    int i, j, m, s, difference = 0, sizeOfCluster = 0;
    double epsilon = 0.001;
    LINK current = NULL;
    double *sum = (double *)calloc(d, sizeof(double));
    assert(sum != NULL && "An Error Has Occurred");

    for (i=0; i<k; i++) {
        current = clusters[i];
        /*printf("currently in cluster %d\n", i);*/
        while (current != NULL && current->datapoint != -1){
            for  (j=0; j<d; j++) {
                sum[j] += inputMatrix[current->datapoint][j];
            }
            current = current->next;
            sizeOfCluster++;
        }
        /*printf("cluster %d finised the while loop\n", i);*/

        for (m=0; m<sizeOfCluster; m++){
            sum[m] /= sizeOfCluster;
        }
        /*printf("cluster %d's new centroid is calculated\n", i);*/
        difference += (calculateNorma(cents[i], sum, d) > epsilon)?1:0;
        /*printf("difference is %d\n", difference);*/
        cents[i] = sum;
        for (s=0; s<d; s++) {
            sum[s] = 0;
        }
        printf("updated  cluster %d/%d\n", (i+1), d);
    }
    restartClusters(clusters, k);
    printf("restarted clusters");
    free(sum);
    printf("finished updating centroids");
    return difference;
}

double calculateNorma(double *old, double * new, int d){
    double sum = 0;
    int i;
    for (i=0; i<d; i++){
        sum += pow(old[i]-new[i],2);
    }
    return sqrt(sum);
}

double calculateDistance(double *datapoint, double *centroid, int d){
    double distance = 0;
    int i;
    for (i = 0; i<d; i++){
        distance+= pow(datapoint[i]-centroid[i],2);
    }
    return distance;
}

int isNatural(char *s){
    while (*s) {
        if(isdigit(*s) == 0){
            return 0;
        }
        s++;
    }
    return 1;
}
int main(int argc, char *argv[]){

    /*check validity of inputs- if not "Invalid Input!"*/
    /*IMPORTANT! - we need to make sure that the input is read corrrectly*/
    int d,size, k, max_iter = 200;
    char *inputFilename =NULL;
    char *outputFilename = NULL;
    FILE *ifp=NULL;
    FILE *outfile=NULL;
    double **datapointMatrix, **centroids;
    LINK *clusters;

    assert((argc > 1) && "Invalid Input!");
    assert((argc <= 5) && "Invalid Input!");
    assert((isNatural(argv[1])) && "Invalid Input!");
    k = atoi(argv[1]);
    assert((k>0) && "Invalid Input!"); /*error if k isn't a natural number */
    if (argc==5){
        max_iter = atoi(argv[2]);
        assert((isNatural(argv[2])) && "Invalid Input!");
        assert((max_iter>0) && "Invalid Input!"); /*error if max_iter isn't a natural number */
        inputFilename=argv[3];
        outputFilename=argv[4];
    }
    else if (argc==4){
        inputFilename=argv[2];
        outputFilename=argv[3];
    }
    else {
        assert(!"Invalid Input!");
    }
    printf("k=%d, max_iter=%d\n\n",k, max_iter);
    ifp = fopen(inputFilename, "r");
    assert(ifp!=NULL);
    printf("Analizing input file...\n");
    d = findDimension(ifp);
    size = findInputSize(ifp);
    assert(size>=k && "Invalid Input!");
    
    printf("Converting the input file into a %dx%d matrix\n",size,d);
    datapointMatrix = fileToDataPoints(ifp, d, size);
    printf("Closing input file\n");
    fclose(ifp);
    
    printf("Inititializing centroids\n");
    centroids = inizializeCentroids(k, d, datapointMatrix); /*this array holds K datapoints*/
    clusters = (LINK *)calloc(k, sizeof(LINK));
    assert(clusters != NULL && "An Error Has Occurred");
    printf("Calculating clusters\n");
    restartClusters(clusters, k);
    kMeans(k, size, d, max_iter, centroids, clusters, datapointMatrix);
    printf("Successfully ran kMeans\n");
    printf("First index of the first centroid is %f\n",centroids[0][0]);
    printf("Opening output file\n");
    printf("Writing to output file...\n");
    outfile = fopen(outputFilename,"w");
    assert(outfile != NULL);
    fclose(outfile);
    printf("Closing output file\n");
    printf("Done!\n");
    return 0;
}