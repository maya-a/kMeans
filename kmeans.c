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
void writeToFile(double **cents, int k, int d, char * file_name);
void validityCheck1(int argc, char *argv[]);
void validityCheck2(int k, int max_iter, int size);

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
    int cnt = 0;
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
    int i,j;
    centroid = (double *)calloc(k*d, sizeof(double));
    assert(centroid != NULL && "An Error Has Occurred");
    matrix = (double **)calloc(k, sizeof(double *)); 
    assert(matrix != NULL && "An Error Has Occurred");
    for (i=0; i<k; i++) {
        matrix[i] = centroid + i*d;
    }
    for (i=0; i<k; i++){
        /*consider using memcpy(matrix[i], datapoints[i], (d+1)* sizeof(double));*/
        for (j=0; j<d; j++){
            matrix[i][j] = datapoints[i][j]; 
        }
    }  
    return matrix; 
}

void restartClusters(LINK *clusters, int k) {
    int i;
    LINK current, next;
    for (i=0; i<k; i++) {
        current = clusters[i];
        /*delete_list(current);*/
        while (current != NULL) {
            next = current->next;
            free(current);
            current = next;
        }
        printf("nullified  cluster %d/%d\n", i+1,k);
        clusters[i] = (ELEMENT*)malloc(sizeof(ELEMENT));
        assert(clusters[i] != NULL && "An Error Has Occurred");
        clusters[i]->datapoint = -1;
        clusters[i]->next = NULL;
    }
}

void delete_list(LINK head) {
    if (head != NULL){
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
       printf("iteration #%d\nassigning datapoints to clusters", iter+1);
       assignToCluster(cents, matrix, clusters, k, size, d);
       printf("updating centroids\n");
       continue_condition = updateCentroids(cents, clusters, matrix, k, d);
       iter++;
    }
}

void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d) {
    double distance;
    double tempDist = 0;
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
        distance = DBL_MAX;
        for (j=0; j<k; j++){
            tempDist = calculateDistance(datapoints[i], cents[j], d);
            if (tempDist < distance) {
                distance = tempDist;
                clusterIndex = j;
            }
        }
        /*printf("\nx%d is closest to centroid %d\n\n",i, clusterIndex);*/
        if (clusters[clusterIndex]->datapoint == -1){
            clusters[clusterIndex]->datapoint = i;
        }else {
            current = clusters[clusterIndex];
            new = (ELEMENT*)malloc( sizeof(ELEMENT));
            assert(new !=NULL && "An Error Has Occurred");
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
        while (current != NULL && current->datapoint != -1){
            /*printf("currently on datapoint %d\n",current->datapoint);*/
            for  (j=0; j<d; j++) {
                sum[j] += inputMatrix[current->datapoint][j];
                /*printf("%f",sum[j]);*/
            }
            /*printf("\n");*/
            current = current->next;
            sizeOfCluster++;
        }
        /*printf("cluster %d finised the while loop\n", i);*/

        for (m=0; m<d; m++){
            sum[m] /= sizeOfCluster;
        }
        /*printf("cluster %d's new centroid is calculated\n", i);*/
        difference += (calculateNorma(cents[i], sum, d) > epsilon)?1:0;
        /*printf("difference is %d\n", difference);*/
        /*consider using memcpy(cents[i],sum, d*sizeof(double));*/
        printf("new centroid %d:(",i+1);
        for (s=0; s<d; s++) {
            cents[i][s] = sum[s];
            sum[s] = 0;
            printf("%f,",cents[i][s]);
        }
        printf(")\n");
        printf("\n*******updated centroid %d/%d*******\n", (i+1), d);
    }
    restartClusters(clusters, k);
    printf("re-initialized clusters\n");
    free(sum);
    printf("centroids have been updated\n");
    /*printf("continue is %d\n",difference > 0);*/
    return difference > 0;
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
    int j;
    for (j = 0; j<d; j++){
        distance+= pow(datapoint[j] - centroid[j],2);
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

void writeToFile(double **cents, int k, int d, char * fileName){
    int i,j;
    FILE *ofp = NULL;
    ofp = fopen(fileName,"w");
    assert(ofp != NULL && "Invalid Input!");
    for (i = 0; i<k; i++){
        for (j=0; j<d; j++){
             fprintf(ofp,"%.4f,",cents[i][j]);
        }
        fprintf(ofp,"\n");
    }
    fclose(ofp);
    printf("Closing output file\n");
}

void validityCheck1(int argc, char *argv[]){
    if ((argc != 4 && argc != 5) || 
        !isNatural(argv[1]) || 
        (argc == 5 && !isNatural(argv[2]))){
        printf("Invalid Input!");
        exit(1);
        }
    else{
        
    }
}
void validityCheck2(int k, int max_iter, int size){
    if (k == 0 ||
        max_iter == 0 ||
        k>size){
        printf("Invalid Input!");
        exit(1);
        }
}

int main(int argc, char *argv[]){
    /*check validity of inputs- if not "Invalid Input!"*/
    /*IMPORTANT! - we need to make sure that the input is read corrrectly*/
    int d, size, k, max_iter = 200;
    char *inputFileName = NULL;
    char *outputFileName = NULL;
    FILE *ifp=NULL;
    double **datapointMatrix, **centroids;
    LINK *clusters;
    validityCheck1(argc,argv);
    k = atoi(argv[1]);

    if (argc==5){
        max_iter = atoi(argv[2]);
        inputFileName=argv[3];
        outputFileName=argv[4];
    }
    else if (argc==4){
        inputFileName=argv[2];
        outputFileName=argv[3];
    }
    /*else {
        assert(!"Invalid Input!");
    }*/
    printf("k=%d, max_iter=%d\n\n",k, max_iter);
    ifp = fopen(inputFileName, "r");
    assert(ifp != NULL && "Invalid Input!");
    printf("Analizing input file...\n");
    d = findDimension(ifp);
    size = findInputSize(ifp);
    validityCheck2(k, max_iter, size);
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
    writeToFile(centroids, k, d, outputFileName);
    free(datapointMatrix);
    free(centroids);
    free(clusters);
    printf("Done!\n");
    return 0;
}