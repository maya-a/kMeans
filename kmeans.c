#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>

typedef struct list ELEMENT;
typedef ELEMENT *LINK;
typedef double *datapoint; 
/* LINK is a pointer to ELEMENT */
/* datapoint is a pointer to a double array */

void anErrorHasOccurred();
void invalidInput();
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

struct list { 
    int datapoint; 
    struct list *next;   
};

void anErrorHasOccurred(){
    printf("An Error Has Occurred");
    exit(1);
}

void invalidInput(){
    printf("Invalid Input!");
    exit(1);
}

double ** fileToDataPoints(FILE *ifp, int d, int size) {
    double currentCoordinate;
    double *vector = NULL;
    double **matrix = NULL;
    int i, j;
    
    vector = (double *)calloc(size*d, sizeof(double));
    if (vector == NULL){
        anErrorHasOccurred();
    }
    matrix = (double **)calloc(size, sizeof(double *)); 
    if (matrix == NULL){
        anErrorHasOccurred();
    }
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
    double *centroid = NULL;
    double **matrix = NULL;
    int i,j;
    centroid = (double *)calloc(k*d, sizeof(double));
    if (centroid == NULL){
        anErrorHasOccurred();
    }
    matrix = (double **)calloc(k, sizeof(double *)); 
    if (matrix == NULL){
        anErrorHasOccurred();
    }
    for (i=0; i<k; i++) {
        matrix[i] = centroid + i*d;
    }
    for (i=0; i<k; i++){
        for (j=0; j<d; j++){
            matrix[i][j] = datapoints[i][j]; 
        }
    }  
    return matrix; 
}

void restartClusters(LINK *clusters, int k) {
    int i;
    LINK current;
    for (i=0; i<k; i++) {
        current = clusters[i];
        delete_list(current);
        clusters[i] = (ELEMENT*)malloc(sizeof(ELEMENT));
        if (clusters[i] == NULL){
            anErrorHasOccurred();
        }
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
    
       assignToCluster(cents, matrix, clusters, k, size, d);
       continue_condition = updateCentroids(cents, clusters, matrix, k, d);
       iter++;
    }
}

void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d) {
    double distance;
    double tempDist = 0;
    int clusterIndex, i, j;
    LINK current, new;
 
    for(i=0; i<size; i++) {
        clusterIndex = -1;
        distance = DBL_MAX;
        for (j=0; j<k; j++){
            tempDist = calculateDistance(datapoints[i], cents[j], d);
            if (tempDist < distance) {
                distance = tempDist;
                clusterIndex = j;
            }
        }
        if (clusters[clusterIndex]->datapoint == -1) {
            clusters[clusterIndex]->datapoint = i;
        }else {
            current = clusters[clusterIndex];
            new = (ELEMENT*)malloc( sizeof(ELEMENT));
            if (new == NULL){
                anErrorHasOccurred();
            }
            new -> datapoint = i;
            new -> next = current;
            clusters[clusterIndex] = new;
        }
    }
}

int updateCentroids(double **cents, LINK *clusters, double** inputMatrix ,int k, int d) {
    int i, j, m, s, difference = 0, sizeOfCluster;
    double epsilon = 0.001;
    LINK current = NULL;
    double *sum = (double *)calloc(d, sizeof(double));
    if (sum == NULL){
        anErrorHasOccurred();
    }
    for (i=0; i<k; i++) {
        current = clusters[i];
        sizeOfCluster = 0;
        while (current != NULL && current->datapoint != -1){
            for  (j=0; j<d; j++) {
                sum[j] += inputMatrix[current->datapoint][j];
            }

            current = current->next;
            sizeOfCluster++;
        }
        if (sizeOfCluster == 0){
            anErrorHasOccurred();
        }

        for (m=0; m<d; m++){
            sum[m] /= sizeOfCluster;
        }
        difference += (calculateNorma(cents[i], sum, d) > epsilon)?1:0;
        for (s=0; s<d; s++) {
            cents[i][s] = sum[s];
            sum[s] = 0;
        }
    }
    restartClusters(clusters, k);
    free(sum);
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
    if (ofp == NULL){
        invalidInput();
    }
    for (i = 0; i<k; i++){
        for (j=0; j<d-1; j++){
            fprintf(ofp,"%.4f,",cents[i][j]);
        }
        fprintf(ofp,"%.4f",cents[i][j]);
        fprintf(ofp,"\n");
    }
    fclose(ofp);
}

void validityCheck1(int argc, char *argv[]){
    if ((argc != 4 && argc != 5) || 
        !isNatural(argv[1]) || 
        (argc == 5 && !isNatural(argv[2]))){
        printf("Invalid Input!");
        exit(1);
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
    
    ifp = fopen(inputFileName, "r");
    if (ifp == NULL){
        invalidInput();
    }
    d = findDimension(ifp);
    size = findInputSize(ifp);
    validityCheck2(k, max_iter, size);
    datapointMatrix = fileToDataPoints(ifp, d, size);
    fclose(ifp);
    
    centroids = inizializeCentroids(k, d, datapointMatrix); /*this array holds K datapoints*/
    clusters = (LINK *)calloc(k, sizeof(LINK));
    if (clusters == NULL){
        anErrorHasOccurred();
    }
    restartClusters(clusters, k);
    kMeans(k, size, d, max_iter, centroids, clusters, datapointMatrix);
    writeToFile(centroids, k, d, outputFileName);
    free(datapointMatrix);
    free(centroids);
    free(clusters);
    return 0;
}