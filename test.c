#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/*int main(int argc, char *argv[]){
    char *inputFilename = argv[1];
   
    FILE *ifp = NULL;
    ifp = fopen(inputFilename, "r");
    assert(ifp != NULL);
    double num;
    while (fscanf(ifp, "%lf %*[,] ", &num ) > 0) {
        printf("%lf",num);
        printf(",");
    }
    fclose(ifp); 
    printf("%s", "done");
}*/
typedef struct list ELEMENT;
typedef ELEMENT *LINK;
typedef double *datapoint; 
/* LINK is a pointer to ELEMENT */
/* datapoint is a pointer to a double array */
void kMeans(int k, int max_iter);
void inizializeCentroids(int k, LINK head, double **cents);
int inizializeDataPoints(FILE *ifp);
void calcClosestCluster(double x_i,double cents);
double** fileToDataPoints(FILE *ifp,int d, int size);
int findDimension(FILE *ifp);

struct list{          /*item in a list*/
    double *datapoint; 
    struct list *next;   
};

double ** fileToDataPoints(FILE *ifp ,int d, int size){
    /*  reading the file into list of dataPoints (each dataPoint is an array dynamicly alocated- of side d
    d is the dimension of the given vectors int the input file*/
    double currentCoordinate;
    double **matrix = calloc(size, sizeof(double *)); 
    double *datapoint = calloc(d, sizeof(double));
    if (fscanf(ifp, "%lf %*[,] ", &currentCoordinate ) > 0) {
        head = ( ELEMENT* )malloc( sizeof( ELEMENT ) );
        assert(head!=NULL);
        /*creating the inner array */
        double *vector = (double *)malloc(d*sizeof( double ));
        vector[0]=currentCoordinate;
        for(int i=1;i<d;i++){
            vector[i]=currentCoordinate;
            fscanf(ifp, "%lf %*[,] ", &currentCoordinate );
        }
        head->datapoint = vector;
        tail = head;
        while (fscanf(ifp, "%lf %*[,] ", &currentCoordinate ) > 0){
            tail->next= ( ELEMENT* )malloc( sizeof( ELEMENT ) );
            tail = tail->next;
            assert(tail!=NULL);
            double *vector = (double *)malloc(d*sizeof( double ));
            vector[0]=currentCoordinate;
            for(int i=1;i<d;i++){
                vector[i] = currentCoordinate;
                fscanf(ifp, "%lf %*[,] ", &currentCoordinate);
            }
            tail->datapoint = vector;
            tail->next=NULL;
        }
    }
    return head;
}

int findDimension(FILE *ifp) {
    char c;
    int commacnt=0;
    int dim=0;
    do {
        c = fgetc(ifp);
        /*printf("%c",c)*/;
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
void inizializeCentroids(int k, LINK head, double **cents){
        /*if there are less vectoer in the input than k than print "invalid input" and terminate*/
        LINK current = head;
        int d = sizeof(head->datapoint)/sizeof((head->datapoint)[0]);   
        for(int i=0; i<k; i++){
            double *vector = cents[i];
            for (int j=0; j<d; j++) { 
                vector[j]=current->datapoint[j];
            }
            current = current->next;        
        }
}

void kMeans(int k, int max_iter){
    int iter = 0;
    int continue_condition = 1;
    while (iter<max_iter && continue_condition){
        /*
        1) for each item in the input list
            2) for each centroid
                3) calculate distance 
                4) update minimum distance and cluster _i_ from which the distance is minimal
            5) add item to the list in cluster _i_

        6) update centroids + empty cluster lists and free the space
        7) update iter
        8) update continue_codition
        */
    }
}

/*void assignToClusters(LINK head, double **cents,LINK **clusters, int K) {
    int distance;
    int tempDist;
    int clusterIndex;
    /*
        1) for each item in the input list
            2) for each centroid
                3) calculate distance 
                4) update minimum distance and cluster _i_ from which the distance is minimal
            5) add item to the list in cluster _i_
        */
/*    LINK current = head;
    while (current != NULL) {
        for (int i=0; i<K; i++){
            tempDist = calculteDistance(current->datapoint,cents[i]);
            if (tempDist<distance) {
                distance = tempDist;
                clusterIndex = i;
            }
        }
        clusters[clusterIndex]->next = current;
        current = current->next;
    }

}
*/

int calculateDistance(double *datapoint, double *centroid, int d){
    int distance = 0;
    for (int i = 0; i<d; i++){
        distance+= pow(datapoint[i]-centroid[i],2);
    }
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
    printf("opened input file successfuly\n");
    int d = findDimension(ifp);
    printf("d = %d\n",d);
    int size = findInputSize(ifp);
    printf("size = %d\n",size);

    double** datapointMatrix = fileToDataPoints(ifp, d, size);
    fclose(ifp);
    printf("closed input file successfully\n");
    
    double **centroids = malloc(K*sizeof(datapoint)); //this array holds K datapoints
    LINK **clusters = malloc( K*sizeof(ELEMENT));
    inizializeCentroids(K, headOfList, centroids);

    void kMeans(int K, int max_iter);
    //double item = (*(centroids)[0])->datapoint[0];
    printf("%d\n",centroids[0]);
    printf("success!\n");
    outfile = fopen(outputFilename,"w");
    assert(outfile != NULL);
    printf("opened outfile successfully\n");
    fclose(outfile);
    printf("closed outfile successfully\n");
    return 0;
}