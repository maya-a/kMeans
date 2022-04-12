#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

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

struct list{          /*this is an item in a list*/
    double *datapoint; 
    struct list *next;   
};

double ** fileToDataPoints(FILE *ifp, int d, int size) {
    double currentCoordinate;
    double *vector = (double *)calloc(size*d, sizeof(double));
    assert(vector != NULL);
    printf("allocated big contiuous vector\n");
    double **matrix = (double **)calloc(size, sizeof(double *)); 
    assert(matrix != NULL);
    printf("allocated pointer vector\n");
    for (int i=0; i<size; i++) {
        matrix[i] = vector + i*d;
    }
    printf("initialized pointers\n");
    for (int i=0; i<size; i++){  
        for(int j=0; j<d; j++){
            fscanf(ifp, "%lf %*[,]", &currentCoordinate);
            matrix[i][j] = currentCoordinate;
        }
    }
    printf("method finished\n");
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
    /*printf("%lf",datapointMatrix[499][2])*/;
    fclose(ifp);
    printf("closed input file successfully\n");
    
    double **centroids = (double **)calloc(K,sizeof(datapoint)); //this array holds K datapoints
    assert(centroids != NULL);
    LINK **clusters = (LINK **)calloc(K,sizeof(ELEMENT));
    assert(clusters != NULL);
    //inizializeCentroids(K, headOfList, centroids);

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