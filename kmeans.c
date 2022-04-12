#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


void kMeans(int k,const char filename, int max_iter);
void inizializeCentroid(int k,FILE *ifp, double *cents);
int inizializeDataPoints(FILE *ifp);
void calcClosestCluster(double x_i,double cents);
LINK fileTolistOfDataPoints(FILE *ifp,int d);

struct list{
    double *datapoint; 
    struct list *next;   
};

typedef struct list ELEMENT;
typedef ELEMENT* LINK;

int main(int argc, char *argv[]){
    //check validity of inputs- if not "Invalid Input!"
    
    int K=argv[1];
    char *inputFilename=argv[2];
    char *outputFilename=argv[3];
    int max_iter=200; 
    if (argc==5){
        max_iter=argv[5];
    }
    FILE *ifp=NULL;
    FILE *outfile=NULL;
    int N=200; //change later 
    double centriods[K];
    char *nameofInput;
    char *nameofOutput;
    
    size_t lenIn = strlen(inputFilename);
    size_t lenOut=strlen(outputFilename);
    // strndup returns a pointer to an array that it created
    nameofInput = strndup(inputFilename, lenIn >= 4 ? lenIn - 4 : 0);
    nameofOutput = strndup(outputFilename, lenOut >= 4 ? lenOut - 4 : 0);

    ifp=fopen(nameofInput, "r");
    assert(ifp!=NULL);
    outfile=fopen(nameofOutput,"w");
    assert(outfile!=NULL);
    return 0;
}

void kMeans(int K,const char filename, int max_iter){


inizializeCentroid(K,ifp,centriods);
double x_i=0; 
for(int i=1 ; i<=N ; i++) {
x_i = fscanf(); //need to read from file the double number 
//calculate euclidian distances from each centriod
//add x_i to the clossest cluster- i thought of creating list for the clusters, what do you think? 
// how do we implement methods for the list? it in possible? 

}
fclose(ifp);
}

LINK fileTolistOfDataPoints(FILE *ifp ,int d){
    /*  reading the file into list of dataPoints (each dataPoint is an array dynamicly alocated- of side d
    d is the dimension of the given vectors int the input file*/
    double currentCoordinate; 
    LINK head = NULL, tail = NULL;
    if (fscanf(ifp, "%lf %*[,] ", &currentCoordinate ) > 0){
        head = ( ELEMENT* )malloc( sizeof( ELEMENT ) );
        assert(head!=NULL);
        /*creating the inner array */
        double* vector = (double*)malloc(d*sizeof( double ));
        for(int i=0;i<d;i++){
            vector[i]=currentCoordinate;
            fscanf(ifp, "%lf %*[,] ", &currentCoordinate );
        }
        head->datapoint = vector;
        tail = head;
        while (fscanf(ifp, "%lf %*[,] ", &currentCoordinate ) > 0){
            tail->next= ( ELEMENT* )malloc( sizeof( ELEMENT ) );
            tail = tail->next;
            assert(tail!=NULL);
            tail->datapoint = vector;
            double* vector = (double*)malloc(d*sizeof( double ));
            for(int i=0;i<d;i++){
                vector[i]=currentCoordinate;
                fscanf(ifp, "%lf %*[,] ", &currentCoordinate );
            }
            tail ->next=NULL;
        }
        
    }
    return head;

}

void inizializeCentroid(int K,FILE *ifp, double *cents){

for( int i=0;i<K;i++){
    fscanf(ifp,"%e",&cents[i]);
}
}


int inizializeDataPoints(FILE *ifp){
    // read first line first then rewind 
    char c = fgetc(ifp);
    int commacnt=0;
    int dim=0;
    do
    {
        // Taking input single character at a time
        if( feof(ifp) ){
            return -1; 
        }

        if (c ==10) // c==\n 
            {
                dim=commacnt+1;
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
    double vector;
    while (fscanf(ifp, "%lf %*[,] ", &vector ) > 0){

    }
    fclose(ifp);
   
}

}
}


