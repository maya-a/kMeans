import sys 
import math 

""" TODO:  update centroids-done
write to the outputfile--done 
take care of continue condition--done 
check validity of input- continue+wait for response of teacher
 """
def an_error_has_occurred():
    sys.exit("An Error Has Occurred")


def file_to_matrix(input_file):
    input_mat = []
    try: 
        f = open(input_file, "r")
        for line in f:
            currentline = line.split(",")
            last = currentline[-1][:-1]
            currentline = currentline[:-1]
            currentline.append(last)
            currentline = [float(i) for i in currentline]
            input_mat.append(currentline)
        #input_mat.pop()
        f.close() 
        return input_mat   
        #print(input_mat)
    except IOError:
        an_error_has_occurred() 

def create_centroids(k, input_mat):
    centroids=[]
    centroids=[input_mat[i] for i in range(k)]
    return centroids

def calculate_distance(centroid, data_point):
    return sum([pow(centroid[i]-data_point[i],2) for i in range(len(centroid))])

def calculate_norma(old_centroid, new_centroid):
    return math.sqrt(sum([pow(old_centroid[i]-new_centroid[i],2) for i in range(len(old_centroid))]))

def assign_to_clusters(input_mat,centroids,k):
    clusters=[[]]*k
    """  
        1) for each item in the input list
            2) for each centroid
                3) calculate distance 
                4) update minimum distance and cluster _i_ from which the distance is minimal
            5) add item to the list in cluster _i_
           """
    #infinity = float("inf")
    #sys.float_info.max

    for idx,data_point in enumerate(input_mat):
        distance= float("inf")
        temp_dist=0
        cluster_index=0
        for i, centroid in enumerate(centroids):
            temp_dist=calculate_distance(centroid, data_point)
            if (temp_dist<distance):
                distance = temp_dist
                cluster_index = i
        if (clusters[cluster_index] == []):
            clusters[cluster_index]=[idx]
        else: 
            clusters[cluster_index].append(idx)
    return clusters


def update_centroids(input_mat,centroids,clusters):
    try:
        d=len(input_mat[0])
        epsilon=0.001
        diff=0
        for i ,current_cluster in enumerate(clusters):
            sum_of_vectors=[0]*d
            size_of_cluster = len(current_cluster)
            #print(current_cluster)
            #print("cluster " + str(i+1) +"'s size is " + str(size_of_cluster))
            #print("cluster "+str(i+1) + "'s datapoints: ")
            for idx in current_cluster:
                for j in range(d):
                    #print(str(input_mat[idx][j]))
                    sum_of_vectors[j]+=input_mat[idx][j] 
            #print("---")
            #if (size_of_cluster == 0):
                #print("cluster " + str(i+1) + " is empty!")
            sum_of_vectors=[float(sum_of_vectors[m]/size_of_cluster) for m in range(d)]
            diff += 1 if calculate_norma(centroids[i],sum_of_vectors)>epsilon else 0
            centroids[i]=sum_of_vectors[:]
           
        return diff>0
    except (ZeroDivisionError): #TypeError
        #print("error in update_centroids")
        an_error_has_occurred()

def write_to_file(file_name, centroids):
    try:
        f=open(file_name,"w")
        for centroid in centroids:
            res=""
            for indx in (range(len(centroid)-1)):
                res+="{:.4f}".format(centroid[indx])+","
            res+="{:.4f}".format(centroid[-1])+"\n"
            f.write(res)
        f.close()
    except IOError: 
        an_error_has_occurred()

def kmeans(input_mat , k ,max_iter):
    
    centroids=create_centroids(k,input_mat)
    iter=0
    next_iteration=True
    while(iter<max_iter and next_iteration):
        """
        1) for each item in the input list
            2) for each centroid
                3) calculate distance 
                4) update minimum distance and cluster _i_ from which the distance is minimal
            5) add item to the list in cluster _i_

        6) update centroids 
        7) update iter
        8) update continue_codition """
        
        clusters=assign_to_clusters(input_mat,centroids,k)
        #print("updating centroids after " +str(iter+1)+ " assignment of datapoints")
        next_iteration=update_centroids(input_mat,centroids,clusters)
        iter+=1
    return centroids

def invalid_input():
    sys.exit("Invalid Input!")

def check_is_natural(num):
    try:
        # Convert it into float
        
        val_f = float(num)
        val_int=int(float(num))
    
    except ValueError:
        invalid_input()
    
    if val_f!=val_int or val_int<=0:
        invalid_input()



def main(): 
    input_file=None
    output_file=None
    check_is_natural(sys.argv[1])
    k=int(sys.argv[1])
    max_iter=200
    #take care of type cheching- is it int? 
    if(len(sys.argv)==4):
        input_file=sys.argv[2]
        output_file=sys.argv[3]
        
    elif (len(sys.argv)==5):
        check_is_natural(sys.argv[2])
        max_iter=int(sys.argv[2])
        input_file=sys.argv[3]
        output_file=sys.argv[4]
    else: 
        invalid_input()
   
    #Reading from file
    data_point_matrix=file_to_matrix(input_file)
    if(k>len(data_point_matrix)):
        invalid_input()
    centroids=kmeans(data_point_matrix,k, max_iter)
    write_to_file(output_file,centroids)
    


main()