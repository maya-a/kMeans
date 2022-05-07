import sys 
import math 

def an_error_has_occurred():
    sys.exit("An Error Has Occurred")

def files_to_dataframe(file_name_1,file_name_2):
    file1 = pd.read_csv(file_name_1)
    size1=file1.shape[1]
    file2 = pd.read_csv(file_name_1)
    size = [str(i) for i in range(file1.shape[1])]
    file1 = pd.read_csv(file_name_1, names=size)
    size = [str(i+size1-1) for i in range(file2.shape[1])]
    size[0] = str(0)
    file2 = pd.read_csv(file_name_2, names=size)

    data = pd.merge(file1, file2, on ='0')
    data=data.drop(['0'],axis=1)
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
            for idx in current_cluster:
                for j in range(d):
                    sum_of_vectors[j]+=input_mat[idx][j] 
            sum_of_vectors=[float(sum_of_vectors[m]/size_of_cluster) for m in range(d)]
            diff += 1 if calculate_norma(centroids[i],sum_of_vectors)>epsilon else 0
            centroids[i]=sum_of_vectors[:]
           
        return diff>0
    except (ZeroDivisionError): #TypeError
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
        clusters=assign_to_clusters(input_mat,centroids,k)
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
    try:
        input_file=None
        output_file=None
        check_is_natural(sys.argv[1])
        k=int(sys.argv[1])
        max_iter=200
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
        if(k > len(data_point_matrix)):
            invalid_input()
        centroids=kmeans(data_point_matrix,k, max_iter)
        write_to_file(output_file,centroids)
    except Exception as e:
        print("An Error Has Occurred")
        exit()

main()