from multiprocessing import Process, Pool
from multiprocessing.managers import SharedMemoryManager
from multiprocessing.shared_memory import SharedMemory
from typing import Tuple
import numpy as np
import cv2
import os 
import pandas as pd
import matplotlib.pyplot as plt
import tifffile
import cv2
from PIL import ImageColor
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool, Value, Array
from time import time
import skimage.segmentation as seg

def top_clusters(centroids, n, cluster_id= "cluster_id", sort_val= "nCluster"):
    clusters= centroids[[cluster_id, sort_val]].drop_duplicates()
    clusters= clusters.iloc[np.argsort(clusters[sort_val])]
    top_clusters= np.array(clusters[cluster_id].tail(n= n))
    return(top_clusters)

def create_np_array_from_shared_mem() -> np.ndarray:
    #shared_mem: SharedMemory, shared_data_dtype: np.dtype, shared_data_shape: Tuple[int, ...]
    arr = np.frombuffer(shared_mem.buf, dtype=SHARED_DATA_DTYPE)
    arr = arr.reshape(SHARED_DATA_SHAPE)
    return arr

def child_process(row):
    currtime= time()
    arr= create_np_array_from_shared_mem()
     # modify the array backed by shared memory
    a= cluster_centroids.Y_centroid.iloc[row]
    b= cluster_centroids.X_centroid.iloc[row]
    arr[a,b,:] = 0
    cv2.floodFill(image= arr, 
                  mask= None, 
                  seedPoint=(a,b), 
                  newVal= rgb_color,
                  flags= cv2.FLOODFILL_FIXED_RANGE 
                  )
        ## This doesn't do the flood fill properly.
        ## I can either try the scipy segmentation floodfill or the 
   

    print("Done with row in ", time() - currtime)

def main():

    path= "/stor/work/Ehrlich/Users/John/projects/AKOYA"
    data_to_share= cv2.imread(os.path.join(path,"scripts/cell_segmentation/MESMER/data/multichannel_thymus_nuc_CD31_CD49f_MHCII/segmentation/fov1.tiff"))
    centroids= pd.read_csv(os.path.join(path,"data/AKOYA_1mo_XY_supGrant_xshift_clusters.csv"))
    centroids= centroids[(centroids.Y_centroid < 6500)  & (centroids.Y_centroid > 6250) &
                         (centroids.X_centroid < 13050) & (centroids.X_centroid > 12750)] 
    
    hex_colors= ["#1B9E77", "#528B54", "#897932", "#C16610", "#C8611F", "#AB6653",
                 "#8D6B86", "#796DB1", "#9B58A5", "#BC4399", "#DD2E8D", "#CC4373",
                 "#A66753", "#808B34", "#70A61B", "#96A713", "#BBA90B", "#E0AA03",
                 "#D59D08", "#C38E10", "#B07E18", "#9D7426", "#8B6F3B", "#786A50",
                 "#666666"] 
    
    global SHARED_DATA_DTYPE
    global SHARED_DATA_SHAPE
    global SHARED_DATA_NBYTES
    
    SHARED_DATA_DTYPE  = data_to_share.dtype
    SHARED_DATA_SHAPE  = data_to_share.shape
    SHARED_DATA_NBYTES = data_to_share.nbytes

    with SharedMemoryManager() as smm:
        global shared_mem
        shared_mem = smm.SharedMemory(size=SHARED_DATA_NBYTES)

        arr = create_np_array_from_shared_mem()
        arr[:] = data_to_share  # load the data into shared memory
        rgb_colors= [ImageColor.getcolor(color, "RGB") for color in hex_colors]

        clusters= top_clusters(centroids= centroids, n= 20)
    
        
        for i in range(len(clusters)):
            currtime= time()
            cluster   = clusters[i]
            print(cluster)
            global cluster_centroids
            cluster_centroids= centroids[(centroids.cluster_id == cluster)]
            print(cluster_centroids.shape)
            global rgb_color
            rgb_color = rgb_colors[i]
            print(rgb_color)
            Pool(processes=30).map(child_process, range(cluster_centroids.shape[0]))

            print(f"DONE WITH CLUSTER {cluster} after {time() - currtime}")

        out_path= os.path.join(path, "data/AKOYA_1mo_thymus_Xshift_cluster_overlay_pool_take2_fixed_range_mop.tiff")
        tifffile.imwrite(out_path, arr)
        del arr  # delete np array so the shared memory can be deallocated

main()

