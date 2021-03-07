import numpy as np
import os
import skimage.morphology
from operator import itemgetter as it
from skimage.filters import threshold_otsu


def condition(x,y,history,labeled):

    if 0<x<np.shape(labeled)[0] and 0<y<np.shape(labeled)[1]:
        if history[x,y]==0:
            return True
        else:
            return False
    else:
        return False

def shuffling(image, output_folder, filename, out='Full'):

#    image = np.load(os.path.join(input_folder, filename+'.npy'))

    
    image = np.nan_to_num(image)
    thresh = threshold_otsu(image) 
    mean_image = np.mean(image)
    
    # Defining binary image to identify spots in pattern. Threshold depends on images. Needs to be calculated separately for each image
    # change in line 93 as well
    binary_image = image > thresh + 0.5*mean_image #@300000
#    binary_image = image > mean_image #@0
    
    labeled_image = skimage.morphology.label(binary_image, connectivity = 2, background=thresh)
    
    max_ = np.max(labeled_image) 
    spot_indices = []
    bkg_indices = []
    for i in range(max_):
        region = np.argwhere(labeled_image==i)
        if len(region)>5000:
            bkg_indices.append(region)
        else:
            spot_indices.append(region)
    
    history_image = np.zeros(np.shape(labeled_image))
    shuffled = np.zeros(np.shape(image))



    for j in range(1,len(spot_indices)):
        spot = spot_indices[j]
        l = len(spot)
        relative_position = []
        for k in range(l):
            relative_position.append(spot[k]-spot[0])
           
        new_spot_position = []   
        while len(new_spot_position) == 0:
            random_place = [np.random.randint(0,np.shape(image)[0]), np.random.randint(0,np.shape(image)[1])]
            history_temp = history_image.copy()
            while condition(random_place[0],random_place[1],history_image,labeled_image) == True:
                new_spot_position = [random_place]
                history_image[random_place[0],random_place[1]] = 1; shuffled[random_place[0],random_place[1]] = image[spot[0][0],spot[0][1]] 
                for m in range(1,len(relative_position)):
                    temp = random_place+relative_position[m]
                    new_spot_position.append(temp)
                    if condition(temp[0],temp[1],history_image,labeled_image) == True:
                        continue
                    else:
                        new_spot_position = []
                        history_image = history_temp
                        break
        for n in range(len(new_spot_position)):
            temp_n = new_spot_position[n]
            history_image[temp_n[0],temp_n[1]] =1;shuffled[temp_n[0],temp_n[1]] = image[spot[n][0],spot[n][1]]

    for p in range(len(bkg_indices)):
        rest_ind_img = bkg_indices[p]
        rest_ind_shuff = np.argwhere(shuffled==0)
        history_new = np.zeros(np.shape(shuffled))    
        ind = 0
        while len(rest_ind_shuff)>=1:
            if len(rest_ind_img)==0:
                break
            x = rest_ind_img[ind][0]; y = rest_ind_img[ind][1]
            n = np.random.randint(0,len(rest_ind_shuff))
            if history_new[rest_ind_shuff[n][0],rest_ind_shuff[n][1]]==0:
        #        print(n)
                shuffled[rest_ind_shuff[n][0],rest_ind_shuff[n][1]] = image[x,y]
                history_new[rest_ind_shuff[n][0],rest_ind_shuff[n][1]]=1
                rest_ind_img = np.delete(rest_ind_img,ind, axis=0)
                rest_ind_shuff = np.delete(rest_ind_shuff,n, axis=0)

    binary_shuffled = shuffled > thresh + 0.5*mean_image  
#    binary_shuffled = shuffled > mean_image 
    
    print(np.sum(image));print(np.sum(shuffled))
    if out=='Full':
        return binary_image,binary_shuffled, shuffled, thresh, mean_image
    else:
        return shuffled

