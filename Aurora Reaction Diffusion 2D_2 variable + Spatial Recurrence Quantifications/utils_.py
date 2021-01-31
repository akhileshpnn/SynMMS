import numpy as np

def removeNans(x, y):
    temp_y = []
    temp_x = []
    temp_idxs = []
    last_non_nan_idx = -1
    i = 0
    max_size_nans = 5
    if np.isnan(y[0]):
        y[0]=np.nanmean(y)
    if np.isnan(y[-1]):
        y[-1] = np.nanmean(y)
    while i < len(y):
        if np.isnan(y[i]):
            j = 0
            while j+i < len(y) and np.isnan(y[j+i]):
                j = j + 1 
            if last_non_nan_idx >= 0 and j == 1:
                temp_idxs.append(last_non_nan_idx + j)
                temp_x.append(x[last_non_nan_idx + j])
                temp_y.append(y[last_non_nan_idx] )
                i = i+1
            elif last_non_nan_idx >= 0 and j <= max_size_nans:
                l = j
                dv = y[i+l] - y[last_non_nan_idx]
                step = dv * 1.0 / l
                for k in range(l):
                    temp_idxs.append(last_non_nan_idx + k)
                    temp_x.append(x[last_non_nan_idx + k])
                    temp_y.append(temp_y[last_non_nan_idx + k] + k * step)
                i = i+l
            elif last_non_nan_idx >= 0 and j > max_size_nans:
                l = j
                for k in range(l):
                    temp_x.append(x[i+k])
                    temp_y.append(np.nanmean(y))
                i = i+l
#        if i == len(y):
#            break
        else:    
            temp_y.append(y[i])
            temp_x.append(x[i]) 
            temp_idxs.append(i)
            last_non_nan_idx = i                       
            i = i + 1
    return temp_x, temp_y


def nan_image(image):
    image_new = np.ones(np.shape(image))
    for i in range(np.shape(image)[0]):
        for j in range(np.shape(image)[1]):
            if np.isnan(image[i,j]):
                image_new[i,j]=np.nan
    return image_new

def windowaveraging_backward(y, y_index,  win_size):
    y_points = [y[y_index - i] for i in range(win_size)] 
    y_mean = np.mean(y_points)
    return y_mean
    
    
def normalizing(y):
#    x = removeNans(x)[0]
    mean = np.mean(y)
    std = np.std(y)
    y = (y-mean)/std   
    return y

def detrending(x,y):
    p = np.polyfit(x,y,6)
    f = np.polyval(p,x)
    y_new = y-f
    return y_new


    
def removespikes(y):
    k =2
    mean = np.mean(y)
    std = np.std(y)
    for ind in range(1, len(y)):
        if y[ind]>=mean+k*std or y[ind]<=mean-k*std:
            y[ind]=mean
            if ind-3<=0 or ind+3>=len(y):
                continue
            else:
                y[ind-1]= mean;y[ind+1] = mean
                y[ind-2]= mean;y[ind+2] = mean
                y[ind-3]= mean;y[ind+3] = mean
    return y
