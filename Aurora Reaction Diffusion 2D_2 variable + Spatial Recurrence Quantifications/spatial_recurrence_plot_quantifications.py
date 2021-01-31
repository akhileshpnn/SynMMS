import numpy as np
import os
from PIL import Image
from utils_ import *
import matplotlib.pyplot as plt


def detrending_3D(image, output_folder, file_name, eps):
#    image = np.load(os.path.join(input_folder,file_name+'.npy'))
    (nx, ny) = np.shape(image)
    image_new = np.zeros((nx,ny))
    x_ = np.arange(nx); y_ = np.arange(ny)
    for i in range(nx):
        y_,nan_removed = removeNans(y_,image[i,:])
        detrended = detrending(y_,nan_removed)
        image_new[i,:] = detrended
    for j in range(ny):
        x_,nan_removed = removeNans(x_,image[:,j])
        detrended = detrending(x_,nan_removed)
        image_new[:,j] = detrended
    np.save(os.path.join(input_folder, file_name+'_detrended.npy'), image_new)
    return image_new

def filling_from_poisson(image_npy):
    
    for i in range(l):
        for j in range(m):
            if A[i,j] < 5:
                A_new[i,j] = np.random.poisson(np.nanmean(A))
    return image_npy

def plot_3D(input_folder, output_folder, filename, eps):
    from mpl_toolkits.mplot3d import Axes3D
    import scipy.misc
#        from cv2_rolling_ball import subtract_background_rolling_ball
    image = np.load(os.path.join(input_folder,filename+'.npy'))
    im = image#np.load(os.path.join(self.input_folder,self.filename+'.npy'))
#        im_unint8 = im.astype('uint8')
#        img_, background = subtract_background_rolling_ball(im_unint8, 30, light_background=True,
#                                     use_paraboloid=False, do_presmooth=True)
#        im = scipy.misc.imresize(im, 0.15, interp='cubic')
    xx, yy = np.mgrid[0:im.shape[0], 0:im.shape[1]]
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(xx, yy, im ,rstride=1, cstride=1, cmap=plt.cm.gray,
            linewidth=0)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Intensity')
    
    plt.savefig(os.path.join(output_folder, filename+'_3D_eps('+str(eps)+').png'), dpi=300)
    # show it
#    plt.show()    

def normalize_image(image):
    mean = np.nanmean(image)
    std = np.nanstd(image)
    image = (image-mean)/std
    return image

def setting_to_nan(M):
    fraction = 0.1
    nx,ny = np.shape(M)
    num_of_Nans = int(nx*ny*fraction)
    folder_temp = '\\\\billy\\abt2\\group\\agkoseska\\Nandan_Akhilesh\\GUV_project\\20190226\\images\\'
    nan_matrix_img = Image.open(os.path.join(folder_temp, 'pplambda_Nans.tif'))
    nan_matrix = np.array(nan_matrix_img)
    nan_matrix = np.pad(nan_matrix, [(0,20),(0,0)],'constant',constant_values=(1))
    nan_matrix = nan_matrix[0:min(np.shape(M)[0],np.shape(nan_matrix)[0]):,0:min(np.shape(M)[1],np.shape(nan_matrix)[1])]
    ones_matrix = np.ones(np.shape(M))
    for i in range(np.shape(nan_matrix)[0]):
        for j in range(np.shape(nan_matrix)[1]):
            if nan_matrix[i,j]==255:
                nan_matrix[i,j] = 1
#    nan_matrix = nan_matrix[nan_matrix == 255]=1
#    nan_matrix = np.ones((nx,ny))
#    points = [[3,5],[3,6],[3,7],[4,6],[4,7],[5,8],[30,85],[30,86],[31,85],[32,85],[33,85],[21,40],[10,40],[11,40],[12,40],[11,41],[30,22],[30,23],[30,24],[]
#    for point in points:
#        nan_matrix[point[0],point[1]] = np.nan
#    nan_matrix = nan_matrix[40:80:,0:80]
    M = M*nan_matrix
    return M
    
    

def spatial_recurrence_matrix(image, output_folder, file_name, eps, detrend = None):
    print('Calculating recurrence matrix')
    
    if not detrend==None:
        image = detrending_3D(image, output_folder, file_name, eps)
        image = img*nan_image_
        plt.figure()
        plt.imshow(image)
        plt.savefig(os.path.join(output_folder, file_name+'_detrended_eps='+str(eps)+'_.png'), dpi=300)


    #img = normalize_image(img)
    nx, ny = np.shape(image)
    R = np.zeros((nx, ny, nx, ny))
    for i in range(nx):
        print(str(i)+'/'+str(nx), end = '\r')
        for j in range(ny):
            x_1 = image[i, j]
            for k in range(nx):  
                for m in range(ny):
                    if np.isnan(image[k, m]-x_1):
#                        R[i,j,k,m] = np.nan
                        continue
                    else:
                        if abs(image[k, m]-x_1)<eps:
                            R[i,j,k,m] = 1
                        elif abs(image[k, m]-x_1)>eps:
                            R[i, j, k, m] = 0
    np.save(os.path.join(output_folder, 'Recurence_matrix_of_'+file_name+'_eps='+str(eps)+'.npy'), R)
    return R

def get_diagonal_lines(image, output_folder, file_name, eps, detrend):
    
    if not os.path.isfile(os.path.join(output_folder, 'Recurence_matrix_of_'+file_name+'_eps='+str(eps)+'.npy')):
        rp_matrix = spatial_recurrence_matrix(image, output_folder, file_name, eps, detrend)
    else:
        rp_matrix = np.load(os.path.join(output_folder, 'Recurence_matrix_of_'+file_name+'_eps='+str(eps)+'.npy'))
        print('Recurrence matrix loaded')
       
    print('Calculating distribution of diagonal hyper-surfaces')
    L = []
    nx1, ny1, nx2, ny2 = np.shape(rp_matrix)
    history = np.zeros(np.shape(rp_matrix))
    for i1 in range(1, nx1-1):
        print(i1, end = '\r')
        for i2 in range(1, ny1-1):
            for j1 in range(1, nx2-1):
                for j2 in range(1, ny2-1):
                    if [i1, i2] == [j1, j2]:
                        break
                    else:
                        l_forward = 1; l_backward = 1
                        if np.isnan(rp_matrix[i1, i2, j1, j2]):
                            continue
                        else:                            
                            if rp_matrix[i1, i2, j1, j2] == 0:
                                history[i1, i2, j1, j2] = 2
                                break
                            elif rp_matrix[i1, i2, j1, j2] == 1 and history[i1, i2, j1, j2] == 2:
                                history[i1, i2, j1, j2] = 2
                                break
                            elif rp_matrix[i1, i2, j1, j2] == 1 and history[i1, i2, j1, j2] == 0:
                                while 0 <= (i1 + l_forward) < nx1 and 0<= (i2 + l_forward) < ny1 and 0 <= (j1 + l_forward) < nx2 and 0 <= (j2 + l_forward) < ny2 :
                                    if rp_matrix[i1+l_forward, i2+l_forward, j1+l_forward, j2+l_forward]==0 or np.isnan(rp_matrix[i1+l_forward, i2+l_forward, j1+l_forward, j2+l_forward]):                                        
                                        break
                                    elif rp_matrix[i1+l_forward, i2+l_forward, j1+l_forward, j2+l_forward]==1:
                                        history[i1+l_forward, i2+l_forward, j1+l_forward, j2+l_forward] = 2
                                        l_forward = l_forward + 1                  
                                while 0 <= (i1 - l_backward) < nx1 and 0<= (i2 - l_backward) < ny1 and 0 <= (j1 - l_backward) < nx2 and 0 <= (j2 - l_backward) < ny2 :
                                    if rp_matrix[i1-l_backward, i2-l_backward, j1-l_backward, j2-l_backward]==0 or np.isnan(rp_matrix[i1-l_backward, i2-l_backward, j1-l_backward, j2-l_backward]):
                                        break
                                    elif rp_matrix[i1-l_backward, i2-l_backward, j1-l_backward, j2-l_backward]==1:
                                        history[i1-l_backward, i2-l_backward, j1-l_backward, j2-l_backward] = 2
                                        l_backward = l_backward + 1
                                history[i1, i2, j1, j2] = 2

                    L.append(l_backward + l_forward - 2)
    np.save(os.path.join(output_folder, 'Distri_of_diag_'+file_name+'_eps='+str(eps)+'.npy'), L)   
    return L

def Remove_zeros(arr,index): 
   
    n = len(index[0])
    temp = [0] * n; 
  
    # arr[i] should be 
        # present at index[i] index 
    for i in range(0,n): 
        temp[i] = arr[index[0][i]] 
  
    return temp
#index = np.nonzero(A)
#A_new = Remove_zeros(A,index[0])

def histogram(a):
    Hist = []
    index_of_nonzeros = np.nonzero(a)
    a = Remove_zeros(a,index_of_nonzeros)
    a_temp = np.unique(a)
    for a_ in a_temp:
#        print(a_, end='\r')
        count = 0
        index_to_delete = []
        for j in range(len(a)):           
            if a[j]==a_:
                count= count+1
                index_to_delete.append(j)
        a = np.delete(a,index_to_delete)
        Hist.append(count)
    return Hist

def get_diagonal_distribution(L, lmin):
    Nl = 0
    for d in L:
        if d>=lmin:
            Nl = Nl+1
    Pl = histogram(L)#np.histogram(ds, range = bins)[0]# range=(min(ds_unique), max(ds_unique)))[0]
#    Pl = Pl[Pl!=0]
    return Pl, Nl


def determinism(image, output_folder, file_name, eps, lmin, detrend):
    if not os.path.isfile(os.path.join(output_folder,'Distri_of_diag_'+file_name+'_eps='+str(eps)+'.npy')):
        L = get_diagonal_lines(image, output_folder, file_name, eps, detrend)
    else:
        L = np.load(os.path.join(output_folder, 'Distri_of_diag_'+file_name+'_eps='+str(eps)+'.npy'))
        print('Diagonal lengths loaded')
    
    Pl, Nl = get_diagonal_distribution(L, lmin)
    index_of_nonzeros = np.nonzero(L)
    L = Remove_zeros(L,index_of_nonzeros)
    L_unique = np.sort(np.unique(L))
    product = L_unique*Pl
    det_top = 0
    for ind, l in enumerate(L_unique):
        if l>=lmin:
            det_top = det_top + product[ind]  
    det_bottom = np.sum(product)
    
    det = det_top/det_bottom
    return det

def entropy(image, output_folder, file_name, eps, lmin, detrend):
    if not os.path.isfile(os.path.join(output_folder,'Distri_of_diag_'+file_name+'_eps='+str(eps)+'.npy')):
        L = get_diagonal_lines(image,output_folder, file_name, eps, detrend)
    else:
        L = np.load(os.path.join(output_folder, 'Distri_of_diag_'+file_name+'_eps='+str(eps)+'.npy'))
        print('Diagonal lengths loaded')
    Pl, Nl = get_diagonal_distribution(L, lmin)
    pl = np.array(Pl)/Nl
    index_of_nonzeros = np.nonzero(L)
    L = Remove_zeros(L,index_of_nonzeros)
    L_unique = np.sort(np.unique(L))
    product = pl*np.log(pl)
    ent = 0
    for ind, l in enumerate(L_unique):
        if l>=lmin:
            ent = ent + product[ind]          
    return -ent
