import numpy as np
import os
import spatial_recurrence_plot_quantifications as rpq
import pandas as pd
import matplotlib.pyplot as plt
from matrix_shuffling import *
import sys

class SpatialRecurrencePlot:
    input_folder = ''
    output_folder = ''
    filename = ''
    filename_basic = ''
    image = None
    
    def __init__(self, input_folder, output_folder, filename, eps):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.filename = filename
        self.filename_basic = filename.split('.')[0]
        self.eps = eps
        self.read_image()
        
        
    def read_image(self):
        
        self.filename_with_extension = self.filename+'.npy'
        self.image = np.load( os.path.join(self.input_folder, self.filename_with_extension))
        Nx,Ny=np.shape(self.image)
        
        span = 50 # Crop the original image to (2*span x 2*span)
        if Nx<2*span:
            print('Cannot crop image with dimension ('+str(Nx)+','+str(Ny)+') to ('+str(2*span)+','+str(2*span)+')') 
            sys.exit()
        else:
            self.image = self.image[int(Nx*0.5)-span:int(Nx*0.5)+span,int(Nx*0.5)-span:int(Nx*0.5)+span]
        
        
        fig, ax = plt.subplots()
        ax.set_xlabel('Coordinate X')
        ax.set_ylabel('Coordinate Y')                
        im=ax.imshow(self.image,cmap=plt.cm.hot,origin='lower',interpolation='none',vmin=np.min(self.image),vmax=np.max(self.image))
        fig.colorbar(im)
#        ax.set_title('Cropped image')
        isdir = os.path.isdir(self.output_folder) 
        if isdir==True:
            plt.savefig(os.path.join(self.output_folder,self.filename+'_eps='+str(self.eps)+'.png'), dpi=300)
            pass
        elif isdir==False:
            print('Create a folder with name : '+save_to+', in '+parent_folder)
            sys.exit()

        
        if shuffle==True:
            binary_image,binary_shuffled,self.image, thresh, mean_image = shuffling(self.image, self.output_folder, self.filename, out='Full')
            self.filename = self.filename+'_'+str(trial_num)
            self.save_npy()
            plt.figure()
            plt.imshow(self.image)
            plt.savefig(os.path.join(self.output_folder,self.filename+'_eps='+str(self.eps)+'.png'), dpi=300)
            plt.figure()
            plt.imshow(binary_image)
            plt.title('threshold = '+str(thresh)+', image mean = '+str(mean_image))
            plt.savefig(os.path.join(self.output_folder,self.filename+'_image_thresholded_eps='+str(self.eps)+'.png'), dpi=300)
            
            plt.figure()
            plt.imshow(binary_shuffled)           
            plt.savefig(os.path.join(self.output_folder,self.filename+'_shuffled_thresholded_eps='+str(self.eps)+'.png'), dpi=300)
    
    def read_csv_file(self):
        self.data = pd.read_csv(os.path.join(self.input_folder, self.filename+'.csv'), delimiter=',')        
        self.set_variables()
                
    def save_npy(self):
        self.filename_with_extension = self.filename+'.npy'
        np.save(os.path.join(self.input_folder, self.filename_with_extension), self.image)

    
    def plot_recurrence_plot(self):
            import matplotlib.pylab as plt
            from mpl_toolkits.mplot3d import Axes3D
            recurrence_matrix = np.load(os.path.join(self.output_folder, 'Recurence_matrix_of_'+self.filename+'_eps='+str(self.eps)+'.npy'))
            
            nx1, ny1, nx2, ny2  = np.shape(recurrence_matrix)
            recurrence_matrix_im_proj = np.zeros((nx1, ny1, nx2))*np.nan
            X, Y, Z = [],[],[]
            for i in range(nx1):
                print(i, end='\r')
                for j in range(ny1):
                    for k in range(nx2):
                        if recurrence_matrix[i,j,k,k] == 1:
                            x_loc = [i]; y_loc = [j]; z_loc = [k] 
                            X.append(x_loc)
                            Y.append(y_loc) 
                            Z.append(z_loc)
            
            fig,ax = plt.subplots()
            ax = plt.axes(projection='3d')
            
            ax.set_xlim(0,nx1-1)
            ax.set_ylim(0,ny1-1)
            ax.set_zlim(0,nx2-1)
            
            ax.set_xlabel('matrix index '+'$i$')
            ax.set_ylabel('matrix index '+'$j$')
            ax.set_zlabel('matrix index '+'$k$')
            
            for j in range(len(X)):
                if X[j][0]==0:
                    ax.scatter3D([X[j]]*len(Z[j]), [Y[j]]*len(Z[j]), Z[j], c = 'k', marker = '.',s = 1)
                if Y[j][0]==ny1-1:
                    ax.scatter3D([X[j]]*len(Z[j]), [Y[j]]*len(Z[j]), Z[j], c = 'k', marker = '.',s = 1)
                if Z[j][0]==0:
                    ax.scatter3D([X[j]]*len(Z[j]), [Y[j]]*len(Z[j]), Z[j], c = 'k', marker = '.',s = 1)
            
            if not output_folder is None:
                plt.savefig(os.path.join(output_folder,self.filename.split('.')[0]+'.png'), dpi=300)
            plt.show()
        
    def calculate_metrics(self, lmin, detrend):

        metrics = [rpq.entropy]
        quantity_metrics = np.zeros((len(metrics),))

        for ind, metric in enumerate(metrics):
            res = metric(self.image, self.output_folder, self.filename, self.eps,lmin ,detrend)
            quantity_metrics[ind] = res

        return quantity_metrics
        
        
    def save_metrics(self, output_folder, quantity_metrics):
        np.savetxt(os.path.join(output_folder, self.filename+'_eps('+str(self.eps)+').csv'), quantity_metrics, delimiter=',')                
    
if __name__=='__main__':
    
    
    load_image_from='saved data'
    save_to='saved recurrence data'
    parent_folder=os.path.abspath(os.getcwd())+'\\'
    input_folder = parent_folder+load_image_from+'\\'
    output_folder = parent_folder+save_to+'\\'


    filename = 'U2'
    eps = 0.2; # epsilon value used for Supplementary Fig. 4h
    # eps = 0.3; # epsilon value used for Fig. 3f
    lmin = 3; # minimum size of diagonal hyper-surface
    shuffle=None # 'shuffle =True' can be used to test whether the Entropy value drops upon pattern shuffling 
    Entropy = []
    
    Entropy_estimation_at=[0, 999000]
    
    if shuffle != True:
        for t in Entropy_estimation_at:
            filename = filename+'_loop_'+str(t)
            rp = SpatialRecurrencePlot(input_folder, output_folder, filename, eps)
            
            print('Estimating Entropy from '+filename)
            quantity_metrics = rp.calculate_metrics(lmin, detrend=None)
            Entropy.append(quantity_metrics[0])
            print('Entropy ='+str(quantity_metrics[0]))
            
            rp.save_metrics(output_folder, quantity_metrics)
            
            filename = 'U2'
            
    if shuffle == True:
        for t in Entropy_estimation_at:
            filename = filename+'_loop_'+str(t)
            Iter = [1]
            
            for trial_num in Iter:
            
                Entropy = []
                print("Trial_num ",trial_num)
                print("File loading....",filename)     
                rp = SpatialRecurrencePlot(input_folder, output_folder, filename, eps)
                print('Estimating Entropy from '+filename)
                quantity_metrics = rp.calculate_metrics(lmin, detrend=None)
                Entropy.append(quantity_metrics[0])
                print('Entropy ='+str(quantity_metrics[0]))
                
                rp.save_metrics(output_folder, quantity_metrics) 
                
                filename = 'U2_loop_'+str(t)
