from model import *
from boundary_conditions import *
from initial_conditions import *
from pheromone import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import sys

class ReactionDiffusion2D:

    size =102 # size of the 2D grid
    dx =2 # space step
    T = 3000# total time
    dt = .001 # time step
    load =False
    
    input_folder = os.path.abspath(os.getcwd())+'\\'
    save_to='saved data'
        
    def __init__(self,model,bc,pher):
        self.model=model
        self.bc=bc
        self.pher=pher
        
    def initialize(self):
        self.ic = InitialConditions(self.model,self.size)
        self.N = int(self.T / self.dt)  # number of iterations
    
    def laplacian(self,Z):
        Ztop = Z[0:-2, 1:-1]
        Zleft = Z[1:-1, 0:-2]
        Zbottom = Z[2:, 1:-1]
        Zright = Z[1:-1, 2:]
        Zcenter = Z[1:-1, 1:-1]
        return (Ztop + Zleft + Zbottom + Zright -4 * Zcenter) / self.dx**2
    
    def handle_close(self, evt):
        self.closed = True
    
    def simulate(self, animate=False, save_data=False):
        self.initialize()
        
        Sp, Sfm, Sfc = self.ic.around_steady_state()
            
        bc.update_boundaries([Sp,Sfm, Sfc])
        if len(np.argwhere(Sp<0))!=0 or len(np.argwhere(Sfm<0))!=0 or len(np.argwhere(Sfc<0))!=0:
            print('negative values')
        print('Total Aurora ',np.round(np.sum(Sp[1:-1,1:-1])+np.sum(Sfm[1:-1,1:-1])+np.sum(Sfc[1:-1,1:-1]),2))
        print(self.N)
        
        Pher = self.pher.add_pheremone(self.model.pheromone_ampl, self.size)
        Pc=Pher[1:-1, 1:-1]
        
        if animate == True:            
            plt.ion()
            fig = plt.figure(figsize=(4, 4))
#            plt.autoscale()
#            fig.canvas.mpl_connect('close_event', self.handle_close)
            ax = fig.add_subplot(111)
            ax.set_xlabel('Cordinate x')
            ax.set_ylabel('Cordinate y')
#            im=ax.imshow(Sp,cmap=plt.cm.hot,origin='lower',interpolation='none',vmin=0,vmax=20)
            im=ax.imshow(Sp,cmap=plt.cm.hot,origin='lower',interpolation='none')
            fig.colorbar(im)
        
        for i in range(self.N):    

            deltaSp = self.laplacian(Sp)
            deltaSfm = self.laplacian(Sfm)
            deltaSfc = self.laplacian(Sfc)
            
            Sp_centre = Sp[1:-1, 1:-1]
            Sfm_centre = Sfm[1:-1, 1:-1]
            Sfc_centre = Sfc[1:-1, 1:-1]
            
            dSp,dSfm,dSfc = self.model.reaction([Sp_centre, Sfm_centre, Sfc_centre],i*self.dt, Pc)
            Sp[1:-1, 1:-1], Sfm[1:-1, 1:-1], Sfc[1:-1, 1:-1] = Sp_centre + self.dt * (self.model.Dm1 * deltaSp + dSp),Sfm_centre + self.dt *(self.model.Dm2 * deltaSfm + dSfm),Sfc_centre + self.dt *(self.model.Dc * deltaSfc + dSfc)
         
            # update boundaries
            bc.update_boundaries([Sp,Sfm,Sfc])
            
            if i%1000 ==0:        
                print (str(i), end="\r")
                if save_data == True:
                    
                    isdir = os.path.isdir(self.input_folder+self.save_to+'\\') 
                    if isdir==True:
                        pass
                    elif isdir==False:
                        print('Create a folder with name : '+self.save_to+', in '+self.input_folder)
                        sys.exit()
                    
                    np.save(self.input_folder+self.save_to+'\\'+'U2_loop_'+str(i)+'.npy',Sp[1:-1,1:-1])
                    np.save(self.input_folder+self.save_to+'\\'+'V2_loop_'+str(i)+'.npy',Sf[1:-1,1:-1])
                    np.save(self.input_folder+self.save_to+'\\'+'Light cue_loop_'+str(i)+'.npy',Pher[1:-1,1:-1])
                if animate==True:
                    im.set_array(Sp)
                    fig.canvas.set_window_title('Time '+str(i)+' Total '+str(np.round(np.sum(Sp[1:-1,1:-1])+np.sum(Sf[1:-1,1:-1]))))
#                    fig.canvas.flush_events()
                    plt.pause(0.1)

if __name__ == '__main__':

    model = AuroraSubSystem3var()    
    nbc = NeumannBoundaryConditions()
    pheromone = Uniform()
    
    rd = ReactionDiffusion2D(model,nbc,pheromone)
    rd.simulate(animate=True, save_data=None)

    


        
