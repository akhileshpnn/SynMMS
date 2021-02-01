from model import *
from boundary_conditions import *
from initial_conditions import *
from external_signal import *
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

class ReactionDiffusion2D:

    size = 1002 # size of the 2D grid (if python plot crashes while updating animation, use lower size. eg. size=302)
    dx = 1 # space step
    T = 100 # total time
    dt = .001 # time step
    
    input_folder = os.path.abspath(os.getcwd())+'\\'
    save_to='saved data'
    
    def __init__(self,model,nbc,light):
        self.model=model
        self.bc=nbc
        self.light=light
        
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
        
        U2, V2 = self.ic.around_steady_state()
#        Sp, Sf = self.ic.random()
        self.bc.update_boundaries([U2,V2])
        
        Light = self.light.add_light(self.model.rho, self.size)
        Pc=Light[1:-1, 1:-1]
        
        if animate == True:            
            plt.ion()
            fig = plt.figure(figsize=(4, 4))
#            fig.canvas.mpl_connect('close_event', self.handle_close)
            ax = fig.add_subplot(111)
#            ax.set_xlabel('coordinate x')
#            ax.set_ylabel('cordinate y')
            #im=ax.imshow(Sp,cmap=plt.cm.hot,origin='lower',interpolation='none',vmin=0,vmax=6)
            im=ax.imshow(U2,cmap=plt.cm.hot,origin='lower',interpolation='none')
            cbar = fig.colorbar(im, ticks=[0, np.max(U2)])
            cbar.ax.set_yticklabels(['0', '100%'])
        
        for n in range(self.N):    
            
            # diffusion in 2D
            deltaU2 = self.laplacian(U2)
            deltaV2 = self.laplacian(V2)
            
            U2_centre = U2[1:-1, 1:-1]
            V2_centre = V2[1:-1, 1:-1]
            
            # reaction
            dU2,dV2 = self.model.reaction([U2_centre, V2_centre],n*self.dt, Pc)
            
            # updating state variables
            U2[1:-1, 1:-1], V2[1:-1, 1:-1] = U2_centre + self.dt * (self.model.Du2 * deltaU2 + dU2),V2_centre + self.dt *(self.model.Dv2 * deltaV2 + dV2)
         
            # update boundaries
            self.bc.update_boundaries([U2,V2])
            
            if n%1000 ==0:    # plot/save at every 1000 frame         
                if save_data == True:
                    
                    isdir = os.path.isdir(self.input_folder+self.save_to+'\\') 
                    if isdir==True:
                        pass
                    elif isdir==False:
                        print('Create a folder with name : '+self.save_to+', in '+self.input_folder)
                        sys.exit()
                    
                    print('saved time point : '+str(n*self.dt)+'(arb. units)/'+str(self.N*self.dt)+'(arb. units)',end='\r')
                    
                    np.save(self.input_folder+self.save_to+'\\'+'U2_loop_'+str(n)+'.npy',U2[1:-1,1:-1])
                    np.save(self.input_folder+self.save_to+'\\'+'V2_loop_'+str(n)+'.npy',V2[1:-1,1:-1])
                    np.save(self.input_folder+self.save_to+'\\'+'Light cue_loop_'+str(n)+'.npy',Light[1:-1,1:-1])
                if animate==True:
                    im.set_array(U2)
                    fig.canvas.set_window_title('Time '+str(n)+' Total '+str(np.round(np.sum(U2[1:-1,1:-1])+np.sum(V2[1:-1,1:-1]))))
#                    fig.canvas.flush_events()
                    plt.pause(0.1)

if __name__ == '__main__':

    
    model = AuroraSubSystem() # defining the model for intergration   
    nbc = NeumannBoundaryConditions() # defining boundary conditions 
    light = Uniform() # defining light activation
    
    rd = ReactionDiffusion2D(model,nbc,light)
    rd.simulate(animate=True, save_data=None) # save==True means the state variables will be saved, animate==True means animated simulataneously

        
