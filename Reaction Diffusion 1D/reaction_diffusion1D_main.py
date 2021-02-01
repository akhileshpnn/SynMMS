
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from itertools import product
import time
import os
import sys

from boundary_conditions import *
from model import *
from external_signal import *

class ReactionDiffusion1D:
    
    h = .1 # grid spacing in 1D
    n = 2+1000 # total number of grid points in 1D 
    N = 1.0e8 # number of iterations ,Total time for integration is N*dt (arb. units)
    dt = 0.0001 # time step for integration

    closed = False
    
    # light stimulus timings
    f=[] # no stimulus   
#    f=[0.7,1] # example of one stimulus between f[0]*N and f[1]*N  
#    f=[0.25,0.5,0.5,1]  #  example of two stimuli in [f[0]*N , f[1]*N] and [f[2]*N , f[3]*N]
          
    input_folder = os.path.abspath(os.getcwd())+'\\'
    save_to='saved data'

    def __init__(self, model, boundary_conditions, light):

        self.model = model
        self.boundary_conditions = boundary_conditions
        self.light = light

    def initialize_system(self):
        
        seed_int=np.random.randint(0,500);

#        seed_int = 200 # Fig. 8c bottom, Supp. 9g
#        seed_int = 170 # Fig. 8b bottom, Fig. 1g, Fig. 8c top
#        seed_int = 170 # Fig. 8b bottom, Fig. 1g
        print('random number seed :'+ str(seed_int))
        np.random.seed(seed_int)
        
        self.Z = np.zeros((self.n,), [('U1', np.double), ('V1', np.double),('U2', np.double),('V2', np.double)])
        
        # initializing the variables
        self.Z['U1'] = self.model.c1/(self.n-2)*np.random.random((self.n,))
        self.Z['V1'] = self.model.c1/(self.n-2) - self.Z['U1']
        self.Z['U2'] = self.model.c2/(self.n-2)*np.random.random((self.n,))
        self.Z['V2'] = self.model.c2/(self.n-2) - self.Z['U2']
        
        #setting boundary conditions
        self.boundary_conditions.update_boundaries(self.Z)
        #defining the space in 1D
        self.X = np.linspace(self.h*.5,(self.n-2.5)*self.h,self.n-2)
    
    def simulate(self, save=None, animate=None):
        
        ylims=[]
        
        self.initialize_system()
        
        vars = ['U1', 'V1','U2', 'V2']
        U1, V1, U2, V2= self.Z['U1'], self.Z['V1'], self.Z['U2'], self.Z['V2']
        self.X = np.linspace(self.h*.5,(self.n-2.5)*self.h,self.n-2)

        
        if animate==True: # for animating the system evolution in realtime
            plt.ion()
            dpi = 72.0
            nrows, ncols = 1, 4
            fig, ax = plt.subplots(nrows, ncols, squeeze=False, dpi=dpi, facecolor="white")
            fig.canvas.mpl_connect('close_event', self.handle_close)
            plt.autoscale(tight=True)
            line = np.empty((1,6), dtype=matplotlib.lines.Line2D)
            
            for (row,col) in product(range(nrows), range(ncols)):
                ind = row*ncols+col
                plt.sca(ax[row,col])
                plt.xlim(0, self.n * self.h)
                if len(ylims)==2:
                    plt.ylim(0, ylims[ind])
                else:
                    ax[row, col].autoscale(True)
                
                line[row, col], = ax[row, col].plot(self.X, self.Z[vars[ind]][1:-1])
                ax[row,col].set_title(vars[ind])
                ax[row,col].set_xlabel('space')
            line[0, 4], = ax[0, 2].plot(self.X, self.light.add_light_aurora(self.X,0))

        for n in range(int(self.N)): # integrating loops
            
            if self.closed:
                return
            
            # laplacian 1D
            Lu1 = (U1[0:-2] - 2*U1[1:-1] + U1[2:])
            Lv1 = (V1[0:-2] - 2*V1[1:-1] + V1[2:])
            Lu2 = (U2[0:-2] - 2*U2[1:-1] + U2[2:])
            Lv2 = (V2[0:-2] - 2*V2[1:-1] + V2[2:])
            
            #light cue (light stimulation)
            rho_profile_aurora = self.light.add_light_aurora(self.X,n)
            
            # reaction term
            fu1v1, gu1v1, fu2v2, gu2v2= self.model.reaction(U1[1:-1], V1[1:-1], U2[1:-1], V2[1:-1],rho_profile_aurora)
            
            # updating state variables with reaction and diffusion 
            U1[1:-1] += self.dt * (self.model.Du1 * Lu1/self.h**2 + fu1v1)
            V1[1:-1] += self.dt * (self.model.Dv1 * Lv1/self.h**2 + gu1v1)
            U2[1:-1] += self.dt * (self.model.Du2 * Lu2/self.h**2 + fu2v2)
            V2[1:-1] += self.dt * (self.model.Dv2 * Lv2/self.h**2 + gu2v2)
            
            #updating boundary
            self.boundary_conditions.update_boundaries(self.Z)
            
            if animate==True and n % 1000 == 0: # animate at every 1000 loop
                
                for (row, col) in product(range(nrows), range(ncols)):
                    ind = row*ncols + col
                    line[row, col].set_ydata(self.Z[vars[ind]][1:-1])
                    if len(ylims) != 2:
                        ax[row, col].relim()
                        ax[row, col].autoscale_view() 
                line[0, 4].set_ydata(rho_profile_aurora)

                total_MTs = np.round(np.sum(self.Z['U1'][1:-1])+np.sum(self.Z['V1'][1:-1]),2)
                total_Aurora = np.round(np.sum(self.Z['U2'][1:-1])+np.sum(self.Z['V2'][1:-1]))

                fig.canvas.set_window_title('Total MTs :' + str(total_MTs) + 'Total Aurora, loop :' + str(total_Aurora) + str(' ') + str(n))
                fig.canvas.flush_events()
                time.sleep(.00001)
                
            if save==True and n% 200000 == 0: # save at every 200000 loop
                
                isdir = os.path.isdir(self.input_folder+self.save_to+'\\') 
                if isdir==True:
                    pass
                elif isdir==False:
                    print('Create a folder with name : '+self.save_to+', in '+self.input_folder)
                    sys.exit()
                
                print('saved time point : '+str(n*self.dt)+'(arb. units)/'+str(self.N*self.dt)+'(arb. units)',end='\r')
                
                np.save(self.input_folder+self.save_to+'\\'+'U1_loop_'+str(n)+'.npy',self.Z['U1'])
                np.save(self.input_folder+self.save_to+'\\'+'V1_loop_'+str(n)+'.npy',self.Z['V1'])
                np.save(self.input_folder+self.save_to+'\\'+'U2_loop_'+str(n)+'.npy',self.Z['U2'])
                np.save(self.input_folder+self.save_to+'\\'+'V2_loop_'+str(n)+'.npy',self.Z['V2'])
 
                np.save(self.input_folder+self.save_to+'\\'+'Light cue_loop_'+str(n)+'.npy',rho_profile_aurora)
         
                
        plt.ioff()
        plt.show()
                     
    def handle_close(self, evt):
        self.closed = True


if __name__ == '__main__':
   
    pbc = PeriodicBoundaryConditions() # periodic boundary condition
    model = SynMMS() # defining the model for intergration
    light = TimeVaryingCue() # light activation of Aurora translocation

    rd = ReactionDiffusion1D(model, pbc, light)

    rd.simulate(save=None, animate=True) # save==True means the state variables will be saved, animate==True means animated simulataneously
