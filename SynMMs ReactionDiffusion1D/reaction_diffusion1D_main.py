
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from itertools import product
from boundary_conditions import *
from model import *
from pheromone import *
import time
import os
import sys


class ReactionDiffusion1D:
    
    h = .1
    n = 2+1000
    t_total = 1.0e8   
    dt = 0.0001
    closed = False
    inc_load = None

    f=[] # no stimulus   
#    f=[0.7,1] # 1 stimulus between f[0]*t_total and f[1]*t_total  
#    f=[0.25,0.5,0.5,1]  # 2 stimuli in [f[0]*t_total , f[1]*t_total] and [f[2]*t_total , f[3]*t_total]
          
    input_folder = os.path.abspath(os.getcwd())+'\\'
    save_to='saved data'

    def __init__(self, model, boundary_conditions, pheromone):

        self.model = model
        self.boundary_conditions = boundary_conditions
        self.pheromone = pheromone

    def initialize_system(self):
        
#        seed_int=np.random.randint(0,500);
#        seed_int = 200 #Fig8c bottom,supp 9g
#        seed_int = 170 #Fig8b bottom, Fig1g, Fig8c top
        seed_int = 170 #Fig8b bottom, Fig1g
        print('random number seed :'+ str(seed_int))
        np.random.seed(seed_int)
        
        self.Z = np.zeros((self.n,), [('U1', np.double), ('V1', np.double),('U2', np.double),('V2', np.double)])
        
        if not self.inc_load==None:
            self.Z['U1'] = np.load(os.path.join(self.input_folder+'\\initial_condition\\','U1_initial_condition.npy'))
            self.Z['V1'] = self.model.c1/(self.n-2) - self.Z['U1']
            self.Z['U2'] = np.load(os.path.join(self.input_folder+'\\initial_condition\\','U2_initial_condition.npy'))
            self.Z['V2'] = self.model.c2/(self.n-2) - self.Z['U2']
            
            print('Inital conditions loaded...')
            print('Total MTs :',np.sum(self.Z['U1'][1:-1])+np.sum(self.Z['V1'][1:-1]))
            print('Total Aurora :',np.sum(self.Z['U2'][1:-1])+np.sum(self.Z['V2'][1:-1]))
        else:
            self.Z['U1'] = self.model.c1/(self.n-2)*np.random.random((self.n,))
            self.Z['V1'] = self.model.c1/(self.n-2) - self.Z['U1']
            self.Z['U2'] = self.model.c2/(self.n-2)*np.random.random((self.n,))
            self.Z['V2'] = self.model.c2/(self.n-2) - self.Z['U2']

        self.boundary_conditions.update_boundaries(self.Z)
        self.X = np.linspace(self.h*.5,(self.n-2.5)*self.h,self.n-2)
    
    def save_initial_conditions(self, save=None):
        
        if not save==None:
            input_folder_ini = self.input_folder+'initial_condition\\'
            np.save(input_folder_ini+'U1 initial condition.npy',self.Z['U1']);
            np.save(input_folder_ini+'V1 initial condition.npy',self.Z['V1']);
            np.save(input_folder_ini+'U2 initial condition.npy',self.Z['U2']);
            np.save(input_folder_ini+'V2 initial condition.npy',self.Z['V2']);
        

    def simulate(self, save=None, animate=None):
        
        ylims=[]
        
        self.initialize_system()
        self.save_initial_conditions(save=None)
        
        vars = ['U1', 'V1','U2', 'V2']
        U1, V1, U2, V2= self.Z['U1'], self.Z['V1'], self.Z['U2'], self.Z['V2']
        self.X = np.linspace(self.h*.5,(self.n-2.5)*self.h,self.n-2)

        
        if animate==True:
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
            line[0, 4], = ax[0, 2].plot(self.X, self.pheromone.add_pheromone_aur(0))
            line[0, 5], = ax[0, 0].plot(self.X, self.pheromone.add_pheromone_mt(0))

        for t in range(int(self.t_total)):
            
            
            
            if self.closed:
                return

            Lu1 = (U1[0:-2] - 2*U1[1:-1] + U1[2:])
            Lv1 = (V1[0:-2] - 2*V1[1:-1] + V1[2:])
            Lu2 = (U2[0:-2] - 2*U2[1:-1] + U2[2:])
            Lv2 = (V2[0:-2] - 2*V2[1:-1] + V2[2:])
            
            pher_profile_aur = self.pheromone.add_pheromone_aur(self.X,t)
            pher_profile_mt = self.pheromone.add_pheromone_mt(self.X,t)
            
            fu1v1, gu1v1, fu2v2, gu2v2= self.model.reaction(U1[1:-1], V1[1:-1], U2[1:-1], V2[1:-1],pher_profile_mt,pher_profile_aur)
            
            U1[1:-1] += self.dt * (self.model.du1 * Lu1/self.h**2 + fu1v1)
            V1[1:-1] += self.dt * (self.model.dv1 * Lv1/self.h**2 + gu1v1)
            U2[1:-1] += self.dt * (self.model.du2 * Lu2/self.h**2 + fu2v2)
            V2[1:-1] += self.dt * (self.model.dv2 * Lv2/self.h**2 + gu2v2)

            self.boundary_conditions.update_boundaries(self.Z)
            
            if animate==True and t % 1000 == 0: # animate at every 1000 loop
                
                for (row, col) in product(range(nrows), range(ncols)):
                    ind = row*ncols + col
                    line[row, col].set_ydata(self.Z[vars[ind]][1:-1])
                    if len(ylims) != 2:
                        ax[row, col].relim()
                        ax[row, col].autoscale_view() 
                line[0, 4].set_ydata(pher_profile_aur)
                line[0, 5].set_ydata(pher_profile_mt)

                total_MTs = np.round(np.sum(self.Z['U1'][1:-1])+np.sum(self.Z['V1'][1:-1]),2)
                total_Aurora = np.round(np.sum(self.Z['U2'][1:-1])+np.sum(self.Z['V2'][1:-1]))

                fig.canvas.set_window_title('Total MTs :' + str(total_MTs) + 'Total Aurora, loop :' + str(total_Aurora) + str(' ') + str(t))
                fig.canvas.flush_events()
                time.sleep(.00001)
                
            if save==True and t% 200000 == 0: # save at every 200000 loop
                
                isdir = os.path.isdir(self.input_folder+self.save_to+'\\') 
                if isdir==True:
                    pass
                elif isdir==False:
                    print('Create a folder with name : '+self.save_to+', in '+self.folder_save)
                    sys.exit()
                
                print('saved loop : '+str(t),end='\r')
                
                np.save(self.input_folder+self.save_to+'\\'+'U1_loop_'+str(t)+'.npy',self.Z['U1']);
                np.save(self.input_folder+self.save_to+'\\'+'V1_loop_'+str(t)+'.npy',self.Z['V1'])
                np.save(self.input_folder+self.save_to+'\\'+'U2_loop_'+str(t)+'.npy',self.Z['U2']);
                np.save(self.input_folder+self.save_to+'\\'+'V2_loop_'+str(t)+'.npy',self.Z['V2'])
 
                np.save(self.input_folder+self.save_to+'\\'+'Pheromone cue_loop_'+str(t)+'.npy',pher_profile_aur);
         
                
        plt.ioff()
        plt.show()
                     
    def handle_close(self, evt):
        self.closed = True


if __name__ == '__main__':

    
    pbc = PeriodicBoundaryConditions()
    model = SynMMS()
    pheromone = TimeVaryingPheromone()

    rd = ReactionDiffusion1D(model, pbc, pheromone)

    rd.simulate(save=True, animate=None)
