import numpy as np
import matplotlib.pyplot as plt
import os
from mt_simulation_attributes import *
from mt_snapshot import *
import sys

class MtSimulation:
    
    Nmtmax = 20; # number of MTs max
    Rcell = 25; # cell radius minimum (µm)
    Nsteps = 2001; # number of steps
    
    Vcell = 1000; # volume of the cell (µm3)
    TubulinTotalConc = 35; # total concentration of tubulin (µM)
    TubulinTotalNumb = TubulinTotalConc*(602)*Vcell; # total # of tubulins

    knuc = 0.0005; # nucleation rate (s-1)
    vg = 0.192; # growth velocity (µm/s)
    Vs = 0.218; # shrinkage velocity (µm/s)

    kc=0.015; # catastrophic frequency (s-1)
    dt = 1; # time step (s)
    time=np.arange(0,Nsteps,dt)
    pnuc = 1-np.exp(-knuc*dt); #  of nucleation
#    A = 0.2; #  (µm/s )
#    k = 0.16; #  per µM
    kres = 0.175; # rescue frequency (s-1)
    
    Lg = vg*dt; # length added in one step (µm)
    Ls = Vs*dt; # length decreased in one step (µm)
    
    ## initializing various arrays for updating
    Deforming_MT_id_ = [[]]
    Attracted_neighbours_ = [[]]
    Theta_=[]
    Force_of_attraction_ = []  
    Lmt = np.zeros(Nmtmax); # initialize MT lengths
    Tub = np.zeros(Nsteps)
    TotalMTlength = sum(Lmt); #µm
    LLmt = np.zeros((Nmtmax,Nsteps)) # entries in the columns are the current length of each microtubule. 
    Kc = np.ones((Nmtmax,Nsteps)); # catastophic frequency update
    Force_from_membrane_ = np.zeros((Nmtmax,Nsteps)); 
    
    state = np.ones(Nmtmax); # initialize MT states
    Pc = np.ones((Nmtmax,Nsteps)); # catastophic probability update
    Vg = np.ones((Nmtmax,Nsteps)); # Growth velocity
    
    X=np.linspace(10,2*np.pi+10,Nmtmax) # Bounded region where MTs are positioned
    X_position = [X]
    TubulinFreeNumb = TubulinTotalNumb-1624*TotalMTlength; # Updating free tubulin
    
    lf=0.05 # scaling factor in σ(Δh) 
    
    
    save=True     
    folder_save = os.path.abspath(os.getcwd())+'\\'
    save_to='saved data4'

    def __init__(self):
        return
    def initialize_system(self):
        
        for i in range(self.Nmtmax):
            self.Lmt[i]=np.random.randint(0,30)*self.Lg 
            
        self.Tub[0] = self.TubulinFreeNumb/(602*self.Vcell); # free tubulin conc (µM) at t=0
        self.Kc[:,0] = self.kc*self.Kc[:,0]        
        self.pc = 1-np.exp(-self.kc*self.dt)
        self.Pc[:,0] = self.pc*np.ones(len(self.Pc[:,0]))
        self.Vg[:,0] = self.vg*np.ones(len(self.Vg[:,0]))
        self.LLmt[:,0] = self.Lmt; #Initiating all the MTs with random lengths less than the cell periphery       
    
    def plot_free_tubulin(self):
        plt.figure()
        plt.plot(self.time[1:],self.Tub[1:],'k-',lw=2.0)
        plt.xlabel('Time')
        plt.ylabel('Free Tubulin Conc.')
        plt.show()
    
    def av_protrusion(self):
        positive_protrusion=np.zeros(self.Nsteps)*np.nan
        for i in range(self.Nsteps):
            indxs=np.argwhere((self.LLmt[:,i]-Rcell)>0)
            if len(indxs)!=0:
                positive_protrusion[i]=np.mean((self.LLmt[:,i]-Rcell)[indxs])
        plt.figure()
        plt.plot(self.time,positive_protrusion,'k-',lw=2.0)
        plt.xlabel('Time')
        plt.ylabel('Av. protrusion length')
        plt.show()
    
    def av_velocity(self):
        v=np.zeros(self.Nsteps)*np.nan
        for i in range(self.Nsteps):
            v[i]=np.mean(self.Vg[:,i])
        plt.figure()
        plt.plot(self.time,v,'k-',lw=2.0)
        plt.xlabel('Time')
        plt.ylabel('Av. velocity')
        # plt.xlim(0,250)
        plt.show()
    
    def cat_freq(self):
        kc=np.zeros(self.Nsteps)*np.nan
        for i in range(self.Nsteps):
            kc[i]=np.mean(self.Kc[:,i])
        plt.figure()
        plt.plot(self.time[1:],kc[1:],'k-',lw=2.0)
        plt.xlabel('Time')
        plt.ylabel('Av. Cat freq')
        # plt.xlim(0,250)
        plt.show()
    
    def simulate(self):
        
        self.initialize_system()
               
        
        if self.save==True:
            
            isdir = os.path.isdir(self.folder_save+self.save_to+'\\') 
            if isdir==True:
                pass
            elif isdir==False:
                print('Create a folder with name : '+self.save_to+', in '+self.folder_save)
                sys.exit()
            
            np.save(self.folder_save+self.save_to+'\\'+'Mts('+str(self.Nmtmax)+')_position_loop_'+str(0)+'.npy',self.X_position[0]);  
            np.save(self.folder_save+self.save_to+'\\'+'Mts('+str(self.Nmtmax)+')_height_loop_'+str(0)+'.npy',self.LLmt[:,0]);  
        

        for j in range(1,self.Nsteps):
        
            pres=1-np.exp(-self.kres*self.dt); # probability of rescue in cell interior
        
            Deforming_MT_id=[]
            Deforming_len=[]
            Attracted_neighbours=[]
            Force_of_attraction=[]
            
            for i in range(self.Nmtmax): # going through microtubules
                
                self.Kc[:,j][i]=cat_freq_update(self.Vg[:,j-1][i],self.Tub[j-1]) # catastophic frequency
                self.pc_new=cat_prob_update(self.Lmt[i],self.Kc[:,j][i], self.Rcell, self.dt)
                self.Pc[:,j][i]=self.pc_new # probability of catastrophy
                
                if self.state[i]==0: #test for nucleation
                    if np.random.rand()<self.pnuc:
                        self.state[i]=1 #put into growing state if nucleated
                        self.Lmt[i]=self.Lmt[i]+self.Vg[:,j-1][i]*self.dt;
                elif self.state[i]==1: #for growing mts
                    self.Lmt[i]=self.Lmt[i]+self.Vg[:,j-1][i]*self.dt;
                    if np.random.rand()<self.Pc[:,j][i]:
                        self.state[i]=2 #put into shrinking state if nucleated
                                                              
                elif self.state[i]==2:
                    self.Lmt[i]=self.Lmt[i]-self.Ls
                    if np.random.rand()<pres:
                        self.state[i]=1;
                
                if self.Lmt[i]<0:
                    self.Lmt[i]=0;
                    self.state[i]=0; 
                    
                    force_on_mt = 0
                    self.Force_from_membrane_[:,j][i]=force_on_mt 
                    vg_new = velocity_g(force_on_mt,self.Tub[j-1]) 
                    self.Vg[:,j][i] = vg_new
                
                elif 0<self.Lmt[i]<=self.Rcell:    
                    force_on_mt = 0
                    self.Force_from_membrane_[:,j][i]=force_on_mt 
                    vg_new = velocity_g(force_on_mt,self.Tub[j-1]) 
                    self.Vg[:,j][i] = vg_new
                    
                elif self.Lmt[i]>self.Rcell: #if microtubule is protruding
                    min_radius = self.lf*(self.Lmt[i]-self.Rcell)

                    neighbours = Attracted_MT(i,min_radius,self.X_position[j-1],self.Nmtmax)
                    if len(neighbours)==0:
                        pass # otherwise empty neighbour list will get appended
                    else:    
                        Deforming_MT_id.append(i)
                        Deforming_len.append(self.Lmt[i]-self.Rcell)
                        Attracted_neighbours.append(neighbours)               
                        for k in range(len(neighbours)): 
                            force_of_attraction=[0,0,0]
                            force=attractive_force(self.X_position[j-1][i]-self.X_position[j-1][neighbours[k]],self.Lmt[i],self.Lmt[i]-self.Rcell,min_radius)
                            force_of_attraction[0]=i
                            force_of_attraction[1]=neighbours[k]
                            force_of_attraction[2]=force
                            Force_of_attraction.append(force_of_attraction)
                    force_on_mt = force_from_membrane(self.Lmt[i]-self.Rcell)
                    self.Force_from_membrane_[:,j][i]=force_on_mt 
                    vg_new = velocity_g(force_on_mt,self.Tub[j-1]) 
                    self.Vg[:,j][i] = vg_new 
                       
                self.LLmt[i,j]=self.Lmt[i];
                
            Attracted_neighbours = updating_neighbours(Attracted_neighbours, self.LLmt, self.Rcell, j)

            ## Order the Deforming_mt_id and attracted neighbours according to the demormation that it creates from low to high
            Deforming_MT_id, Attracted_neighbours = sorting_mts(Attracted_neighbours, Deforming_MT_id, Deforming_len)
                
            ## checking whether one microtubule is being attracted by multiple deforming mts at the same time
            if len(Attracted_neighbours)==0:
                pass # pass and goes out of the loop
            else:
                Attracted_neighbours_updated=[]
                for k in range(len(Attracted_neighbours)): 
                    if len(Attracted_neighbours[k])==0.0:
                        Attracted_neighbours_updated.append(Attracted_neighbours[k])
                        continue
                    else:
                        index_for_deleting = []
                        for q in range(len(Attracted_neighbours[k])):
        
                            l=int(Attracted_neighbours[k][q]) ## choose the mt
                            domi=finding_attracting_mt(l,Force_of_attraction, self.Nmtmax) ## find the mt that attracts l
                            if domi == self.Nmtmax*10:
                                index_for_deleting.append(q) # not attracted by any other
                            elif Deforming_MT_id[k]==domi:
                                pass
                            elif Deforming_MT_id[k]!=domi:
                                index_for_deleting.append(q) # not attracted by the Deforming_MT_id[k]    
                    # Attracted_neighbours[k][0] = np.delete(Attracted_neighbours[k],index_for_deleting) 
                    Attracted_neighbours_updated.append(np.delete(Attracted_neighbours[k],index_for_deleting))
                Attracted_neighbours=Attracted_neighbours_updated
            
            ## checking whether two deforming mts are attracting each other . If yes give preference to the one that attracts with strongest force
            if len(Attracted_neighbours)==0:
                pass # pass ang goes out of the loop
            else: 
                Index_for_deleting = []
                for c in range(len(Deforming_MT_id)):
                    temp = []
                    for d in range(len(Attracted_neighbours)):
                        if checking_for_element(Deforming_MT_id[c],Attracted_neighbours[d])==False:
                            temp.append(1)
                        else:                   
                            continue
                    if len(temp)==len(Attracted_neighbours):
                        Index_for_deleting.append(c)
                    else:
                        continue
                Deforming_MT_id_temp = np.delete(Deforming_MT_id,Index_for_deleting) ## new list of deforming mts that are being attracted by atleast by another deforming mts
 
                Pairwise_= Pairwise_sorting(Deforming_MT_id_temp) # non repeating pair the deforming mts
                if Pairwise_ == False:
                    pass
                else:
                    for e in range(len(Pairwise_)):
                        a1 = Pairwise_[e][0] # first mt
                        a2 = Pairwise_[e][1] # second mt
                        a1_index = Finding_index(a1,Deforming_MT_id) ## index of the first mt in Deforming_MT_id
                        a2_index = Finding_index(a2,Deforming_MT_id) ## index of the second mt in Deforming_MT_id
                        if checking_for_element(a1,Attracted_neighbours[a2_index])==True and checking_for_element(a2,Attracted_neighbours[a1_index])==True: ## if they are attracting each other
                            f12 = finding_force_of_attraction(a1_index,a2_index,Force_of_attraction,Deforming_MT_id) ## force by a1 on a2
                            f21 = finding_force_of_attraction(a2_index,a1_index, Force_of_attraction, Deforming_MT_id) ## force by a2 on a1
                            f_max= max(f12,f21)
                            if f_max ==0.0:
                                a1_index_in_neighbour_list = Finding_index(a1,Attracted_neighbours[a2_index])
                                Attracted_neighbours[a2_index]=np.delete(Attracted_neighbours[a2_index],a1_index_in_neighbour_list) # delete first mt from the attracting mt list of the second
                                a2_index_in_neighbour_list = Finding_index(a2,Attracted_neighbours[a1_index])
                                Attracted_neighbours[a1_index]=np.delete(Attracted_neighbours[a1_index],a2_index_in_neighbour_list)  # delete second mt from the attracting mt list of the first
                                    
                            elif f_max == f12: # first mt is the strongest
                                a1_index_in_neighbour_list = Finding_index(a1,Attracted_neighbours[a2_index])
                                Attracted_neighbours[a2_index]=np.delete(Attracted_neighbours[a2_index],a1_index_in_neighbour_list) # delete first mt from the attracting mt list of the second
                                Attracted_neighbours[a1_index]=np.unique(np.concatenate([Attracted_neighbours[a1_index],Attracted_neighbours[a2_index]])) 
                                Attracted_neighbours[a2_index]=   [] 
                            elif f_max == f21: # second mt is the strongest
                                a2_index_in_neighbour_list = Finding_index(a2,Attracted_neighbours[a1_index])
                                Attracted_neighbours[a1_index]=np.delete(Attracted_neighbours[a1_index],a2_index_in_neighbour_list)  # delete second mt from the attracting mt list of the first
                                Attracted_neighbours[a2_index]=np.unique(np.concatenate([Attracted_neighbours[a1_index],Attracted_neighbours[a2_index]]))
                                Attracted_neighbours[a1_index]=   [] 
                        
           #updating the spacing between the MTs
           
            ll = len(Attracted_neighbours)
            if ll==0:
                new_position = np.sort(self.X_position[j-1])
            else:
                new_position = self.X_position[j-1]
                for i in range(ll):
                    
                    if len(Attracted_neighbours[i])==0:
                        pass
                    else:
                        deforming_mt = Deforming_MT_id[i] 
                        new_position = MT_position_update(deforming_mt,new_position,Attracted_neighbours[i])
                new_position = np.sort(new_position)
            
            self.X_position.append(new_position)
            
            if self.save==True:            
                if j % 1 == 0:       
                    np.save(self.folder_save+self.save_to+'\\'+'Mts('+str(self.Nmtmax)+')_position_loop_'+str(j)+'.npy',self.X_position[j]);  
                    np.save(self.folder_save+self.save_to+'\\'+'Mts('+str(self.Nmtmax)+')_height_loop_'+str(j)+'.npy',self.LLmt[:,j]);  
        
            # update free tubulin concentration
            TotalMTlength=sum(self.Lmt); #µm
            self.TubulinFreeNumb=self.TubulinTotalNumb-1624*TotalMTlength;
            self.Tub[j]=self.TubulinFreeNumb/(602*self.Vcell); # free tubulin conc (µM)
            
            self.Deforming_MT_id_.append(Deforming_MT_id)
            self.Attracted_neighbours_.append(Attracted_neighbours)
            self.Force_of_attraction_.append(Force_of_attraction)

            print (j, end="\r")
            
        return self.X_position, self.Rcell, self.LLmt

if __name__=='__main__':
    
#    seed_num=np.random.randint(0,100000)
#    seed_num = 96152 # N=20
    seed_num=80532 # N=80
    np.random.seed(seed_num)
    print('random seed '+str(seed_num))
    
    
    ms = MtSimulation()
    X_position, Rcell, LLmt = ms.simulate()  

#    ms.plot_free_tubulin() # Free tubulin Vs time
#    ms.av_protrusion() # Average protrusion length of MTs Vs time
#    ms.av_velocity() # Average velocity of MTs Vs time

    t=20
    plot_bundle(t,X_position,LLmt,xrange=0,Rcell=Rcell,label=True) 

    t=ms.Nsteps-ms.dt
    plot_bundle(t,X_position,LLmt,xrange=0.1,Rcell=Rcell,label=True)
